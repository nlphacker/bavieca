/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
 *                                                                                             *
 * www.bavieca.org is the website of the Bavieca Speech Recognition Toolkit                    *
 *                                                                                             *
 * Licensed under the Apache License, Version 2.0 (the "License");                             *
 * you may not use this file except in compliance with the License.                            *
 * You may obtain a copy of the License at                                                     *
 *                                                                                             *
 *         http://www.apache.org/licenses/LICENSE-2.0                                          *
 *                                                                                             *
 * Unless required by applicable law or agreed to in writing, software                         *
 * distributed under the License is distributed on an "AS IS" BASIS,                           *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.                    *
 * See the License for the specific language governing permissions and                         *
 * limitations under the License.                                                              *
 *---------------------------------------------------------------------------------------------*/


#include "AudioFile.h"
#include "BatchFile.h"
#include "ExceptionBase.h"
#include "FeatureExtractor.h"
#include "FeatureFile.h"
#include "Numeric.h"
#include "Waveform.h"

#include "ConfigurationFeatures.h"

namespace Bavieca {

// note: the bank of filters used comprises N filters that are equally spaced on the mel scale
// all the filters have the same height and each filter bounday coincides with the center
// frequency of the adjacent filter

// constructor
FeatureExtractor::FeatureExtractor(ConfigurationFeatures *configurationFeatures, float fWarpFactor, int iCepstralBufferSize, int iCepstralNormalizationMode, int iCepstralNormalizationMethod) {

	m_iSamplingRate = configurationFeatures->getIntParameterValue("waveform.samplingRate");
	m_iSampleSize = configurationFeatures->getIntParameterValue("waveform.sampleSize");
	m_bDCRemoval = configurationFeatures->getBoolParameterValue("dcRemoval");
	m_bPreemphasis = configurationFeatures->getBoolParameterValue("preemphasis");
	m_iWindowWidth = configurationFeatures->getIntParameterValue("window.width");
	m_iWindowShift = configurationFeatures->getIntParameterValue("window.shift");
	m_iWindowTapering = getTaperingWindow(configurationFeatures->getStrParameterValue("window.tapering"));
	m_iType = getFeatureType(configurationFeatures->getStrParameterValue("features.type"));
	
	m_iFilterbankFrequencyMin = configurationFeatures->getIntParameterValue("filterbank.frequency.min");
	m_iFilterbankFrequencyMax = configurationFeatures->getIntParameterValue("filterbank.frequency.max");
	m_iFilterbankFilters = configurationFeatures->getIntParameterValue("filterbank.filters");
	m_iCepstralCoefficients = configurationFeatures->getIntParameterValue("cepstralCoefficients");
	m_bEnergy = configurationFeatures->getBoolParameterValue("energy");
	m_fWarpFactor = fWarpFactor; 
	m_iDerivativesOrder = configurationFeatures->getIntParameterValue("derivatives.order");
	m_iDerivativesDelta = configurationFeatures->getIntParameterValue("derivatives.delta");
	
	//m_iSplicedSize = configurationFeatures->getIntParameterValue("spliced.size");
	m_iSplicedSize = -1;
	
	m_iCepstralBufferSize = iCepstralBufferSize;
	m_iCepstralNormalizationMode = iCepstralNormalizationMode;
	m_iCepstralNormalizationMethod = iCepstralNormalizationMethod;	
	m_bCepstralBufferFull = false;
	
	if (m_bEnergy) {
		m_iCoefficients = m_iCepstralCoefficients+1;
	} else {
		m_iCoefficients = m_iCepstralCoefficients;
	}
	m_iCoefficientsTotal = m_iCoefficients*(m_iDerivativesOrder+1);
		
	// filterbanks
	m_fCenterFrequencyMel = NULL;
		
	// fft
	m_iFFTPoints = 0;
	m_iFFTPointBin = NULL;
	
	// cepstral buffer	
	m_iCepstralBufferPointer = 0;
	m_fCepstralBuffer = NULL;
}

// destructor
FeatureExtractor::~FeatureExtractor()
{
	if (m_fCenterFrequencyMel != NULL) {
		delete [] m_fCenterFrequencyMel;
	}
	if (m_iFFTPointBin != NULL) {
		delete [] m_iFFTPointBin;
	}
	if (m_fFFTPointGain != NULL) {
		delete [] m_fFFTPointGain;
	}
	if (m_fCepstralBuffer != NULL) {
		delete [] m_fCepstralBuffer;
	}
	if (m_sSamplesUsefulPrev) {
		delete [] m_sSamplesUsefulPrev;
	}
}

// initialization
void FeatureExtractor::initialize() {

	// check parameters
	if (m_iCepstralNormalizationMode != CEPSTRAL_NORMALIZATION_MODE_NONE) {
		assert((m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMN) ||
			(m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMVN));		
	}
	
	assert(m_iSamplingRate % 1000 == 0);
	m_iSamplesFrame = (m_iSamplingRate/1000)*m_iWindowWidth;
	m_iSamplesSkip = (m_iSamplingRate/1000)*m_iWindowShift;

	// compute the number of points in the FFT (smaller power of two above the number of samples in the window)
	m_iFFTPoints = 2;
	while(m_iFFTPoints < m_iSamplesFrame) {
		m_iFFTPoints *= 2; 
	}

	// build the filterbank
	buildFilterBank();
	
	// cepstral buffer
	if (m_iCepstralBufferSize != -1) {
		m_fCepstralBuffer = new float[m_iCepstralCoefficients*m_iCepstralBufferSize];
	} else {
		m_fCepstralBuffer = NULL;
	}
	
	// stream mode feature extraction (might not be used)
	m_sSamplesUsefulPrev = new short[m_iSamplesFrame];
	m_iSamplesUsefulPrev = 0;
	m_iSamplesStream = 0;		
}

// build the filter bank (it takes into account the warp factor)
void FeatureExtractor::buildFilterBank() {

	// (1) compute the center frequency of each filter in the filterbank
	
	// allocate memory to store the center frequencies
	m_fCenterFrequencyMel = new float[m_iFilterbankFilters+2];
	
	// min and max allowed frequencies (for band-limited feature extraction)
	float fFrequencyMinMel = linearToMel((float)m_iFilterbankFrequencyMin);
	float fFrequencyMaxMel = linearToMel((float)m_iFilterbankFrequencyMax);	
		
	// filterbank centers are equally spaced in the mel scale
	float fStepMel = (fFrequencyMaxMel-fFrequencyMinMel)/(m_iFilterbankFilters+1);	
	
	// note: an extra center frequency is computed for an imaginary n+1 filter, this center 
	// is necessary for the computations during binning
	for(int i=1 ; i <= m_iFilterbankFilters+1 ; ++i) {
		m_fCenterFrequencyMel[i] = fFrequencyMinMel+(i*fStepMel);
	}
	
	// warped features
	if (m_fWarpFactor != 1.0) {

		// convert back to linear frequencies
		for(int i=1 ; i <= m_iFilterbankFilters+1 ; ++i) {
			float fLinear = melToLinear(m_fCenterFrequencyMel[i]);
			float fWarped = getWarpedFrequency(fLinear);
			//assert((fWarped > ((float)m_iFilterbankFrequencyMin)) && (fWarped < ((float)m_iFilterbankFrequencyMax)));
			m_fCenterFrequencyMel[i] = linearToMel(fWarped);
		}
	}
	
	// (2) 
	computeFFTPointBin();
}

// return a warped frequency (VTLN) using a piece-wise linear function
float FeatureExtractor::getWarpedFrequency(float fFrequency) {
	
	// extended filter bank (towards high frequencies)
	double dFrequencyMax = (double)m_iFilterbankFrequencyMax;
	double dW0;
	if (m_fWarpFactor > 1.0) {
		dW0 = (7.0/(8.0*m_fWarpFactor))*dFrequencyMax;	
	} 
	// compressed filter bank (towards low frequencies)
	else {
		dW0 = (7.0/8.0)*dFrequencyMax;	
	}
	// apply the piecewise function
	if (fFrequency <= dW0) {
		return m_fWarpFactor*fFrequency;
	} else {
		return (float)((m_fWarpFactor*dW0)+((dFrequencyMax-(m_fWarpFactor*dW0))/(dFrequencyMax-dW0))*(fFrequency-dW0));
	}	
}

// extract features in batch mode, each entry in the path file is a pair [rawFile featureFile]
bool FeatureExtractor::extractFeaturesBatch(const char *strFileBatch, bool bHaltOnFailure) {

	// load the batch file
	BatchFile batchFile(strFileBatch,"raw|features");
	batchFile.load();
	
	if (batchFile.size() < 1) {
		return false;
	}
	
	VUtteranceData vUtteranceData;
	for(unsigned int iUtterance = 0 ; iUtterance < batchFile.size() ; ++iUtterance) {
	
		// load the raw audio		
		int iSamples = -1;
		short int *sSamples = NULL;
		try {
			sSamples = AudioFile::load(batchFile.getField(iUtterance,"raw"),&iSamples);
		} catch(ExceptionBase) {
			if (bHaltOnFailure) {
				return false;
			}	
		}
	
		UtteranceData utteranceData;
		utteranceData.samples.sSamples = sSamples;
		utteranceData.samples.iSamples = iSamples;
		utteranceData.features.fFeatures = NULL;
		utteranceData.features.iFeatures = -1;
		vUtteranceData.push_back(utteranceData);
	}
	
	extractFeaturesSession(vUtteranceData,bHaltOnFailure);
	
	for(unsigned int iUtterance = 0 ; iUtterance < batchFile.size() ; ++iUtterance) {
	
		// write the features to file
		const char *strFileFeatures = batchFile.getField(iUtterance,"features");
		FeatureFile featureFile(strFileFeatures,MODE_WRITE,FORMAT_FEATURES_FILE_DEFAULT,m_iCoefficientsTotal);
		try {
			featureFile.store(vUtteranceData[iUtterance].features.fFeatures,
			vUtteranceData[iUtterance].features.iFeatures);
		} catch(ExceptionBase) {	
			if (bHaltOnFailure) {
				return false;
			}	
		}
		delete [] vUtteranceData[iUtterance].samples.sSamples;
		delete [] vUtteranceData[iUtterance].features.fFeatures;
	}	

	return true;
}

// extract features in  batch mode
void FeatureExtractor::extractFeaturesSession(VUtteranceData &vUtteranceData, bool bHaltOnFailure) {

	// utterance-based normalization
	if (m_iCepstralNormalizationMode == CEPSTRAL_NORMALIZATION_MODE_UTTERANCE) {
	
		// compute the static coefficients
		int iUtterance = 0;
		for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it, ++iUtterance) {
		
			it->features.fFeatures = extractFeatures(it->samples.sSamples,it->samples.iSamples,&it->features.iFeatures);
			if (it->features.fFeatures == NULL) {
				if (bHaltOnFailure) {
					BVC_ERROR << "unable to extract features from utterance: " << iUtterance;
				} else {
					BVC_WARNING << "unable to extract features from utterance: " << iUtterance;
				}
			}
		}
	} 
	// session-based normalization
	else if (m_iCepstralNormalizationMode == CEPSTRAL_NORMALIZATION_MODE_SESSION) {	
		
		// compute the static coefficients
		int iUtterance = 0;
		for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it, ++iUtterance) {
			
			// extract static features
			it->features.fFeatures = extractStaticFeatures(it->samples.sSamples,it->samples.iSamples,&it->features.iFeatures);
			if (it->features.fFeatures == NULL) {
				it->features.fFeatures = NULL;
				it->features.iFeatures = 0;
				if (bHaltOnFailure) {
					BVC_ERROR << "unable to extract static features from utterance: " << iUtterance;
				} else {
					BVC_WARNING << "unable to extract static features from utterance: " << iUtterance;
				}
			}
		}
		
		// allocate memory for the static coefficients
		FeaturesUtterance *featuresUtterance = new FeaturesUtterance[vUtteranceData.size()];
		iUtterance = 0;
		for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it, ++iUtterance) {
			featuresUtterance[iUtterance].fFeatures = it->features.fFeatures;
			featuresUtterance[iUtterance].iFeatures = it->features.iFeatures;
		}	
	
		// CMN
		if (m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMN) {
			applyCMN(featuresUtterance,vUtteranceData.size(),m_iCoefficients,m_iCepstralCoefficients);
		}
		// CMVN
		else { 
			assert(m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMVN);
			applyCMVN(featuresUtterance,vUtteranceData.size(),m_iCoefficients,m_iCepstralCoefficients);
		}
		
		delete [] featuresUtterance;
		
		// extract derivatives and create the final features
		for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it) {
		
			if (it->features.fFeatures == NULL) {
				continue;
			}
			
			float *fFeatures = computeDerivatives(it->features.fFeatures,it->features.iFeatures);
			if (fFeatures == NULL) {
				if (bHaltOnFailure) {
					BVC_ERROR << "unable to compute derivatives";					
				} else {
					BVC_WARNING << "unable to compute derivatives";		
					it->features.iFeatures = 0;
				}
			}
			
			delete [] it->features.fFeatures;
			it->features.fFeatures = fFeatures;
		}
	}
}


// extract features from a RAW audio file and store them into the given file (utterance-level cepstral normalization)
void FeatureExtractor::extractFeatures(const char *strFileRaw, const char *strFileFea) {

	// extract the features
	int iFeatures = -1;
	float *fFeatures = extractFeatures(strFileRaw,&iFeatures);
	assert(fFeatures);
	
	// write the features to file
	FeatureFile featureFile(strFileFea,MODE_WRITE,FORMAT_FEATURES_FILE_DEFAULT,m_iCoefficientsTotal);
	featureFile.store(fFeatures,iFeatures);
	delete [] fFeatures;
}

// extract features from a RAW audio file (utterance-level cepstral normalization)
float *FeatureExtractor::extractFeatures(const char *strFile, int *iFeatures) {

	// load the audio
	int iSamples = -1;
	short *sSamples = AudioFile::load(strFile,&iSamples);
	
	// extract features
	float *fFeatures = extractFeatures(sSamples,iSamples,iFeatures);
	if (fFeatures == NULL) {
		delete [] sSamples;
		return NULL;
	}
	
	delete [] sSamples;
	
	return fFeatures;
}

// extract static features in stream mode (make use of "left context" speech samples)
float *FeatureExtractor::extractFeaturesStream(short *sSamples, int iSamples, int *iFeatures) {

	float *fFeatures = NULL; 
	//int iSamplesLeft = m_iSamplesFrame-m_iSamplesSkip;
	if (m_iSamplesStream > 0) {
		short *sSamplesAll = new short[iSamples+m_iSamplesUsefulPrev];	
		memcpy(sSamplesAll,m_sSamplesUsefulPrev,m_iSamplesUsefulPrev*sizeof(short));
		memcpy(sSamplesAll+m_iSamplesUsefulPrev,sSamples,iSamples*sizeof(short));
		fFeatures = extractFeatures(sSamplesAll,iSamples+m_iSamplesUsefulPrev,iFeatures);
		delete [] sSamplesAll;
	} else {
		fFeatures = extractFeatures(sSamples,iSamples,iFeatures);
	}	
	// update useful data from previous (current) chunk
	m_iSamplesUsefulPrev = (((m_iSamplesUsefulPrev+iSamples)-m_iSamplesFrame)%m_iSamplesSkip)+(m_iSamplesFrame-m_iSamplesSkip);
	assert(iSamples > m_iSamplesUsefulPrev);
	memcpy(m_sSamplesUsefulPrev,sSamples+(iSamples-m_iSamplesUsefulPrev),m_iSamplesUsefulPrev*sizeof(short));
	m_iSamplesStream += iSamples;
	
	return fFeatures;
}

// extract features (utterance or stream-based normalization)
float *FeatureExtractor::extractFeatures(short *sSamples, int iSamples, int *iFeatures) {

	// extract the static features (cepstral normalization takes place inside)
	float *fFeaturesStatic = extractStaticFeatures(sSamples,iSamples,iFeatures);
	if (fFeaturesStatic == NULL) {
		return NULL;
	}
	
	// compute derivatives
	float *fFeatures = computeDerivatives(fFeaturesStatic,*iFeatures);
	if (fFeatures == NULL) {
		delete [] fFeaturesStatic;
		return NULL;
	}
	delete [] fFeaturesStatic;

	return fFeatures;
}

// extract static features
float *FeatureExtractor::extractStaticFeatures(short *sSamples, int iSamples, int *iFeatures) {

	// mfcc
	if (m_iType == FEATURE_TYPE_MFCC) {
		return extractStaticFeaturesMFCC(sSamples,iSamples,iFeatures);
	} 
	// plp
	else if (m_iType == FEATURE_TYPE_PLP) {
		return extractStaticFeaturesPLP(sSamples,iSamples,iFeatures);
	} 
	// not supported
	else {
		return NULL;
	}
}

// extract MFCC features (only the static coefficients)
float *FeatureExtractor::extractStaticFeaturesMFCC(short *sSamplesOriginal, int iSamples, int *iFeatures) {

	// (1) make a copy of the samples so they do not get modified
	short *sSamples = new short[iSamples];
	memcpy(sSamples,sSamplesOriginal,iSamples*sizeof(short));

	// (2) remove the DC mean
	if (m_bDCRemoval) {
		removeDC(sSamples,iSamples);
	}
	
	// (3) pre-emphasis
	if (m_bPreemphasis) {
		applyPreemphasis(sSamples,iSamples);
	}
	
	// check that the number of frames is positive, otherwise there are not enough samples available
	if (iSamples < m_iSamplesFrame) {
		BVC_WARNING << "insufficient number of audio samples (available: " << iSamples 
			<< " required: " << m_iSamplesFrame << ")";
		return NULL;
	}
	// compute the number of frames (feature vectors) to extract
	int iFrames = 1+(iSamples-m_iSamplesFrame)/m_iSamplesSkip;
	// allocate memory for the frames
	float *fMFCC = new float[iFrames*m_iCoefficients];	
	*iFeatures = iFrames;	

	// allocate memory for the Fast Fourier Transform output (there has to be room for the real and the imaginary part)
	double *dFFT = new double[2*m_iFFTPoints+2];
	double *dFFTMagnitude = new double[m_iFFTPoints];
	
	// allocate memory for the bins
	double *dBins = new double[m_iFilterbankFilters+1];
	
	// auxiliar buffer
	double *dSamplesFrame = new double[m_iSamplesFrame];
	
	// process each window
	for(int i = 0 ; i < iFrames ; ++i) {
	
		for(int j=0 ; j < m_iSamplesFrame ; ++j) {
			dSamplesFrame[j] = (double)sSamples[(i*m_iSamplesSkip)+j];
		}
	
		// do window tapering if needed
		applyWindowTapering(dSamplesFrame,m_iSamplesFrame);
		
		// compute the frame energy
		double dFrameEnergy = computeLogEnergy(dSamplesFrame,m_iSamplesFrame);
		
		// apply the Fourier transform and get the frequencies
		for(int j = 0 ; j < m_iSamplesFrame ; ++j) {
			dFFT[j] = dSamplesFrame[j]; 
		}
		for(int j = m_iSamplesFrame ; j < m_iFFTPoints+2 ; ++j) {
			dFFT[j] = 0.0; 
		}	
		Numeric::fft(dFFT,m_iFFTPoints);
		
		// compute the magnitude from the real and imaginary parts
		int iIndex = 2;
		for(int j = 2 ; j < m_iFFTPoints/2 ; ++j, iIndex+=2) {
			double dReal = dFFT[iIndex];
			double dImaginary = dFFT[iIndex+1];
			dFFTMagnitude[j] = (float)sqrt(dReal*dReal+dImaginary*dImaginary);
		}
		
		// apply the bank of filters
		for(int j=0 ; j < m_iFilterbankFilters+1 ; ++j) {
			dBins[j] = 0.0;
		}
		for(int j=2 ; j < m_iFFTPoints/2 ; ++j) {
			int iBin = m_iFFTPointBin[j];
			float fValue = (float)(m_fFFTPointGain[j]*dFFTMagnitude[j]);
			if (iBin > 0) {
				dBins[iBin] += fValue;
			}
			if (iBin < m_iFilterbankFilters) {
				dBins[iBin+1] += dFFTMagnitude[j]-fValue;
			}	
		}
		
		// apply log compression
		for(int h=0 ; h < m_iFilterbankFilters+1 ; ++h) {
			if (dBins[h] <= 1.0) {
				dBins[h] = 0.0;
			} else {
				dBins[h] = log(dBins[h]);
			}
		}
		
		// Discrete Cosine Transform (DCT): compute cepstral features from the log filterbank amplitudes
		// note: we are ignoring c0
		float *fMFCCAux = fMFCC+i*m_iCoefficients;
		for(int j=1 ; j <= m_iCepstralCoefficients ; ++j) {
			double dAux = 0.0;
			for(int k=1 ; k < m_iFilterbankFilters+1 ; ++k) {	
				dAux += dBins[k]*cos((PI_NUMBER*j*(((float)k)-0.5))/((float)m_iFilterbankFilters));
			}
			dAux *= sqrt(2.0/((float)m_iFilterbankFilters));
			fMFCCAux[j-1] = (float)dAux;
		}
		fMFCCAux[m_iCepstralCoefficients] = (float)dFrameEnergy;
	}
	
	// apply cepstral normalization (optional)
	// utterance
	if (m_iCepstralNormalizationMode == CEPSTRAL_NORMALIZATION_MODE_UTTERANCE) {
		// CMN
		if (m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMN) {
			applyCMN(fMFCC,*iFeatures,m_iCoefficients,m_iCepstralCoefficients);
		} 
		// CMVN
		else {
			assert(m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMVN);
			applyCMVN(fMFCC,*iFeatures,m_iCoefficients,m_iCepstralCoefficients);
		}
	} 
	// stream
	else if (m_iCepstralNormalizationMode == CEPSTRAL_NORMALIZATION_MODE_STREAM) {
		// CMN
		if (m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMN) {
			applyCMNStream(fMFCC,*iFeatures,m_iCoefficients,m_iCepstralCoefficients);
		} 
		// CMVN
		else {
			assert(m_iCepstralNormalizationMethod == CEPSTRAL_NORMALIZATION_METHOD_CMVN);
			BVC_ERROR << "CMVN not supported for stream mode";
		}
	}	
	// none (session-based normalization will be performed later)
	else {
		assert((m_iCepstralNormalizationMode == CEPSTRAL_NORMALIZATION_MODE_SESSION) ||
			(m_iCepstralNormalizationMode == CEPSTRAL_NORMALIZATION_MODE_NONE));
	}
	
	float fMinFrameEnergy = 10000.0;
	float fMaxFrameEnergy = 0.0;
	
	// compute normalized log energy
	fMaxFrameEnergy = max(fMFCC[m_iCepstralCoefficients], fMaxFrameEnergy);
	fMinFrameEnergy = min(fMFCC[m_iCepstralCoefficients], fMinFrameEnergy);	
	for(int i=1 ; i<iFrames ; i++) {	
		fMaxFrameEnergy = max(fMFCC[i*m_iCoefficients+m_iCepstralCoefficients], fMaxFrameEnergy);
		fMinFrameEnergy = min(fMFCC[i*m_iCoefficients+m_iCepstralCoefficients], fMinFrameEnergy);	
	}		
	
   for(int i=0 ; i<iFrames ; i++){
      fMFCC[i*m_iCoefficients+m_iCepstralCoefficients] += 1.0f - fMaxFrameEnergy;
      if (fMFCC[i*m_iCoefficients+m_iCepstralCoefficients] < -5.0) {
      	fMFCC[i*m_iCoefficients+m_iCepstralCoefficients] = -5.0;
      }
   }	
	
	// clean-up
	delete [] dBins;
	delete [] dFFT;
	delete [] dFFTMagnitude;
	delete [] dSamplesFrame;
	delete [] sSamples;

	return fMFCC;
}

// extract PLP features (only the static coefficients)
float *FeatureExtractor::extractStaticFeaturesPLP(short *sSamples, int iSamples, int *iFeatures) {

	return NULL;
}


// remove the DC-mean
void FeatureExtractor::removeDC(short *sSamples, int iSamples) {

	double *dSamples = new double[iSamples];

	// compute the mean
	double dAcc = 0.0;
	for(int i=0 ; i < iSamples ; ++i) {
		dAcc += (double)sSamples[i];
	}
	dAcc /= ((double)iSamples);
	
	// substract the mean
	for(int i=0 ; i < iSamples ; ++i) {
		dSamples[i] = ((double)sSamples[i])-dAcc;
	}
	
	// convert back to short
	convert(dSamples,sSamples,iSamples);
	
	delete [] dSamples;
}

// apply pre-emphasis to all the samples (note that some samples at the end may be discarded)
void FeatureExtractor::applyPreemphasis(short *sSamples, int iSamples) {

	for(int i=iSamples-1 ; i > 0 ; --i) {
		sSamples[i] = convert(((float)sSamples[i])-((float)PREEMPHASIS_COEFFICIENT*((float)sSamples[i-1])));
	}
	sSamples[0] = convert(((float)sSamples[0])-((float)PREEMPHASIS_COEFFICIENT*((float)sSamples[0])));
}

// apply Window-tapering
void FeatureExtractor::applyWindowTapering(double *dSamples, int iSamples) {

	// Hamming window
	if (m_iWindowTapering == WINDOW_TAPERING_METHOD_HAMMING) {
		applyHammingWindow(dSamples,iSamples);
	}
	// Hann window (original)
	else if (m_iWindowTapering == WINDOW_TAPERING_METHOD_HANN) {
		applyHannWindowOriginal(dSamples,iSamples);
	}
	// Hann window (modified)
	else if (m_iWindowTapering == WINDOW_TAPERING_METHOD_HANN_MODIFIED) {
		applyHannWindowModified(dSamples,iSamples);
	} 
	// no window tapering
	else {
		assert(m_iWindowTapering == WINDOW_TAPERING_METHOD_NONE);
	}
}

// apply the Hamming window
void FeatureExtractor::applyHammingWindow(double *dSamples, int iSamples) {

	double dConst = (2.0*PI_NUMBER)/((double)(iSamples-1));

	for(int i = 0 ; i < iSamples ; ++i) {
		dSamples[i] = dSamples[i]*(0.54-0.46*cos(dConst*i));
	}
}

// apply the Hann window (original)
// note: this is the original Hann window (initial and final samples are lost)
void FeatureExtractor::applyHannWindowOriginal(double *dSamples, int iSamples) {

	double dConst = (2.0*PI_NUMBER)/((double)(iSamples-1));

	for(int i = 0 ; i < iSamples ; ++i) {
		dSamples[i] = dSamples[i]*(0.5-0.5*cos(dConst*(i)));
	}
}

// apply the Hann window (modified)
// note: this is a modified version of the Hann window so the initial and final samples are not lost
// conventional Hann window uses: dConst = 2*PI/(iSamples-1), which makes the first and last sample to be 0
void FeatureExtractor::applyHannWindowModified(double *dSamples, int iSamples) {

	double dConst = (2.0*PI_NUMBER)/((double)(iSamples+1));

	for(int i = 0 ; i < iSamples ; ++i) {
		dSamples[i] = dSamples[i]*(0.5-0.5*cos(dConst*(i+1)));
	}
}

// computes the log energy of a speech frame
double FeatureExtractor::computeLogEnergy(double *dFrame, int iSamples){

   double dEnergy = 0.0; 
   for(int i = 0 ; i < iSamples ; i++) {
		dEnergy += dFrame[i]*dFrame[i];
   }

   dEnergy /= ((double)iSamples);
   if (dEnergy > 1.0) {
      dEnergy = log10(dEnergy);
   } else {
      dEnergy = 0.0;
   }

   return dEnergy;
}

// compute the lower bin connected to each FFT point
void FeatureExtractor::computeFFTPointBin() {

	m_iFFTPointBin = new int[m_iFFTPoints];
	float fStep = ((float)m_iSamplingRate)/((float)m_iFFTPoints);
	int iBin = 1;
	for(int i=1 ; i <= m_iFFTPoints/2 ; ++i) {
		while((m_fCenterFrequencyMel[iBin] < linearToMel((i-1)*fStep)) && (iBin < (m_iFilterbankFilters+1)))  {
			++iBin;
		}
		m_iFFTPointBin[i] = iBin-1;
	}
	
	m_fFFTPointGain = new float[m_iFFTPoints];
	for(int i=1 ; i <= m_iFFTPoints/2 ; ++i) {
		int iBin = m_iFFTPointBin[i];
		if (iBin > 0) {
			m_fFFTPointGain[i] = (m_fCenterFrequencyMel[iBin+1]-linearToMel((i-1)*fStep))/(m_fCenterFrequencyMel[iBin+1]-m_fCenterFrequencyMel[iBin]);
		} else {
			m_fFFTPointGain[i] = (m_fCenterFrequencyMel[1]-linearToMel((i-1)*fStep))/m_fCenterFrequencyMel[1];
		}
	}	
}

// print the features (debugging)
void FeatureExtractor::print(float *fFeatures, int iFeatures) {

	for(int i=0 ; i < iFeatures ; ++i) {
		printf("%5d",i);
		for(int j=0 ; j < m_iCoefficientsTotal ; ++j) {
			printf(" %7.4f",fFeatures[i*(m_iCoefficients*3)+j]);
		}
		printf("\n");
	}	
}

// compute derivatives
float *FeatureExtractor::computeDerivatives(float *fFeatures, int iFeatures) {

	int iOrder = m_iDerivativesOrder;
	int iDelta = m_iDerivativesDelta;

	assert(iOrder >= 1);
	int iFeaturesMinimum = ((iDelta*2)+1);
	if (iFeatures < iFeaturesMinimum) {
		BVC_WARNING << "insufficient number of feature vectors to compute derivatives (available: " 
			<< iFeatures << " required: " << iFeaturesMinimum << ")";	
		return NULL;
	}

	// allocate memory for the features
	float *fFeaturesDerivatives = new float[iFeatures*m_iCoefficientsTotal];
	for(int i=0 ; i < iFeatures*m_iCoefficientsTotal ; ++i) {
		fFeaturesDerivatives[i] = FLT_MAX;
	}	
	
	// precompute the constant
	float fDiv = 0.0;
	for(int i=1 ; i <= iDelta ; ++i) {
		fDiv += i*i;
	}
	fDiv *= 2.0;
	fDiv *= 0.5;
	
	// regression window	
	float **fDerivativeWindow = new float*[iOrder*2];

	// copy the original features and set the derivative positions to zero
	for(int i = 0 ; i < iFeatures ; ++i) {
		for(int j = 0 ; j < m_iCoefficients ; ++j) {
			fFeaturesDerivatives[i*m_iCoefficientsTotal+j] = fFeatures[i*m_iCoefficients+j]; 
		}
		for(int j = m_iCoefficients ; j < m_iCoefficientsTotal ; ++j) {
			fFeaturesDerivatives[i*m_iCoefficientsTotal+j] = 0.0; 
		}	
	}
	
	// compute the derivatives: n order derivatives are computed after n-1 order derivatives
	for(int h = 1 ; h <= iOrder ; ++h) {		
		for(int i = 0 ; i < iFeatures ; ++i) {
			// build the window
			for(int j = iDelta ; j >= 1 ; --j) {
				// replicate the first vector if needed
				if (i-j < 0) {
					fDerivativeWindow[iDelta-j] = fFeaturesDerivatives+m_iCoefficients*(h-1);
				} else {
					fDerivativeWindow[iDelta-j] = fFeaturesDerivatives+((i-j)*m_iCoefficientsTotal)+m_iCoefficients*(h-1);
				}
				// replicate the last vector if needed
				if (i+j >= iFeatures) {
					fDerivativeWindow[iDelta+j-1] = fFeaturesDerivatives+((iFeatures-1)*m_iCoefficientsTotal)+m_iCoefficients*(h-1);
				} else {
					fDerivativeWindow[iDelta+j-1] = fFeaturesDerivatives+((i+j)*m_iCoefficientsTotal)+m_iCoefficients*(h-1);
				}
			}	
			// compute derivatives from vectors in the window
			float *fDerivative = fFeaturesDerivatives+(i*m_iCoefficientsTotal)+m_iCoefficients*h;
			for(int j = iDelta ; j >= 1 ; --j) {
				for(int k = 0 ; k < m_iCoefficients ; ++k) {
					fDerivative[k] += j*(fDerivativeWindow[iDelta+j-1][k]-fDerivativeWindow[iDelta-j][k]);
				}
			}	
			for(int k = 0 ; k < m_iCoefficients ; ++k) {
				fDerivative[k] /= fDiv;
			}
		}
	}
	
	delete [] fDerivativeWindow;
	
	return fFeaturesDerivatives;
}

// apply utterance-based cepstral mean normalization
void FeatureExtractor::applyCMN(float *fFeatures, int iFeatures, int iCoefficients, int iCoefficientsNormalization) {

	double *dObservation = new double[iCoefficientsNormalization];
	double *dMean = new double[iCoefficientsNormalization];
	
	// initialization
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dObservation[i] = 0.0;
	}
 
 	// compute the mean
	for(int i=0 ; i < iFeatures ; ++i) {
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
			dObservation[j] += fFeatures[i*iCoefficients+j];
		}
	} 
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dMean[i] = (dObservation[i]/((double)iFeatures));
	}
	
	// substract the mean
	for(int i=0 ; i < iFeatures ; ++i) {
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
			fFeatures[i*iCoefficients+j] -= (float)dMean[j];
		}
	}
 
 	delete [] dObservation;
 	delete [] dMean;
}

// apply utterance-based cepstral mean normalization
void FeatureExtractor::applyCMVN(float *fFeatures, int iFeatures, int iCoefficients, int iCoefficientsNormalization) {

	double *dObservation = new double[iCoefficientsNormalization];
	double *dObservationSquare = new double[iCoefficientsNormalization];
	double *dMean = new double[iCoefficientsNormalization];
	double *dStandardDeviation = new double[iCoefficientsNormalization];
	
	// initialization
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dObservation[i] = 0.0;
		dObservationSquare[i] = 0.0;
		dMean[i] = 0.0;
		dStandardDeviation[i] = 0.0;
	}
 
 	// accumulate statistics
	for(int i=0 ; i < iFeatures ; ++i) {
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
			dObservation[j] += fFeatures[i*iCoefficients+j];
			dObservationSquare[j] += fFeatures[i*iCoefficients+j]*fFeatures[i*iCoefficients+j];
		}
	} 
	// compute the mean and standard deviation
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dMean[i] = (dObservation[i]/((double)iFeatures));
		dStandardDeviation[i] = sqrt((dObservationSquare[i]/((double)iFeatures))-(dMean[i]*dMean[i]));
	}
	
	// substract the mean
	for(int i=0 ; i < iFeatures ; ++i) {
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
			fFeatures[i*iCoefficients+j] = (float)((fFeatures[i*iCoefficients+j]-dMean[j])/dStandardDeviation[j]);
		}
	}	
 
 	delete [] dObservation;
 	delete [] dObservationSquare;
 	delete [] dMean;
 	delete [] dStandardDeviation;
}

// apply stream-based cepstral mean normalization
void FeatureExtractor::applyCMNStream(float *fFeatures, int iFeatures, int iCoefficients, int iCoefficientsNormalization) {
	
	double *dAcc = new double[iCoefficientsNormalization];
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dAcc[i] = 0.0;
	}
	
	// (1) compute the mean from the new vectors and those already in the buffer	
	int iCepstralBufferElements = -1;
	if (m_bCepstralBufferFull) {
		iCepstralBufferElements = m_iCepstralBufferSize;
	} else {
		iCepstralBufferElements = m_iCepstralBufferPointer;
	}
	for(int i=0 ;  i < iCepstralBufferElements ; ++i) {
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
			dAcc[j] += m_fCepstralBuffer[i*iCoefficientsNormalization+j];
		}
	}	
	for(int i=0 ;  i < iFeatures ; ++i) {
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
			dAcc[j] += fFeatures[i*iCoefficients+j];
		}
	}	
	int iElements = iCepstralBufferElements+iFeatures;
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dAcc[i] /= ((double)iElements);
	}	
	
	// (2) update the cepstral buffer (unnormalized features)
	for(int i=0 ; i < iFeatures ; ++i) {	
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {	
			m_fCepstralBuffer[m_iCepstralBufferPointer*iCoefficientsNormalization+j] = fFeatures[i*iCoefficients+j];
		}
		// circular buffer		
		if (++m_iCepstralBufferPointer >= m_iCepstralBufferSize) {
			m_bCepstralBufferFull = true;
			m_iCepstralBufferPointer = 0;
		}		
	}
	
	// (3) substract the mean
	for(int i=0 ; i < iFeatures ; ++i) {	
		for(int j=0 ; j < iCoefficientsNormalization ; ++j) {	
			fFeatures[i*iCoefficients+j] -= (float)dAcc[j];
		}		
	}
		
	// clean-up
	delete [] dAcc;
}

// compute CMN over a set of utterances (session mode)
void FeatureExtractor::applyCMN(FeaturesUtterance *featuresUtterance, int iUtterances, int iCoefficients, int iCoefficientsNormalization) {

	assert(featuresUtterance != NULL);
	assert(iUtterances >= 1);

	long lFeaturesTotal = 0;
	double *dObservation = new double[iCoefficientsNormalization];
	double *dMean = new double[iCoefficientsNormalization];
	
	// initialize
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dObservation[i] = 0.0;
		dMean[i] = 0.0;
	}
	
	// accumulate statistics
	for(int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) { 	
		float *fFeatures = featuresUtterance[iUtterance].fFeatures;
		if (fFeatures == NULL) {
			continue;
		}
		// accumulate statistics
		for(int i=0 ; i < featuresUtterance[iUtterance].iFeatures ; ++i) {
			for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
				dObservation[j] += fFeatures[i*iCoefficients+j];
			}
		} 
		lFeaturesTotal += featuresUtterance[iUtterance].iFeatures;
	}
	
	// compute the mean and the covariance
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dMean[i] = dObservation[i]/((double)lFeaturesTotal);
	}	
	
	// substract the mean and divide by the covariance
	for(int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) {
		float *fFeatures = featuresUtterance[iUtterance].fFeatures;
		if (fFeatures == NULL) {
			continue;
		}
		// accumulate statistics
		for(int i=0 ; i < featuresUtterance[iUtterance].iFeatures ; ++i) {
			for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
				fFeatures[i*iCoefficients+j] = (float)(fFeatures[i*iCoefficients+j]-dMean[j]);
			}	
		}
	}	
 
 	delete [] dObservation;
 	delete [] dMean;
}

// compute CMVN over a set of utterances (session mode)
void FeatureExtractor::applyCMVN(FeaturesUtterance *featuresUtterance, int iUtterances, int iCoefficients, int iCoefficientsNormalization) {

	assert(featuresUtterance != NULL);
	assert(iUtterances >= 1);

	long lFeaturesTotal = 0;
	double *dObservation = new double[iCoefficientsNormalization];
	double *dObservationSquare = new double[iCoefficientsNormalization];
	double *dMean = new double[iCoefficientsNormalization];
	double *dStandardDeviation = new double[iCoefficientsNormalization];	
	
	// initialize
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dObservation[i] = 0.0;
		dObservationSquare[i] = 0.0;
		dMean[i] = 0.0;
		dStandardDeviation[i] = 0.0;
	}
	
	// accumulate statistics
	for(int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) { 	
		float *fFeatures = featuresUtterance[iUtterance].fFeatures;
		if (fFeatures == NULL) {
			continue;
		}	
		// accumulate statistics
		for(int i=0 ; i < featuresUtterance[iUtterance].iFeatures ; ++i) {
			for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
				dObservation[j] += fFeatures[i*iCoefficients+j];
				dObservationSquare[j] += fFeatures[i*iCoefficients+j]*fFeatures[i*iCoefficients+j];
			}
		} 
		lFeaturesTotal += featuresUtterance[iUtterance].iFeatures;
	}
	
	// compute the mean and the covariance
	for(int i=0 ; i < iCoefficientsNormalization ; ++i) {
		dMean[i] = dObservation[i]/((double)lFeaturesTotal);
		dStandardDeviation[i] = sqrt((dObservationSquare[i]/((double)lFeaturesTotal))-(dMean[i]*dMean[i]));
	}	
	
	// substract the mean and divide by the covariance
	for(int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) {
		float *fFeatures = featuresUtterance[iUtterance].fFeatures;
		if (fFeatures == NULL) {
			continue;
		}
		// accumulate statistics
		for(int i=0 ; i < featuresUtterance[iUtterance].iFeatures ; ++i) {
			for(int j=0 ; j < iCoefficientsNormalization ; ++j) {
				fFeatures[i*iCoefficients+j] = (float)((fFeatures[i*iCoefficients+j]-dMean[j])/dStandardDeviation[j]);
			}	
		}
	}	
 
 	delete [] dObservation;
 	delete [] dObservationSquare;
 	delete [] dMean;
 	delete [] dStandardDeviation;
}

// splice features (concatenates static coefficients)
float *FeatureExtractor::spliceFeatures(float *fFeatures, int iFeatures, int iElements) {

	// it must be an odd number
	assert(iElements%2 == 1);

	float *fFeaturesSpliced = new float[iFeatures*m_iCoefficients*iElements];	
	int iContextSize = (iElements-1)/2;
	int *iWindow = new int[iElements];
	
	for(int i=0 ; i < iFeatures ; ++i) {
		// build the window
		for(int j=1 ; j <= iContextSize ; ++j) {
			// left context (replicate first vector if necessary)
			if (i-j >= 0) {
				iWindow[iContextSize-j] = i-j;
			} else {
				iWindow[iContextSize-j] = 0;
			}
			// right context (replicate last vector if necessary)
			if (i+j < iFeatures) {
				iWindow[iContextSize+j] = i+j;
			} else {
				iWindow[iContextSize+j] = iFeatures-1;
			}	
		}
		iWindow[iContextSize] = i;
		// create the stacked vector
		float *fVector = fFeaturesSpliced+i*(m_iCoefficients*iElements);
		for(int j=0 ; j < iElements ; ++j) {
			for(int k=0 ; k < m_iCoefficients ; ++k) {
				fVector[j*m_iCoefficients+k] = fFeatures[iWindow[j]*m_iCoefficients+k];
			}
		}
	}	
	
	delete [] iWindow;

	return fFeaturesSpliced;
}

}; // end-of-namespace



