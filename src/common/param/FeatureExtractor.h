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


#ifndef FEATUREEXTRACTOR_H
#define FEATUREEXTRACTOR_H

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "Global.h"
#include "Transform.h"

namespace Bavieca {

class ConfigurationFeatures;

// feature type
#define FEATURE_TYPE_MFCC		0	
#define FEATURE_TYPE_PLP		1

// feature extraction parameters
#define PREEMPHASIS_COEFFICIENT				0.97		// pre-emphasis coefficient
#define WINDOW_SIZE					20		// window size in milliseconds
#define SKIP_RATE					10		// skip rate in milliseconds

// window tapering method
#define WINDOW_TAPERING_METHOD_NONE			0		// no tapering
#define WINDOW_TAPERING_METHOD_HANN			1		// Hann window
#define WINDOW_TAPERING_METHOD_HAMMING			2		// Hamming window
#define WINDOW_TAPERING_METHOD_HANN_MODIFIED		3		// Hann window (modified)

// cepstral normalization mode
#define CEPSTRAL_NORMALIZATION_MODE_NONE		0		// unnormalized features
#define CEPSTRAL_NORMALIZATION_MODE_UTTERANCE		1		// training / transcription mode decoding
#define CEPSTRAL_NORMALIZATION_MODE_SESSION		2		// training / transcription mode decoding
#define CEPSTRAL_NORMALIZATION_MODE_STREAM		3		// live decoding

// string version
#define STR_CEPSTRAL_NORMALIZATION_MODE_NONE		"none"
#define STR_CEPSTRAL_NORMALIZATION_MODE_UTTERANCE	"utterance"
#define STR_CEPSTRAL_NORMALIZATION_MODE_SESSION		"session"
#define STR_CEPSTRAL_NORMALIZATION_MODE_STREAM		"stream"

// cepstral normalization method
#define CEPSTRAL_NORMALIZATION_METHOD_CMN		0		// cepstral mean normalization
#define CEPSTRAL_NORMALIZATION_METHOD_CMVN		1		// cepstral mean variance normalization

// string version
#define STR_CEPSTRAL_NORMALIZATION_METHOD_CMN		"CMN"
#define STR_CEPSTRAL_NORMALIZATION_METHOD_CMVN		"CMVN"

typedef vector<Matrix<float>* > VUtteranceFeatures; 

typedef struct {
	short *sSamples;
	unsigned int iSamples;
} SamplesUtterance;

typedef struct {
	SamplesUtterance samples;	
	Matrix<float> *mFeatures;	
} UtteranceData;

typedef vector<UtteranceData> VUtteranceData;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FeatureExtractor {

	private:
	
		int m_iSamplingRate;			// sampling rate
		int m_iSampleSize;			// sample size
		int m_bDCRemoval;			// whether to remove the DC offset (if any)
		int m_bPreemphasis;			// whether to preemphasize the waveform
		
		// windowing
		int m_iWindowWidth;			// width of the analysis window (in milliseconds)
		int m_iWindowShift;			// shift between adjacent analysis windows (in milliseconds)
		int m_iWindowTapering;			// window tapering mode (Hamming window, Hann window, etc)
		
		// feature type
		int m_iType;
		
		// filterbank
		int m_iFilterbankFrequencyMin;		// minimum frequency used to build the filterbank (useful for narrowband)
		int m_iFilterbankFrequencyMax;		// maximum frequency used to build the filterbank (useful for narrowband)
		int m_iFilterbankFilters;		// number of filters in the filterbank
		
		// DCT
		int m_iCepstralCoefficients;		// number of cepstral coefficients when applying the Discrete Cosine Transform
		
		bool m_bEnergy;				// whether to include normalized energy as a feature
		float m_fWarpFactor;			// warp factor for Vocal Tract Lenght Normalization
		
		// derivatives
		int m_iDerivativesOrder;		// derivatives order (2 for speed and acceleration, etc)
		int m_iDerivativesDelta;		// number of frames on each side of the regression window
		
		// spliced features
		int m_iSplicedSize;
		
		// cepstral normalization
		int m_iCepstralBufferSize;		// cepstral buffer size in frames
		int m_iCepstralNormalizationMode;	// cepstral normalization mode (utterance/stream)
		int m_iCepstralNormalizationMethod;	// cepstral normalization method
		bool m_bCepstralBufferFull;		// whether the cepstral buffer is full (circular buffer)
		
		int m_iCoefficients;			// number of coefficients (including cepstral coefficients and energy without derivatives)
		int m_iCoefficientsTotal;		// total number of coefficients in the whole feature vectors
		
			
		double *m_dCenterFrequencyMel;		// center frequencies of the filters in the Mel scale 		
		
		// FFT
		unsigned int m_iFFTPoints;			// number of points in the Fast Fourier Transform
		int *m_iFFTPointBin;					// keeps the lower bin of a given point in the FFT (each point is connected to two bins)
		double *m_dFFTPointGain;
		
		// cepstral buffer: it is used for different purposes
		// - cepstral normalization (mean/variance)
		// - cepstral context in live mode (last feature vectors from the previous speech chunk are used)	
		int m_iCepstralBufferPointer;		// pointer to the last used position in the cepstral buffer
		Matrix<float> *m_mCepstralBuffer;		// cepstral buffer
		
		// auxiliar variables
		unsigned int m_iSamplesFrame;
		unsigned int m_iSamplesSkip;
		
		// stream mode feature extraction (left context)
		int m_iSamplesStream;			// samples in current stream
		unsigned int m_iSamplesUsefulPrev;		// number of useful samples from previous chunk (size >= iSamplesFrame-iSamplesSkip)
		short *m_sSamplesUsefulPrev;		// useful samples from previous chunk of samples
		
		// remove DC-mean
		void removeDC(short *sSamples, int iSamples);
				
		// apply pre-emphasis
		void applyPreemphasis(short *sSamples, int iSamples);
		
		// apply Window-tapering
		void applyWindowTapering(double *dSamples, int iSamples);
		
		// apply the Hamming window
		void applyHammingWindow(double *dSamples, int iSamples);	
		
		// apply the Hann window (original)
		void applyHannWindowOriginal(double *dSamples, int iSamples);	
		
		// apply the Hann window (modified)
		void applyHannWindowModified(double *dSamples, int iSamples);	
		
		// convert a float value to a short value
		short convert(float f) {
		
			if (f > ((float)SHRT_MAX)) {
				return SHRT_MAX;
			} else if (f < ((float)SHRT_MIN)) {
				return SHRT_MIN;
			} else {
				return (short)f;
			}
		}
		
		// convert a vector from double to short
		void convert(double *dV, short *sV, int iLength) {
		
			int i;
			
			for(i = 0 ; i < iLength ; ++i) {
				if (dV[i] > ((double)SHRT_MAX)) {
					sV[i] = SHRT_MAX;
				} else if (dV[i] < ((double)SHRT_MIN)) {
					sV[i] = SHRT_MIN;
				} else {
					sV[i] = (short)dV[i];
				}
			}
		}		
		
		// computes the log energy of a speech frame
		double computeLogEnergy(double *dFrame, int iWindowLength);
		
		// convert Mel frequency to linear frequency
		double melToLinear(double dFrequency) {
		
			return 700.0*(exp(dFrequency/1127.0)-1.0);
		}

		// convert Mel frequency to linear frequency
		double linearToMel(double dFrequency) {
		
			return 1127.0*log((1.0)+(dFrequency/700.0));
		}
		
		// compute the lower bin connected to each FFT point
		void computeFFTPointBin();
		
		// extract MFCC features (only the static coefficients)
		Matrix<float> *extractStaticFeaturesMFCC(short *sSamples, unsigned int iSamples);
		
		// extract PLP features (only the static coefficients)
		Matrix<float> *extractStaticFeaturesPLP(short *sSamples, unsigned int iSamples);
		
		// compute derivatives
		Matrix<float> *computeDerivatives(MatrixBase<float> &mStatic);
		
		// stack features
		Matrix<float> *spliceFeatures(MatrixBase<float> &mFeatures, unsigned int iElements);	
		
		// compute CMN over a set of utterances (session mode)
		void applyCMN(VUtteranceFeatures &vUtteranceFeatures, unsigned int iCoeffNorm);	
		
		// compute CMVN over a set of utterances (session mode)
		void applyCMVN(VUtteranceFeatures &vUtteranceFeatures, unsigned int iCoeffNorm);	
		
		// aplly utterance-based cepstral mean normalization
		void applyCMN(Matrix<float> &mFeatures, unsigned int iCoeffNorm);
		
		// apply utterance-based cepstral mean variance normalization
		void applyCMVN(Matrix<float> &mFeatures, unsigned int iCoeffNorm);	
		
		// apply stream-based cepstral mean normalization
		void applyCMNStream(Matrix<float> &mFeatures, unsigned int iCoeffNorm);	
		
		// build the filter bank (it takes into account the warp factor)
		void buildFilterBank();	
		
		// return a warped frequency (VTLN) using a piece-wise linear function
		float getWarpedFrequency(float fFrequency);
		
		// return the feature type
		int getFeatureType(const char *strType) {
			
			// mfcc
			if (strcmp(strType,"mfcc") == 0) {
				return FEATURE_TYPE_MFCC;
			} 
			// plp
			else {
				assert(strcmp(strType,"plp") == 0);
				return FEATURE_TYPE_PLP;
			} 
		}
		
		// return the tapering window
		int getTaperingWindow(const char *strTapering) {
		
			if (strcmp(strTapering,"none") == 0) {
				return WINDOW_TAPERING_METHOD_NONE;
			} else if (strcmp(strTapering,"Hann") == 0) {
				return WINDOW_TAPERING_METHOD_HANN;
			} else if (strcmp(strTapering,"Hamming") == 0) {
				return WINDOW_TAPERING_METHOD_HAMMING;
			} else {
				assert(strcmp(strTapering,"HannModified") == 0);
				return WINDOW_TAPERING_METHOD_HANN_MODIFIED;
			}
		}

	public:
	
		// constructor
		FeatureExtractor(ConfigurationFeatures *configurationFeatures, float fWarpFactor, int iCepstralBufferSize, int iCepstralNormalizationMode, int iCepstralNormalizationMethod);
		
		// destructor
		~FeatureExtractor();
		
		// initialization
		void initialize();
		
		// extract features in batch mode, each entry in the path file is a pair [rawFile featureFile]
		bool extractFeaturesBatch(const char *strFileBatch, bool bHaltOnFailure = false);
		
		// extract features in  batch mode
		void extractFeaturesSession(VUtteranceData &vUtteranceData, bool bHaltOnFailure = false);
		
		// extract features from a RAW audio file and store them into the given file
		void extractFeatures(const char *strFileRaw, const char *strFileFea);
			
		// extract features (utterance-based cepstral normalization)
		Matrix<float> *extractFeatures(const char *strFile);
		
		// extract static features in stream mode (make use of "left context" speech samples)
		Matrix<float> *extractFeaturesStream(short *sSamples, unsigned int iSamples);
		
		// extract features (utterance or stream-based cepstral normalization)
		Matrix<float> *extractFeatures(short *sSamples, unsigned int iSamples);	
		
		// extract static features 
		Matrix<float> *extractStaticFeatures(short *sSamples, unsigned int iSamples);	
		
		// return the feature dimensionality
		unsigned int getFeatureDim() {
		
			return m_iCoefficientsTotal;
		}
		
		// return the dimensionality of a feature container (considers memory alignment)
		inline unsigned int getFeatureDimContainer() {
		#if defined __AVX__ || defined __SSE__
			// round up the dimensionality to the nearest multiple	
			unsigned int iReminder = m_iCoefficientsTotal%(ALIGN_BOUNDARY/sizeof(float));
			unsigned int iDimContainer = (iReminder == 0) ? m_iCoefficientsTotal : m_iCoefficientsTotal + (ALIGN_BOUNDARY/sizeof(float))-iReminder;
			assert(iDimContainer % ALIGN_BOUNDARY == 0);
			return iDimContainer;
		#else 
			return getFeatureDim();
		#endif	
		}
		
		// reset the stream for stream based feature extraction
		void resetStream() {
		
			m_iSamplesStream = 0;
			m_iSamplesUsefulPrev = 0;
		}
		
		// return the cepstral normalization mode
		static int getNormalizationMode(const char *strNormalizationMode) {
			
			if (strcmp(strNormalizationMode,STR_CEPSTRAL_NORMALIZATION_MODE_NONE) == 0) {
				return CEPSTRAL_NORMALIZATION_MODE_NONE;
			} else if (strcmp(strNormalizationMode,STR_CEPSTRAL_NORMALIZATION_MODE_UTTERANCE) == 0) {
				return CEPSTRAL_NORMALIZATION_MODE_UTTERANCE;
			} else if (strcmp(strNormalizationMode,STR_CEPSTRAL_NORMALIZATION_MODE_SESSION) == 0) {
				return CEPSTRAL_NORMALIZATION_MODE_SESSION;
			} else {
				assert(strcmp(strNormalizationMode,STR_CEPSTRAL_NORMALIZATION_MODE_STREAM) == 0);
				return CEPSTRAL_NORMALIZATION_MODE_STREAM;
			}
		}

		// return the cepstral normalization method
		static int getNormalizationMethod(const char *strNormalizationMethod) {
			
			if (strcmp(strNormalizationMethod,STR_CEPSTRAL_NORMALIZATION_METHOD_CMN) == 0) {
				return CEPSTRAL_NORMALIZATION_METHOD_CMN;
			} else {
				assert(strcmp(strNormalizationMethod,STR_CEPSTRAL_NORMALIZATION_METHOD_CMVN) == 0);
				return CEPSTRAL_NORMALIZATION_METHOD_CMVN;
			}
		}

};

#endif

}; // end-of-namespace
