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


#include "SpectralStabilityEstimator.h"

namespace Bavieca {

// constructor
SpectralStabilityEstimator::SpectralStabilityEstimator(PhoneSet *phoneSet, LexiconManager *lexiconManager) {

	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
}

// destructor
SpectralStabilityEstimator::~SpectralStabilityEstimator() {

}

// process a batch file
bool SpectralStabilityEstimator::processBatchFile(const char *strFileBatch, const char *strFileOutput) {

	assert(0);
	// load the batch file
	/*BatchFile *batchFile = new BatchFile(strFileBatch,BATCH_FILE_TYPE_FEATURES_ALIGNMENT);
	if (batchFile->load() == false) {
		return false;
	}
	
	// process all the entries in the batch file
	for(int i=0 ; i<batchFile->size() ; ++i) {	
		if (computeSpectralStability(batchFile->getAlignmentFile(i),batchFile->getFeatureFile(i)) == false) {
			return false;
		}
	}
	
	
	delete batchFile;*/

	return true;
}

// compute the spectral stability of a feature frames for the given utterance
bool SpectralStabilityEstimator::computeSpectralStability(const char *strFileAlignment, const char *strFileFeatures) {

	/*AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,m_lexiconManager);	
	FeatureFile *featureFile = new FeatureFile(strFileFeatures,MODE_READ);
	if (featureFile->load() == false) {
		return false;
	}
	
	// get the feature vectors
	int iFeatureVectors = -1;
	float *fFeatureVectors = featureFile->getFeatureVectors(&iFeatureVectors);
	
	// get the phonetic alignment
	VPhoneAlignment *vPhoneAlignment = alignmentFile->load(strFileAlignment);
	if (vPhoneAlignment == NULL) {
		return false;
	}
	
	// # of vectors to compute the mean (excluding the central vector)
	int iN = NUMBER_VECTORS_CONTEXT;
	float fMean[CEPSTRAL_PARAMETERS_NUMBER];
	float fAux[CEPSTRAL_PARAMETERS_NUMBER];
	
	// stats
	int iFramesSpeech = 0;
	int iFramesFiller = 0;
	int iFramesSilence = 0;
	float fAccSpeech = 0.0;
	float fAccFiller = 0.0;
	float fAccSilence = 0.0;
	int iBinsSpeech[NUMBER_BINS];
	int iBinsFiller[NUMBER_BINS];
	int iBinsSilence[NUMBER_BINS];
	for(int i=0 ; i < NUMBER_BINS ; ++i) {
		iBinsSpeech[i] = 0;
		iBinsFiller[i] = 0;
		iBinsSilence[i] = 0;
	}
	float *fInstability = new float[iFeatureVectors];
	
	unsigned char iPhonePrev = m_phoneSet->getPhoneIndexSilence();
	for(VPhoneAlignment::iterator it = vPhoneAlignment->begin() ; it != vPhoneAlignment->end() ; ++it) {
		// get the next phone
		unsigned char iPhoneNext;
		VPhoneAlignment::iterator jt = it;
		++jt;
		if (jt == vPhoneAlignment->end()) {
			iPhoneNext = m_phoneSet->getPhoneIndexSilence();
		} else {
			iPhoneNext = (*jt)->iPhone;
		}
		// compute the cepstral instability of each of the feature vectors aligned to the phone
		for(int i=(*it)->iStateBegin[0] ; i<=(*it)->iStateEnd[NUMBER_HMM_STATES-1] ;  ++i) {
			// initialize to zero
			for(int j = 0 ; j < CEPSTRAL_PARAMETERS_NUMBER ; ++j) {
				fMean[j] = 0.0;
			}
			// compute the mean
			int iFrameStart = max(0,i-(iN/2));
			int iFrameEnd = min(i+(iN/2),iFeatureVectors-1);
			for(int j = iFrameStart ; j <= iFrameEnd ; ++j) {
				for(int k = 0 ; k < CEPSTRAL_PARAMETERS_NUMBER ; ++k) {
					fMean[k] += fFeatureVectors[j*FEATURE_DIMENSIONALITY_DEFAULT+k];
				}
			}
			for(int k = 0 ; k < CEPSTRAL_PARAMETERS_NUMBER ; ++k) {
				fMean[k] /= (iFrameEnd-iFrameStart+1);
			}
			// compute the instability
			fInstability[i] = 0.0;
			for(int j = iFrameStart ; j <= iFrameEnd ; ++j) {
				for(int k = 0 ; k < CEPSTRAL_PARAMETERS_NUMBER ; ++k) {
					fAux[k] = fFeatureVectors[j*FEATURE_DIMENSIONALITY_DEFAULT+k]-fMean[k];
				}
				fInstability[i] += computeEuclideanNorm(fAux); 
			}
			fInstability[i] /= computeEuclideanNorm(&fFeatureVectors[i*FEATURE_DIMENSIONALITY_DEFAULT]);
			printf("%5d instability: %12.4f (mean: %12.4f)\n",i,fInstability[i],fMean[0]);
			
			
			if ((*it)->lexUnit == m_lexiconManager->getLexUnitSilence()) {
				++iFramesSilence;
				fAccSilence += fInstability[i];			
				int iBin = fInstability[i]*5;
				if (iBin >= NUMBER_BINS) {
					iBin = NUMBER_BINS-1;
				}
				iBinsSilence[iBin]++; 				
			} else if (m_lexiconManager->isFiller((*it)->lexUnit->iLexUnit) == true) {
				++iFramesFiller;
				fAccFiller += fInstability[i];
				int iBin = fInstability[i]*5;
				if (iBin >= NUMBER_BINS) {
					iBin = NUMBER_BINS-1;
				}
				iBinsFiller[iBin]++; 
			} else {
				++iFramesSpeech;
				fAccSpeech += fInstability[i];
				int iBin = fInstability[i]*5;
				if (iBin >= NUMBER_BINS) {
					iBin = NUMBER_BINS-1;
				}
				iBinsSpeech[iBin]++; 
			}
		}	
		iPhonePrev = (*it)->iPhone;
	}
	
	float fAverageSpeech = fAccSpeech/((float)iFramesSpeech);
	float fAverageFiller = fAccFiller/((float)iFramesFiller);
	float fAverageSilence = fAccSilence/((float)iFramesSilence);
	
	printf("speech:  %6d %12.4f %12.4f\n",iFramesSpeech,fAccSpeech,fAverageSpeech);
	printf("filler:  %6d %12.4f %12.4f\n",iFramesFiller,fAccFiller,fAverageFiller);
	printf("silence: %6d %12.4f %12.4f\n",iFramesSilence,fAccSilence,fAverageSilence);
	
	printf("      speech filler silence\n");
	for(int i=0 ; i < NUMBER_BINS ; ++i) {
		printf("%3d %6d %6d %6d\n",i,iBinsSpeech[i],iBinsFiller[i],iBinsSilence[i]);
	}
	
	hypothesizeFilledPauses(fInstability,iFeatureVectors,1.0,16);
	
	// clean-up
	AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
	delete featureFile;
	delete alignmentFile;*/

	return true;
}

// hypothesizes regions where there is a filled-pause, it uses two thresholds
// - a frame is considered to be likely a filled-pause frame if its instability falls below fMaxInstability
// - a segment of N consecutive likely frames is hypothesized as a filled pause if (N >= iMinFrames)
void SpectralStabilityEstimator::hypothesizeFilledPauses(float *fInstability, int iFrames, float fMaxInstability, int iMinFrames) {
	
	int iFramesHypothesized = 0;
	int iSegmentsHypothesized = 0;
	int iSegmentStart = -1;
	int iLikelySegmentSize = 0;
	for(int i=0 ; i<iFrames ; ++i) {
		if (fInstability[i] <= fMaxInstability) {
			++iLikelySegmentSize;
			if (iSegmentStart == -1) {
				iSegmentStart = i;
			}
		} else {
			if (iLikelySegmentSize >= iMinFrames) {
				printf("%d <filler> %d\n",iSegmentStart,i-1);
				++iSegmentsHypothesized;
				iFramesHypothesized += i-iSegmentStart;
			}
			iLikelySegmentSize = 0;
			iSegmentStart = -1;
		}
	}
	float fFramesHypothesizedPercent = (100.0*iFramesHypothesized)/((float)iFrames);
	printf("%d segments were hypothesized (%5.2f%% of total frames)\n",iSegmentsHypothesized,fFramesHypothesizedPercent);
}

};	// end-of-namespace
