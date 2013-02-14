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


#ifndef SPECTRALSTABILITYESTIMATOR_H
#define SPECTRALSTABILITYESTIMATOR_H

#include "AlignmentFile.h"
#include "FeatureFile.h"
#include "BatchFile.h"
#include "LexiconManager.h"
#include "PhoneSet.h"

namespace Bavieca {

// number of bins used to bin the 
#define NUMBER_BINS		50

// number of vectors used to compute the cepstral mean (excluding the central vector)
#define NUMBER_VECTORS_CONTEXT	6

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class SpectralStabilityEstimator {

	private:
	
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;

	public:
		
		// constructor
		SpectralStabilityEstimator(PhoneSet *phoneSet, LexiconManager *lexiconManager);
		
		// destructor
		~SpectralStabilityEstimator();
		
		// process a batch file
		bool processBatchFile(const char *strFileBatch, const char *strFileOutput);
		
		// compute the spectral stability of a feature frames for the given utterance
		bool computeSpectralStability(const char *strFileAlignment, const char *strFileFeatures);
		
		// compute the euclidean norm of a feature vector
		float computeEuclideanNorm(float *fFeatureVector) {
		
			/*float fAcc = 0.0;
			for(int i=0 ; i < CEPSTRAL_PARAMETERS_NUMBER ; ++i) {			// just the MFCC not the energy
				fAcc += fFeatureVector[i]*fFeatureVector[i];
			}
			fAcc = sqrt(fAcc);
			
			return fAcc;*/
			
			return 0.0f;
		}
		
		// hypothesizes regions where there is a filled-pause, it uses two thresholds
		// - a frame is considered to be likely a filled-pause frame if its instability falls below fMaxInstability
		// - a segment of N consecutive likely frames is hypothesized as a filled pause if (N >= iMinFrames)
		void hypothesizeFilledPauses(float *fInstability, int iFrames, float fMaxInstability, int iMinFrames);

};

};	// end-of-namespace

#endif
