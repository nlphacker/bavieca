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


#ifndef VITERBI_H
#define VITERBI_H

#include "AlignmentFile.h"
#include "HMMManager.h"
#include "LexiconManager.h"

using namespace std;

#include <string>
#include <vector>
#include <deque>

#include <float.h>

namespace Bavieca {

class BestPath;
class HypothesisLattice;
class PhoneSet;

// alignment error codes
#define ALIGNER_ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS			0
#define ALIGNER_ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY		1

// Viterbi beam search
#define VITERBI_BEAM_WIDTH		25*40.0

// Viterbi alignment node for the trellis used to implement the Viterbi search
typedef struct _VNode {
	double dScore;			// path score (we keep the best)
	int iHMMStatePrev;	// HMM-state associated to the previous node in the Viterbi Path (the time is t-1)	
} VNode;

/**
	@author root <root@localhost.localdomain>
*/
class Viterbi {

	private:
	
		int m_iFeatureDimensionality;
		PhoneSet *m_phoneSet;						// phonetic symbol set
		HMMManager *m_hmmManager;					// acoustic models manager		
		LexiconManager *m_lexiconManager;		// lexicon manager
		float m_fBeamWidth;							// beam-width used when doing Viterbi
		VNode *nodeTrellis;							// trellis to implement the Viterbi Search
		VNode *m_nodeTrellisCache;					// size of the trellis in the cache	
		int m_iTrellisCacheSize;					// to cache memory across alignments (if possible)
		float *m_fScoreCache;						// on position per frame and physical HMM-state [frame][HMM-state]
		int m_iHMMStatesPhysical;					// keeps the number of physical HMM-states (for cache access)
		bool m_bUseCache;								// whether the cache is used
		
		// cache recycling
		int m_iEntriesCacheOld;						// size of the old cache
		float *m_fScoreCacheOld;					// old cache
		
		// note: the use of a score cache is very helpful when we align difference sequences of lexical units against the
		// same segment of feature vectors, it reduces computation time considerably
		
		// return an HMM-composite from a sequence of lexical units
		void getHMMStateDecodingComposite(VLexUnit &vLexUnit, VHMMStateDecoding &vHMMStateDecodingComposite, LexUnit *lexUnitLeft = NULL, LexUnit *lexUnitRight = NULL);
		
		// align a sequence of HMM-states to the audio
		VPhoneAlignment *alignHMMStates(float *fFeatureVectors, int iFeatureVectors, VHMMStateDecoding &vHMMStateDecodingComposite, VLexUnit &vLexUnit, float *fLikelihood, int iOffset, unsigned char &iErrorCode);	
		
		// print a sequence of lexical units
		void printLexUnitSequence(VLexUnit &vLexUnitText);
		
		// compute emission probability or retrieves it from the cache
		float computeEmissionProbability(HMMStateDecoding *hmmStateDecoding, float *fFeatureVector, int iFeatureVector);
		
		// return the number of HMM-states of a sequence of lexical units
		int getHMMStatesNumber(VLexUnit &vLexUnit) {
		
			int iHMMStates = 0;
			for(VLexUnit::iterator it = vLexUnit.begin() ; it != vLexUnit.end() ; ++it) {
				iHMMStates += (*it)->vPhones.size();
			}
			
			return iHMMStates;
		}	

	public:

		// constructor
		Viterbi(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager, float fBeamWidth = VITERBI_BEAM_WIDTH);

		// destructor
		~Viterbi();
		
		// align a sequence of lexical units agains the features
		VPhoneAlignment *align(VLexUnit &vLexUnit, float *fFeatures, int iFeatureVectors, float *fLikelihood);	
		
		// return a state-level alignment given the BestPath
		VPhoneAlignment *align(float *fFeatures, int iFeatures, BestPath *bestPath);
		
		// align each of the lexical units in the lattice to the given set of feature vectors and store the
		// time alignment information into the edges 
		bool align(float *fFeatures, int iFeatures, HypothesisLattice *hypothesisLattice);	
		
		// print the HMM-state composite
		void printHMMStateDecodingComposite(VHMMStateDecoding &vHMMStateDecodingComposite) {
		
			int i=0;
			for(VHMMStateDecoding::iterator it = vHMMStateDecodingComposite.begin() ; it != vHMMStateDecodingComposite.end() ; ++it, ++i) {
				printf("%4d: %6d\n",i,(*it)->getId());
			}
		}
		
		// activates or deactivates cache use
		void useCache(bool bUseCache) {
		 
			m_bUseCache = bUseCache; 
		}
		
		// allocate memory for the cache
		void allocateCache(int iFeatures);		
};

};	// end-of-namespace

#endif
