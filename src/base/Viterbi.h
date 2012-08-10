/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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
#include "BestPath.h"
#include "HMMManager.h"
#include "HypothesisLattice.h"
#include "LexiconManager.h"
#include "PhoneSet.h"

// alignment error codes
#define ALIGNER_ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS			0
#define ALIGNER_ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY		1

// Viterbi beam search
#define VITERBI_BEAM_WIDTH		25*40.0

using namespace std;

#include <string>
#include <vector>
#include <deque>

#include <float.h>

#define INVALID	-DBL_MAX		// invalid path

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
		inline float computeEmissionProbability(HMMStateDecoding *hmmStateDecoding, float *fFeatureVector, int iFeatureVector);
		
		// return the number of HMM-states of a sequence of lexical units
		inline int getHMMStatesNumber(VLexUnit &vLexUnit) {
		
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
		
		// get the optimal alignment between a sequence of lexical units (typically words) and a speech segment by optionally:
		// - inserting filler models if needed
		// - select the phonetic pronunciation that best matches the audio
		VPhoneAlignment *align(VLexUnit &vLexUnit, float *fFeatures, int iFeatureVectors, float *fLikelihood);	
		
		// return a lexical unit alignment from a phone alignment
		static VLexUnitAlignment *getVLexUnitAlignment(VPhoneAlignment &vPhoneAlignment);	
		
		// return a given lexical unit alignment from a phone alignment
		LexUnitAlignment *getLexUnitAlignment(VPhoneAlignment &vPhoneAlignment, int iIndex);	
		
		// return the lexical units in the phonetic alignment
		void getLexUnitsFromPhoneAlignment(VPhoneAlignment *vPhoneAlignment, VLexUnit *vLexUnits);
		
		// return a state-level alignment given the BestPath
		VPhoneAlignment *align(float *fFeatures, int iFeatures, BestPath *bestPath);
		
		// align each of the lexical units in the lattice to the given set of feature vectors and store the
		// time alignment information into the edges 
		bool align(float *fFeatures, int iFeatures, HypothesisLattice *hypothesisLattice);	
		
		// prints the alignment
		inline void print(VPhoneAlignment &vPhoneAlignment) {
		
			if (vPhoneAlignment.empty() == true) {
				printf("<empty alignment>\n");
				return;
			}
		
			// determine the width of each alignment number
			int iDigits = 0;
			for(int i=10 ; i<INT_MAX ; i*=10) {
				++iDigits;
				if (vPhoneAlignment.back()->iStateEnd[NUMBER_HMM_STATES-1] < i) {
					break;
				}
			}
			
			// determine the width of the phone field
			int iPhoneCharactersMax = 0;
			for(int i=0 ; i < m_phoneSet->getSize() ; ++i) {
				int iLength = strlen(m_phoneSet->getStrPhone(i));
				if (iLength > iPhoneCharactersMax) {
					iPhoneCharactersMax = iLength;
				}
			}
			
			for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {	
				for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
					printf("%*d %*d ",iDigits,(*it)->iStateBegin[iState],iDigits,(*it)->iStateEnd[iState]);
				}
				string strLexUnit = "<undefined>";
				if (m_lexiconManager != NULL) {
					strLexUnit = m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit);
				} 
				printf("%*s %.4f (%d) %-s\n",iPhoneCharactersMax,m_phoneSet->getStrPhone((*it)->iPhone),(*it)->fLikelihood,(*it)->iPosition,strLexUnit.c_str());
			}
		
			return;
		}		
		
		// prints the lexical unit alignment
		inline void print(VLexUnitAlignment &vLexUnitAlignment) {
		
			float fLikelihoodTotal = 0.0;
		
			for(VLexUnitAlignment::iterator it = vLexUnitAlignment.begin() ; it != vLexUnitAlignment.end() ; ++it) {	
				printf("%5d %5d %-20s %10.2f\n",(*it)->iBegin,(*it)->iEnd,m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit),(*it)->fLikelihood);
				fLikelihoodTotal += (*it)->fLikelihood;
			}
			printf("total likelihood: %.2f\n",fLikelihoodTotal);
		
			return;
		}		
		
		// prints the lexical unit alignment
		inline void print(LLexUnitAlignment &lLexUnitAlignment) {
		
			float fLikelihoodTotal = 0.0;
		
			for(LLexUnitAlignment::iterator it = lLexUnitAlignment.begin() ; it != lLexUnitAlignment.end() ; ++it) {	
				printf("%5d %5d %-20s %10.2f\n",(*it)->iBegin,(*it)->iEnd,m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit),(*it)->fLikelihood);
				fLikelihoodTotal += (*it)->fLikelihood;
			}
			printf("total likelihood: %.2f\n",fLikelihoodTotal);
		
			return;
		}		
		
		// print the HMM-state composite
		inline void printHMMStateDecodingComposite(VHMMStateDecoding &vHMMStateDecodingComposite) {
		
			int i=0;
			for(VHMMStateDecoding::iterator it = vHMMStateDecodingComposite.begin() ; it != vHMMStateDecodingComposite.end() ; ++it, ++i) {
				printf("%4d: %6d\n",i,(*it)->getId());
			}
		
			return;
		}
		
		// activates or deactivates cache use
		inline void useCache(bool bUseCache) {
		 
			m_bUseCache = bUseCache; 
		}
		
		// allocate memory for the cache
		void allocateCache(int iFeatures);	
	
};

#endif
