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


#ifndef VITERBIX_H
#define VITERBIX_H

#include "HMMGraph.h"
#include "HMMManager.h"
#include "LexiconManager.h"

namespace Bavieca {

class HMMGraph;
class PhoneSet;

// trellis node for the Viterbi alignment
typedef struct {
	float fScore;					// log-likelihood
	double dViterbi;				// Viterbi log-likelihood
} VTrellisNode;

// return values for processUtterance
#define UTTERANCE_PROCESSED_SUCCESSFULLY													0
// error codes 
#define ERROR_CODE_EMPTY_TRANSCRIPTION														-1
#define ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS								-2
#define ERROR_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED				-3
#define ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY							-4
#define ERROR_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES						-5
#define ERROR_CODE_UNABLE_TO_CREATE_HMM_GRAPH											-6

// error codes in string format
#define ERROR_CODE_EMPTY_TRANSCRIPTION_STR											"no lexical units found in the transcription"
#define ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS_STR						"insufficient number of feature vectors to perform the alignment"
#define ERROR_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED_STR		"maximum trellis size exceeded: utterance too long"
#define ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY_STR					"not enough memory: utterance too long?"
#define ERROR_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES_STR				"numerical inaccuracies in fwd/bwd: utterance too long?"


/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ViterbiX {

	private: 
	
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		HMMManager *m_hmmManager;
		float m_fPruningBeam;
		
		int m_iFeatureDimensionality;
		int m_iTrellisMaxSizeBytes;
		
		// trellis cache (the idea is to keep the largest trellis alocated and reuse it across utterances)
		bool m_bTrellisCache;						// whether to cache the trellis
		int m_iTrellisCacheMaxSizeBytes;			// maximum size of a cached trellis
		VTrellisNode *m_trellisCache;			// cached trellis
		int m_iTrellisCacheSize;					// size of the cached trellis			

	public:

		// constructor
		ViterbiX(PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMManager *hmmManager, 
			float fPruningBeam, int iTrellisMaxSizeMB, bool bTrellisCache, int iTrellisCacheMaxSizeMB);

		// destructor
		~ViterbiX();
		
		// prints the alignment
		void print(VPhoneAlignment &vPhoneAlignment) {
		
			if (vPhoneAlignment.empty()) {
				cout << "<empty alignment>" << endl;
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
			unsigned int iPhoneCharactersMax = 0;
			for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
				unsigned int iLength = (unsigned int)strlen(m_phoneSet->getStrPhone(i));
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
				printf("%*s %.4f (%d) %-s\n",iPhoneCharactersMax,m_phoneSet->getStrPhone((*it)->iPhone),
					(*it)->fLikelihood,(*it)->iPosition,strLexUnit.c_str());
			}
		}		
		
		
		// process the given utterance
		// - handles multiple pronunciations
		// - handles optional symbols (typically silence+fillers)
		Alignment *processUtterance(VLexUnit &vLexUnitTranscription, bool bMultiplePronunciations, 
			VLexUnit &vLexUnitOptional, MatrixBase<float> &mFeatures, double *dUtteranceLikelihood, int &iErrorCode);
			
		// forward/backward computation on the trellis
		VTrellisNode *viterbi(MatrixBase<float> &mFeatures, int iNodes, FBNodeHMM **nodes, int iEdges, FBNodeHMM *nodeInitial, FBNodeHMM *nodeFinal, float fBeam, int *iErrorCode);
		
		// allocate memory for the trellis
		VTrellisNode *newTrellis(int iFeatures, int iEdges, int *iErrorCode);
		
		// delete a trellis
		void deleteTrellis(VTrellisNode *trellis);
		
		// print the trellis
		void print(VTrellisNode *node, int iRows, int iColumns);
		
};

};	// end-of-namespace

#endif
