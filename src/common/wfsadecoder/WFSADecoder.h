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


#ifndef WFSADECODER_H
#define WFSADECODER_H

#include "Global.h"

using namespace std;

#include <deque>
#include <list>
#include <string>
#include <vector>

#include "ActiveStateTable.h"

namespace Bavieca {

class BestPath;
class HMMManager;
class LexiconManager;
class PhoneSet;
class FeatureExtractor;

#define NUMBER_BINS   25          // number of bins created for histogram based pruning

#include "WFSAcceptor.h"

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class WFSADecoder {

	private:
	
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManager;
		LexiconManager *m_lexiconManager;
		LMManager *m_lmManager;
		WFSAcceptor *m_wfsAcceptor;
		HMMStateDecoding *m_hmmStatesDecoding;				// HMM-states
		unsigned int m_iHMMStatesDecoding;					// # HMM-states
		unsigned int m_iFeatureVectors;						// feature vectors
		
		float m_fScoreBest;										// best partial score
		
		ActiveStateTable *m_activeStateTable;				// table of active states
		
		// pruning
		float m_fPruningLikelihood;
		int m_iPruningMaxActiveStates;
		
		// lattice generation
		bool m_bLatticeGeneration;								// whether to generate a lattice
		int m_iMaxWordSequencesState;							// maximum number of word sequences arriving at any state

	public:

		// constructor
		WFSADecoder(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager, WFSAcceptor *wfsAcceptor, 
			int iPruningMaxActiveStates, float fPruningLikelihood, bool bLatticeGeneration, int iMaxWordSequencesState);

		// destructor
		~WFSADecoder();
		
		// initialization
		void initialize();
		
		// initialize the Viterbi search
		void initializeViterbi();	
		
		// perform Viterbi search
		void viterbi(Matrix<float> &mFeatures);
		
		// apply beam pruning to the set of active states
		void beamPruning(ActiveState *activeStates, unsigned int iActiveStates);
			
		// create a new active state
		inline ActiveState *newActiveState(StateX *state, float fScore, HMMStateDecoding *hmmStateDecoding) {
		
			ActiveState *activeState = new ActiveState;
			activeState->state = state;
			activeState->fScore = fScore;
			activeState->hmmStateDecoding = hmmStateDecoding;
		
			return activeState;
		}
		
		// return the best path
		BestPath *getBestPath();
		
		// return the hypothesis lattice
		HypothesisLattice *getHypothesisLattice();
		
};

};	// end-of-namespace

#endif
