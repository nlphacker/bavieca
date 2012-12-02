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


#ifndef FORWARDBACKWARDX_H
#define FORWARDBACKWARDX_H

using namespace std;

#include "HMMGraph.h"
#include "HMMManager.h"
#include "LexiconManager.h"

namespace Bavieca {

class Alignment;
class HypothesisLattice;
class PhoneSet;

typedef struct {
	float fScore;					// log-likelihood
	double dForward;				// forward log-likelihood
	double dBackward;				// backward log-likelihood
} FBTrellisNode;

// return codes
#define FB_RETURN_CODE_SUCCESS																"success"
#define FB_RETURN_CODE_EMPTY_TRANSCRIPTION												"no lexical units found in the transcription"
#define FB_RETURN_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS							"insufficient number of feature vectors to perform the alignment"
#define FB_RETURN_CODE_UNABLE_TO_CREATE_HMM_GRAPH										"unable to create the graph of hmm-states"
#define FB_RETURN_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED			"maximum trellis size exceeded, utterance too long"
#define FB_RETURN_CODE_UTTERANCE_TOO_LONG_UNABLE_TO_CREATE_TRELLIS				"unable to allocate memory for the trellis, utterance too long"
#define FB_RETURN_CODE_UTTERANCE_UNKNOWN_LATTICE_HMMSTATE							"hmm-state in lattice was not found in hmm-set"
#define FB_RETURN_CODE_UNABLE_TO_GET_BEST_PATH_FROM_LATTICE							"unable to get the best path from the lattice"
#define FB_RETURN_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES					"numerical inaccuracies, utterance too long"

// default values for trellis and trellis cache size
#define MAX_TRELLIS_SIZE_MB_DEFAULT					500
#define MAX_TRELLIS_CACHE_SIZE_MB_DEFAULT			500

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ForwardBackwardX {

	private:
	
		PhoneSet *m_phoneSet;								// phone set
		LexiconManager *m_lexiconManager;				// lexicon manager
		HMMManager *m_hmmManagerAlignment;				// models used to do the forward/backward alignment
		HMMManager *m_hmmManagerAccumulation;			// models for which data is accumulated
		
		int m_iFeatureDimensionalityAlignment;
		int m_iFeatureDimensionalityAccumulation;
		
		// pruning 
		float m_fForwardPruningBeam;
		float m_fBackwardPruningBeam;
		
		// maximum trellis size (utterances that need bigger trellis size are discarded)
		int m_iTrellisMaxSizeBytes;
		
		// trellis cache (the idea is to keep the largest trellis alocated and reuse it across utterances)
		bool m_bTrellisCache;						// whether to cache the trellis
		int m_iTrellisCacheMaxSizeBytes;			// maximum size of a cached trellis
		FBTrellisNode *m_trellisCache;			// cached trellis
		int m_iTrellisCacheSize;					// size of the cached trellis		
		
		// forward/backward computation on the trellis
		FBTrellisNode *forwardBackward(int iFeatures, float *fFeatures, int iNodes, FBNodeHMM **nodes, 
			int iEdges, FBNodeHMM *nodeInitial, FBNodeHMM *nodeFinal, float fBeamForward, float fBeamBackward, const char **strReturnCode);
			
		// compute the number of unused positions in the trellis (used for analysis)
		void countUnusedPositions(FBTrellisNode *trellis, int iFeatures, int iNodes);	
			
		// print the trellis
		void print(FBTrellisNode *node, int iRows, int iColumns);
		
		// allocate memory for the trellis
		FBTrellisNode *newTrellis(int iFeatures, int iEdges, const char **strReturnCode);
		
		// delete a trellis
		void deleteTrellis(FBTrellisNode *trellis);
		
	public:

		// constructor
		ForwardBackwardX(PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMManager *hmmManagerAlignment,
			HMMManager *hmmManagerAccumulation, float fForwardPruningBeam, float fBackwardPruningBeam, 
			int iTrellisMaxSizeMB, bool bTrellisCache, int iTrellisCacheMaxSizeMB);

		// destructor
		~ForwardBackwardX();
		
		// process the given utterance
		Alignment *processUtterance(VLexUnit &vLexUnitTranscription, bool bMultiplePronunciations, 
			VLexUnit &vLexUnitOptional, float *fFeaturesAlignment, float *fFeaturesAccumulation, int iFeatures, 
			double *dUtteranceLikelihood, const char **strReturnCode);
			
		// aligns a lattice against the given feature vectors using Posterior Probabilities as the edge occupation
		// probability (discriminative training)
		Alignment *processLattice(HypothesisLattice *lattice, float *fFeatures, int iFeatures, const char **strReturnCode);
};

};	// end-of-namespace

#endif
