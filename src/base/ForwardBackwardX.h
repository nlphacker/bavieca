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

#ifndef FORWARDBACKWARDX_H
#define FORWARDBACKWARDX_H

using namespace std;

#include "Alignment.h"
#include "FeatureExtractor.h"
#include "HMMManager.h"
#include "HMMGraph.h"
#include "HypothesisLattice.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "Time.h"

typedef struct {
	float fScore;					// log-likelihood
	double dForward;				// forward log-likelihood
	double dBackward;				// backward log-likelihood
} FBTrellisNode;

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
		Log *m_log;
		
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
			int iEdges, FBNodeHMM *nodeInitial, FBNodeHMM *nodeFinal, float fBeamForward, float fBeamBackward, int *iErrorCode);
			
		// compute the number of unused positions in the trellis (used for analysis)
		void countUnusedPositions(FBTrellisNode *trellis, int iFeatures, int iNodes);	
			
		// print the trellis
		void print(FBTrellisNode *node, int iRows, int iColumns);
		
		// allocate memory for the trellis
		FBTrellisNode *newTrellis(int iFeatures, int iEdges, int *iErrorCode);
		
		// delete a trellis
		void deleteTrellis(FBTrellisNode *trellis);
		
	public:

		// constructor
		ForwardBackwardX(PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMManager *hmmManagerAlignment,
			HMMManager *hmmManagerAccumulation, float fForwardPruningBeam, float fBackwardPruningBeam, 
			int iTrellisMaxSizeMB, bool bTrellisCache, int iTrellisCacheMaxSizeMB, Log *log);

		// destructor
		~ForwardBackwardX();
		
		// process the given utterance
		Alignment *processUtterance(VLexUnit &vLexUnitTranscription, bool bMultiplePronunciations, 
			VLexUnit &vLexUnitOptional, float *fFeaturesAlignment, float *fFeaturesAccumulation, int iFeatures, 
			double *dUtteranceLikelihood, int &iErrorCode);
			
		// aligns a lattice against the given feature vectors using Posterior Probabilities as the edge occupation
		// probability (discriminative training)
		Alignment *processLattice(HypothesisLattice *lattice, float *fFeatures, int iFeatures, int &iErrorCode);
			
		// return the error message associated to the given error code
		inline const char *getErrorMessage(int iErrorCode) {
		
			switch(iErrorCode) {
				case ERROR_CODE_EMPTY_TRANSCRIPTION: {				
					return ERROR_CODE_EMPTY_TRANSCRIPTION_STR;
				}
				case ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS: {				
					return ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS_STR;
				}
				case ERROR_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED: {
					return ERROR_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED_STR;
				}
				case ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY: {
					return ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY_STR;
				}
				case ERROR_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES: {
					return ERROR_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES_STR;
				}	
				default: {
					assert(0);
				}
			}
			
			return NULL;
		}			

};

#endif
