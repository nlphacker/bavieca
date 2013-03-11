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


#ifndef FORWARDBACKWARD_H
#define FORWARDBACKWARD_H

#include "HMMManager.h"
#include "HypothesisLattice.h"
#include "LexiconManager.h"

namespace Bavieca {

class Alignment;
class PhoneSet;

// keeps information related to a HMM-state in the trellis used for training
typedef struct {
	float fScore;
	double dForward;				// forward log-likelihood
	double dBackward;				// backward log-likelihood
} NodeTrellis;

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

// maximum size (in MB) for the cache of likelihood computations
#define MAX_LIKELIHOOD_CACHE_SIZE_MB				500

// ad-hoc functions to use Duple as the key in a hash_map data structure
struct MTimeStateFunctions {

#if defined __linux__ || defined __APPLE__
	
	// comparison function (used for matching, comparison for equality)
	bool operator()(const pair<int,int> timeState1, const pair<int,int> timeState2) const {
	
		return ((timeState1.first == timeState2.first) && (timeState1.second == timeState2.second));
	}
	
#elif _WIN32
		
	static const size_t bucket_size = 4;
	static const size_t min_buckets = 8;

	// comparison function (used for ordering)
	bool operator()(const pair<int,int> timeState1, const pair<int,int> timeState2) const {
		
		if (timeState1.first == timeState2.first) {
			return (timeState1.second < timeState2.second);
		} else {
			return (timeState1.first < timeState2.first);
		}
	}

#endif
	
	// hash function
	size_t operator()(const pair<int,int> timeState) const {
	
		unsigned int iAcc = timeState.first;
		iAcc <<= 16;
		iAcc += timeState.second;
	
		return iAcc;
	}
};

#if defined __linux__ || defined __APPLE__
// maps (timeFrame+hmmId) to occupation probability (posterior prob)
typedef hash_map<pair<int,int> , double,MTimeStateFunctions,MTimeStateFunctions> MOccupation;
// maps (timeFrame+hmmId) to likelihood values
typedef hash_map<pair<int,int> , float,MTimeStateFunctions,MTimeStateFunctions> MLikelihood;
#elif _WIN32
// maps (timeFrame+hmmId) to occupation probability (posterior prob)
typedef hash_map<pair<int,int> , double,MTimeStateFunctions> MOccupation;
// maps (timeFrame+hmmId) to likelihood values
typedef hash_map<pair<int,int> , float,MTimeStateFunctions> MLikelihood;
#endif

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ForwardBackward {

	private:
	
		int m_iFeatureDimensionality;
		int m_iFeatureDimensionalityUpdate;
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManagerEstimation;				// models used to estimate occupation at the Gaussian level
		HMMManager *m_hmmManagerUpdate;					// models that are updated after the FB estimation (can be the same)
		float m_fForwardPruningBeam;
		float m_fBackwardPruningBeam;
		int m_iTrellisSizeMaxBytes;
		
		// trellis cache (the idea is to keep the largest trellis and reuse it across utterances)
		bool m_bTrellisCacheEnabled;				// whether a cache will be used in order to reuse allocated trellis
		int m_iTrellisCacheSizeMaxBytes;			// maximum cache size in bytes
		NodeTrellis *m_nodeTrellisCache;			// cached trellis of nodes 
		int m_iTrellisCacheSize;					// size of the cached trellis
		
		// likelihood cache
		MLikelihood *m_mLikelihoodCache;
		
		// counts the utterances processed
		int m_iUtterancesProcessed;
		int m_iLatticesProcessed;
		
		// allocate memory for the trellis
		NodeTrellis *newTrellis(int iHMMStates, int iFeatures);
		
		// delete a trellis
		void deleteTrellis(NodeTrellis *nodeTrellis);
		
		// given a trellis estimates forward-backward scores
		void computeFBTrellis(NodeTrellis *nodeTrellis, VHMMState &vHMMStateComposite, int iHMMStates, 
			float *fFeatureVectors, int iFeatureVectors, float fBackwardThreshold, int iOffset = 0);
		
		// print the given trellis (debugging purposes)
		void printTrellis(NodeTrellis *node, int iRows, int iColumns);	
		
		// compute observation prob
		float computeLikelihood(HMMState *hmmState, float *fFeatures, int iT);
		
	public:

		// constructor
		ForwardBackward(PhoneSet *phoneSet, HMMManager *hmmManagerEstimation, HMMManager *hmmManagerUpdate, float fForwardPruningBeam, float fBackwardPruningBeam, int iTrellisMaxSize, bool bTrellisCacheEnabled, int iTrellisCacheMaxSize);

		// destructor
		~ForwardBackward();
		
		// perform the forward-backward algorithm over an utterance
		Alignment *processUtterance(VLexUnit &vLexUnitTranscription, float *fFeaturesAlignment, 
			float *fFeaturesAccumulation, int iFeatures, double *dUtteranceLikelihood, const char **strReturnCode);
		
		// forward-backward alignment preserving phone-boundaries
		Alignment *processPhoneAlignment(float *fFeatures, int iFeatures, 
			VLPhoneAlignment *vLPhoneAlignment, double &dLikelihood, const char **strReturnCode);
		
		// aligns a lattice against the given feature vectors using Posterior Probabilities as the edge occupation
		// probability (discriminative training)
		MOccupation *processLattice(HypothesisLattice *lattice, float *fFeatures, int iFeatures, 
			float fScaleAM, float fScaleLM, double &dLikelihood, bool bMMI, float fBoostingFactor, const char **strReturnCode);
		
		// data conversion
		static Alignment *getAlignment(MOccupation *mOccupation, int iFrames);		
		
		// update the pruning beams
		void updatePruningBeams(float fForwardPruningBeam, float fBackwardPruningBeam) {
		
			m_fForwardPruningBeam = fForwardPruningBeam;
			m_fBackwardPruningBeam = fBackwardPruningBeam;
		}

};

};	// end-of-namespace

#endif
