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


#include "PhoneSet.h"
#include "TimeUtils.h"
#include "ViterbiX.h"

namespace Bavieca {

// constructor
ViterbiX::ViterbiX(PhoneSet *phoneSet, LexiconManager *lexiconManager, 
			HMMManager *hmmManager, float fPruningBeam, int iTrellisMaxSizeMB, bool bTrellisCache, 
			int iTrellisCacheMaxSizeMB)
{
	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
	m_hmmManager = hmmManager;
	m_iFeatureDimensionality = hmmManager->getFeatureDimensionality();	
	
	// pruning
	m_fPruningBeam = fPruningBeam;
	
	// max trellis size
	m_iTrellisMaxSizeBytes = iTrellisMaxSizeMB*1024*1024;
	
	// trellis cache
	m_bTrellisCache = bTrellisCache;
	m_iTrellisCacheMaxSizeBytes = iTrellisCacheMaxSizeMB*1024*1024;
	m_trellisCache = NULL;
	m_iTrellisCacheSize = 0;
}

// destructor
ViterbiX::~ViterbiX()
{
	if (m_trellisCache != NULL) {
		delete [] m_trellisCache;
	}
}

// process the given utterance
// - handles multiple pronunciations
// - handles optional symbols (typically silence+fillers)
Alignment *ViterbiX::processUtterance(VLexUnit &vLexUnitTranscription, bool bMultiplePronunciations, 
	VLexUnit &vLexUnitOptional, float *fFeatures, int iFeatures, double *dUtteranceLikelihood, int &iErrorCode) {
	
	//double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// make sure there is at least one lexical unit to align to
	if (vLexUnitTranscription.empty() == true) {
		iErrorCode = ERROR_CODE_EMPTY_TRANSCRIPTION;
		return NULL;	
	}	
	
	// create the HMM-graph
	int iNodes = -1;
	int iEdges = -1;
	FBNodeHMM *nodeInitial = NULL;
	FBNodeHMM *nodeFinal = NULL;
	HMMGraph *hmmGraph = new HMMGraph(m_phoneSet,m_lexiconManager,
		m_hmmManager,m_hmmManager,bMultiplePronunciations,vLexUnitOptional);
	FBNodeHMM **nodes = hmmGraph->create(vLexUnitTranscription,&iNodes,&iEdges,&nodeInitial,&nodeFinal);
	if (nodes == NULL) {
		delete hmmGraph;
		iErrorCode = ERROR_CODE_UNABLE_TO_CREATE_HMM_GRAPH;
		return NULL;
	}
	delete hmmGraph;

	assert(nodeInitial->iDistanceEnd % NUMBER_HMM_STATES == 0);
	
	// there can't be fewer feature vectors than HMM-states in the composite
	if (iFeatures < nodeInitial->iDistanceEnd) {
		HMMGraph::destroy(nodes,iNodes);
		iErrorCode = ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS;
		return NULL;
	}	
	
	// reset the emission probability computation (to avoid using cached computations that are outdated)
	VHMMStateDecoding vHMMState;
	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			assert(edge->hmmStateEstimation != NULL);
			vHMMState.push_back((HMMStateDecoding*)edge->hmmStateEstimation);
		}
	}
	m_hmmManager->resetHMMEmissionProbabilityComputation(vHMMState);		
	
	// do the actual Viterbi pass
	int iErrorCodeFB = -1;
	VTrellisNode *trellis = viterbi(iFeatures,fFeatures,iNodes,nodes,iEdges,
		nodeInitial,nodeFinal,m_fPruningBeam,&iErrorCodeFB);
	if (trellis == NULL) {
		HMMGraph::destroy(nodes,iNodes);
		iErrorCode = iErrorCodeFB;
		return NULL;
	}
	
	// get the maximum viterbi score from the terminal edges
	double dViterbiBest = -DBL_MAX;
	FBEdgeHMM *edgeBest = NULL; 
	for(FBEdgeHMM *edge = nodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		if (trellis[(iFeatures-1)*iEdges+edge->iEdge].dViterbi > dViterbiBest) {
			edgeBest = edge;
			dViterbiBest = trellis[(iFeatures-1)*iEdges+edge->iEdge].dViterbi;
		}
	}
	assert(edgeBest != NULL);
	
	// utterance likelihood
	assert(dViterbiBest != -DBL_MAX);
	assert(edgeBest != NULL);
	double dLikelihoodUtterance = dViterbiBest;
	
	//countUnusedPositions(trellis,iFeatures,iNodes);
	
	Alignment *alignment = new Alignment(ALIGNMENT_TYPE_VITERBI);	
	FBEdgeHMM *edgeTmp = edgeBest;
	int iState = 0;
	int iFrameEnd = -1;
	int iStatesLeft = 1;
	LexUnit *lexUnit = NULL;
	for(int t = iFeatures-1 ; t >= 0 ; --t) {
		FBEdgeHMM *edgePrevBest = NULL;
		double dViterbiPrevBest = -DBL_MAX;	
		// self-transition
		if (trellis[t*iEdges+edgeTmp->iEdge].dViterbi > dViterbiPrevBest) {
			dViterbiPrevBest = trellis[t*iEdges+edgeTmp->iEdge].dViterbi;
			edgePrevBest = edgeTmp; 
		}	
		// previous nodes
		for(FBEdgeHMM *edge = edgeTmp->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
			if (trellis[t*iEdges+edge->iEdge].dViterbi != -FLT_MAX) {
				if (trellis[t*iEdges+edge->iEdge].dViterbi > dViterbiPrevBest) {
					dViterbiPrevBest = trellis[t*iEdges+edge->iEdge].dViterbi;
					edgePrevBest = edge; 
				}
			}
		}	
		HMMStateDecoding *hmmStateDecoding = (HMMStateDecoding*)edgePrevBest->hmmStateUpdate;
		// state-level alignment
		VStateOcc *vStateOcc = new VStateOcc;	
		vStateOcc->push_back(Alignment::newStateOcc(hmmStateDecoding->getId(),1.0));
		alignment->addFrameAlignmentFront(vStateOcc);
		// word-level alignment
		if (hmmStateDecoding->getState() != iState) {
			--iStatesLeft;
		}
		iState = hmmStateDecoding->getState();
		if (iStatesLeft == 0) {
			iStatesLeft = edgePrevBest->lexUnit->vPhones.size()*NUMBER_HMM_STATES;
			if (lexUnit != NULL) {
				alignment->addLexUnitAlignmentFront(t+1,iFrameEnd,lexUnit);
			}
			iFrameEnd = t;
			lexUnit = edgePrevBest->lexUnit;
		}	
		edgeTmp = edgePrevBest;
	}	
	if ((iStatesLeft == 1) && (lexUnit != NULL)) {
		alignment->addLexUnitAlignmentFront(0,iFrameEnd,lexUnit);	
	}	
	
	// TODO if multiple pronunciations are allowed the alternatives at each edge might be wrong due to the
	// path recombination in HMMGraph, this needs to be addressed
	
	//double dTimeEnd = TimeUtils::getTimeMilliseconds();		
	//double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	
	//printf("alignment: %12.4f (%12.4f seconds)\n",dViterbiBest,dTimeSeconds);
	
	// clean-up
	deleteTrellis(trellis);
	HMMGraph::destroy(nodes,iNodes);	
	
	*dUtteranceLikelihood = dLikelihoodUtterance;

	iErrorCode = UTTERANCE_PROCESSED_SUCCESSFULLY;
	
	return alignment;
}


// forward/backward computation on the trellis
VTrellisNode *ViterbiX::viterbi(int iFeatures, float *fFeatures, int iNodes, FBNodeHMM **nodes, int iEdges, FBNodeHMM *nodeInitial, FBNodeHMM *nodeFinal, float fBeam, int *iErrorCode) {

	//double dTimeBegin = TimeUtils::getTimeMilliseconds();

	// allocate memory for the trellis
	VTrellisNode *trellis = newTrellis(iFeatures,iEdges,iErrorCode);
	if (trellis == NULL) {
		return NULL;
	}
	
	for(int i=0 ; i < iFeatures*iEdges ; ++i) {
		trellis[i].fScore = -FLT_MAX;
		trellis[i].dViterbi = -FLT_MAX;
	}

	// active states
	FBEdgeHMM **edgeActiveCurrent = new FBEdgeHMM*[iEdges];
	FBEdgeHMM **edgeActiveNext = new FBEdgeHMM*[iEdges];
	int iActiveCurrent = 0;
	int iActiveNext = 0;
	
	// mark all the edges as inactive except the initial edge
	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			edge->iActive = -1;
			//printf("%2d (%2d) %2d [%6d]\n",edge->nodePrev->iNode,edge->iEdge,edge->nodeNext->iNode,edge->hmmState->getId());
		}
	}
	
	// compute forward score of initial edges and activate their successors
	for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		trellis[edge->iEdge].fScore = ((HMMStateDecoding*)edge->hmmStateEstimation)->computeEmissionProbabilityNearestNeighborPDE(fFeatures,0);
		trellis[edge->iEdge].dViterbi = trellis[edge->iEdge].fScore;
		for(FBEdgeHMM *edge1 = edge->nodeNext->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
			edgeActiveCurrent[iActiveCurrent++] = edge1;
			edge1->iActive = 1;
		}
		// self loop
		edgeActiveCurrent[iActiveCurrent++] = edge;
		edge->iActive = 1;
	}
	
	// fill the trellis	
	for(int t = 1 ; t < iFeatures ; ++t) {
		float *fFeatureVector = fFeatures+(t*m_iFeatureDimensionality);
		// (a) compute Viterbi scores
		for(int i = 0 ; i < iActiveCurrent ; ++i) {		
			int s = edgeActiveCurrent[i]->iEdge;
			VTrellisNode *nodeTrellis = &trellis[t*iEdges+s];
			// accumulate probability mass from previous time-frame
			nodeTrellis->dViterbi = -FLT_MAX;
			// (1) same node	
			if (t > edgeActiveCurrent[i]->nodePrev->iDistanceStart) {
				if (trellis[(t-1)*iEdges+s].dViterbi != -FLT_MAX) {
					nodeTrellis->dViterbi = trellis[(t-1)*iEdges+s].dViterbi;	
				}
			}
			// (2) previous nodes
			for(FBEdgeHMM *edge = edgeActiveCurrent[i]->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
				if (t > edge->nodePrev->iDistanceStart) {
					//if (trellis[(t-1)*iEdges+edge->iEdge].dViterbi != -FLT_MAX) {
						//nodeTrellis->dViterbi = Numeric::logAddition(trellis[(t-1)*iEdges+edge->iEdge].dViterbi,nodeTrellis->dViterbi);
					//}
					if (trellis[(t-1)*iEdges+edge->iEdge].dViterbi != -FLT_MAX) {
						nodeTrellis->dViterbi = max(nodeTrellis->dViterbi,trellis[(t-1)*iEdges+edge->iEdge].dViterbi);
					}
				}
			}
			if (nodeTrellis->dViterbi != -FLT_MAX) {	
				HMMStateDecoding *hmmStateDecoding = (HMMStateDecoding*)(edgeActiveCurrent[i]->hmmStateEstimation);
				assert(hmmStateDecoding->getGaussianComponents() > 0);
				nodeTrellis->fScore = hmmStateDecoding->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,t);
				nodeTrellis->dViterbi += nodeTrellis->fScore;
			}
		}
		// (b) beam pruning
		// - get the best scoring node
		double dBestScore = -DBL_MAX;
		for(int i = 0 ; i < iActiveCurrent ; ++i) {
			if (trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dViterbi > dBestScore) {
				dBestScore = trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dViterbi;
			}
		}	
		// - compute the threshold (minimum backward value for a node not to be pruned)
		double dThreshold = dBestScore-m_fPruningBeam;
		assert(dThreshold > -DBL_MAX);	
		// (c) determine active states for next time frame
		iActiveNext = 0;
		for(int i=0 ; i<iActiveCurrent ; ++i) {
			// actual pruning	
			if (trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dViterbi < dThreshold) {
				continue;
			}		
			// self-loop
			if ((iFeatures-2-t) >= edgeActiveCurrent[i]->nodeNext->iDistanceEnd) {
				if (edgeActiveCurrent[i]->iActive < t+1) {
					edgeActiveCurrent[i]->iActive = t+1;
					edgeActiveNext[iActiveNext++] = edgeActiveCurrent[i];
				}
			}
			// successors
			for(FBEdgeHMM *edge = edgeActiveCurrent[i]->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
				if ((iFeatures-2-t) >= edge->nodeNext->iDistanceEnd) {
					if (edge->iActive < t+1) {
						edge->iActive = t+1;
						edgeActiveNext[iActiveNext++] = edge;
					}
				}
			}	
		}
		// swap
		iActiveCurrent = iActiveNext;
		FBEdgeHMM **edgeAux = edgeActiveNext;
		edgeActiveNext = edgeActiveCurrent;
		edgeActiveCurrent = edgeAux;
	}

	
	//double dTimeEnd = TimeUtils::getTimeMilliseconds();
	//double dTimeTotal = (dTimeEnd-dTimeBegin)/1000.0;
	//printf("%8.4f %8.4f total: %8.4f \n",dTimeBwd,dTimeFwd,dTimeTotal);
	
	//print(trellis,iFeatures,iEdges);
	
	// clean-up
	delete [] edgeActiveCurrent;
	delete [] edgeActiveNext;	

	return trellis;
}


// allocate memory for the trellis
VTrellisNode *ViterbiX::newTrellis(int iFeatures, int iEdges, int *iErrorCode) {

	// trellis size: # elements
	int iTrellisSize = iFeatures*iEdges;
	
	// check if the memory to be allocated does not exceed the maximum allowed
	int iTrellisSizeBytes = iTrellisSize*sizeof(VTrellisNode);
	if (iTrellisSizeBytes > m_iTrellisMaxSizeBytes) {	
		*iErrorCode = ERROR_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED;
		return NULL;
	}	
	
	// allocate memory for the trellis (if enabled, try to reuse a cached trellis)
	VTrellisNode *trellis = NULL;
	if ((m_bTrellisCache == false) || (m_trellisCache == NULL) || (m_iTrellisCacheSize < (iFeatures*iEdges))) {
		// try to allocate memory for the trellis
		try {
			trellis = new VTrellisNode[iFeatures*iEdges];
		}
		catch (const std::bad_alloc&) {
			*iErrorCode = ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY;
			return NULL;
		}
		// if caching is enabled and the new trellis it not too large, cache it 
		// note: not caching large trellis prevents a more efficient use of the main memory
		if ((m_bTrellisCache == true) && (iTrellisSizeBytes <= m_iTrellisCacheMaxSizeBytes)) {	
			// empty the cache if necessary
			if (m_trellisCache != NULL) {
				delete [] m_trellisCache;
			}
			// cache the new trellis
			m_trellisCache = trellis;
			m_iTrellisCacheSize = iFeatures*iEdges;
		}
	} else {
		trellis = m_trellisCache;
	} 
	
	return trellis;
}

// delete a trellis
void ViterbiX::deleteTrellis(VTrellisNode *trellis) {

	// clean-up
	if (m_trellisCache != trellis) {
		delete [] trellis;
	}	
}


// print the trellis
void ViterbiX::print(VTrellisNode *node, int iRows, int iColumns) {

	printf("\n");
	for(int i=0 ; i < iRows ; ++i) {
		for(int j=0 ; j < iColumns ; ++j) {
			float fScore = node[i*iColumns+j].fScore;
			if (fScore == -FLT_MAX) {
				fScore = -100000.0000;
			}
			double dViterbi = node[i*iColumns+j].dViterbi;
			if (dViterbi == -FLT_MAX) {
				dViterbi = -100000.0000;
			}
			printf("%12.4f %12.4f   ",fScore,dViterbi);
		}
		printf("\n");
	}
}

};	// end-of-namespace

