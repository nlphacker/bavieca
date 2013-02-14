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


#include "Alignment.h"
#include "ForwardBackwardX.h"
#include "HypothesisLattice.h"
#include "Numeric.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor 
ForwardBackwardX::ForwardBackwardX(PhoneSet *phoneSet, LexiconManager *lexiconManager, 
	HMMManager *hmmManagerAlignment, HMMManager *hmmManagerAccumulation, float fForwardPruningBeam, 
	float fBackwardPruningBeam, int iTrellisMaxSizeMB, bool bTrellisCache, 
	int iTrellisCacheMaxSizeMB)
{
	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
	m_hmmManagerAlignment = hmmManagerAlignment;
	m_hmmManagerAccumulation = hmmManagerAccumulation;
	m_iFeatureDimensionalityAlignment = m_hmmManagerAlignment->getFeatureDimensionality();	
	m_iFeatureDimensionalityAccumulation = m_hmmManagerAccumulation->getFeatureDimensionality();	
	
	// pruning
	m_fForwardPruningBeam = fForwardPruningBeam;
	m_fBackwardPruningBeam = fBackwardPruningBeam;
	
	// max trellis size
	m_iTrellisMaxSizeBytes = iTrellisMaxSizeMB*1024*1024;
	
	// trellis cache
	m_bTrellisCache = bTrellisCache;
	m_iTrellisCacheMaxSizeBytes = iTrellisCacheMaxSizeMB*1024*1024;
	m_trellisCache = NULL;
	m_iTrellisCacheSize = 0;
}

// destructor
ForwardBackwardX::~ForwardBackwardX()
{
	if (m_trellisCache != NULL) {
		delete [] m_trellisCache;
	}
}

// process the given utterance
// - handles multiple pronunciations
// - handles optional symbols (typically silence+fillers)
Alignment *ForwardBackwardX::processUtterance(VLexUnit &vLexUnitTranscription, bool bMultiplePronunciations, 
	VLexUnit &vLexUnitOptional, float *fFeaturesAlignment, float *fFeaturesAccumulation, int iFeatures, 
	double *dUtteranceLikelihood, const char **strReturnCode) {
	
	//double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// estimation properties
	bool bSingleGaussian = m_hmmManagerAlignment->isSingleGaussian();
	bool bAccumulatorsLogical = m_hmmManagerAccumulation->areAccumulatorsLogical();
	
	// make sure there is at least one lexical unit to align to
	if (vLexUnitTranscription.empty()) {
		*strReturnCode = FB_RETURN_CODE_SUCCESS;
		return NULL;	
	}	
	
	// create the HMM-graph
	int iNodes = -1;
	int iEdges = -1;
	FBNodeHMM *nodeInitial = NULL;
	FBNodeHMM *nodeFinal = NULL;
	HMMGraph *hmmGraph = new HMMGraph(m_phoneSet,m_lexiconManager,m_hmmManagerAlignment,
		m_hmmManagerAccumulation,bMultiplePronunciations,vLexUnitOptional);
	FBNodeHMM **nodes = hmmGraph->create(vLexUnitTranscription,&iNodes,&iEdges,&nodeInitial,&nodeFinal);
	delete hmmGraph;
	if (nodes == NULL) {
		*strReturnCode = FB_RETURN_CODE_UNABLE_TO_CREATE_HMM_GRAPH;
		return NULL;
	}

	assert(nodeInitial->iDistanceEnd % NUMBER_HMM_STATES == 0);
	
	//double dTimeMiddle1 = TimeUtils::getTimeMilliseconds();
	
	// there can't be fewer feature vectors than HMM-states in the composite
	if (iFeatures < nodeInitial->iDistanceEnd) {
		HMMGraph::destroy(nodes,iNodes);
		*strReturnCode = FB_RETURN_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS;
		return NULL;
	}	
	
	// reset the emission probability computation (to avoid using cached computations that are outdated)
	VHMMState vHMMState;
	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			assert(edge->hmmStateEstimation != NULL);
			vHMMState.push_back(edge->hmmStateEstimation);
		}
	}
	m_hmmManagerAlignment->resetHMMEmissionProbabilityComputation(vHMMState);
	
	// do the actual forward-backward
	FBTrellisNode *trellis = forwardBackward(iFeatures,fFeaturesAlignment,iNodes,nodes,iEdges,
		nodeInitial,nodeFinal,m_fForwardPruningBeam,m_fBackwardPruningBeam,strReturnCode);
	if (trellis == NULL) {
		HMMGraph::destroy(nodes,iNodes);
		return NULL;
	}
	
	// get the total forward score from the terminal edges
	double dForwardTotal = -DBL_MAX;
	for(FBEdgeHMM *edge = nodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		dForwardTotal = Numeric::logAddition(dForwardTotal,trellis[(iFeatures-1)*iEdges+edge->iEdge].dForward+
			trellis[(iFeatures-1)*iEdges+edge->iEdge].dBackward);
	}
	// get the total backward score from the initial edges
	double dBackwardTotal = -DBL_MAX;
	for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		dBackwardTotal = Numeric::logAddition(dBackwardTotal,trellis[edge->iEdge].dBackward+trellis[edge->iEdge].dForward);
	}
	// utterances that are too long can produce numerical inaccuracies during forward/backward
	//printf("%12.4f\n",fabs(dForwardTotal-dBackwardTotal));
	if (fabs(dForwardTotal-dBackwardTotal) > 0.1) {
		HMMGraph::destroy(nodes,iNodes);
		deleteTrellis(trellis);
		*strReturnCode = FB_RETURN_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES;
		return NULL;
	}
	
	// utterance likelihood
	assert(dForwardTotal != -DBL_MAX);
	double dLikelihoodUtterance = dForwardTotal;
	
	//double dTimeMiddle2 = TimeUtils::getTimeMilliseconds();
	
	//countUnusedPositions(trellis,iFeatures,iNodes);
	
	Alignment *alignment = new Alignment(ALIGNMENT_TYPE_FORWARD_BACKWARD);
	
	// single gaussian estimation (state-occupation, accumulate statistics in the logical HMM-accumulator)
	double dOccupationTotal = 0.0;
	if (bSingleGaussian) {
	
		FBEdgeHMM **edgeActiveCurrent = new FBEdgeHMM*[iEdges];
		FBEdgeHMM **edgeActiveNext = new FBEdgeHMM*[iEdges];
		int iActiveCurrent = 0;
		int iActiveNext = 0;
		
		// mark all the edges as inactive except the initial edge
		for(int i=0 ; i < iNodes ; ++i) {
			for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
				edge->iActive = -1;
			}
		}
		
		// activate the initial edges
		for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			edgeActiveCurrent[iActiveCurrent++] = edge;
			edge->iActive = 0;
		}	
		
		// traverse the trellis
		for(int t = 0 ; t < iFeatures ; ++t) {
			float *fFeatureVectorAccumulation = fFeaturesAccumulation+(t*m_iFeatureDimensionalityAccumulation);
			for(int i = 0 ; i < iActiveCurrent ; ++i) {
				int s = edgeActiveCurrent[i]->iEdge;
				FBTrellisNode *nodeTrellis = &trellis[t*iEdges+s];
				if ((nodeTrellis->dForward == -DBL_MAX) || (nodeTrellis->dBackward == -DBL_MAX)) {
					continue;
				}
				double dOccLikelihood = nodeTrellis->dForward+nodeTrellis->dBackward-dLikelihoodUtterance;
				double dOccProbability = exp(dOccLikelihood);
				assert(finite(dOccProbability));
				assert(!edgeActiveCurrent[i]->vAccumulator.empty());
				// accumulate statistics
				edgeActiveCurrent[i]->vAccumulator[0]->accumulateObservation(fFeatureVectorAccumulation,dOccProbability);
				dOccupationTotal += dOccProbability;
			}
			// determine active states for next time frame
			iActiveNext = 0;
			for(int i=0 ; i<iActiveCurrent ; ++i) {
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
		delete [] edgeActiveCurrent;
		delete [] edgeActiveNext;		
		
		//double dDifference = fabs(((double)iFeatures)-dOccupationTotal)/((double)iFeatures);
		//printf("%12d %12.2f -> %12.2f\n",iFeatures,dOccupationTotal,dDifference);
	}
	// multiple gaussian estimation (gaussian-occupation, accumulate statistics in the physical HMM-accumulator)
	else {
	
		FBEdgeHMM **edgeActiveCurrent = new FBEdgeHMM*[iEdges];
		FBEdgeHMM **edgeActiveNext = new FBEdgeHMM*[iEdges];
		int iActiveCurrent = 0;
		int iActiveNext = 0;
		
		// mark all the edges as inactive except the initial edge
		for(int i=0 ; i < iNodes ; ++i) {
			for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
				edge->iActive = -1;
			}
		}
		
		// activate the initial edges
		for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			edgeActiveCurrent[iActiveCurrent++] = edge;
			edge->iActive = 0;
		}	
		
		// traverse the trellis	
		double dOccupationTotal = 0.0;
		for(int t = 0 ; t < iFeatures ; ++t) {
			FrameAlignment *frameAlignment = new FrameAlignment;
			float *fFeatureVectorAlignment = fFeaturesAlignment+(t*m_iFeatureDimensionalityAlignment);
			float *fFeatureVectorAccumulation = fFeaturesAccumulation+(t*m_iFeatureDimensionalityAccumulation);
			for(int i = 0 ; i < iActiveCurrent ; ++i) {
				int s = edgeActiveCurrent[i]->iEdge;
				FBTrellisNode *nodeTrellis = &trellis[t*iEdges+s];
				if ((nodeTrellis->dForward == -DBL_MAX) || (nodeTrellis->dBackward == -DBL_MAX)) {
					continue;
				}
				double dOccupationLikelihood = 0.0;	
				// add forward score of predecessor edges and previous time frame (for t=0 the log(fwd score) is 0.0)
				if (t > 0) {
					dOccupationLikelihood = -DBL_MAX;
					// same edge
					if (trellis[(t-1)*iEdges+edgeActiveCurrent[i]->iEdge].dForward != -DBL_MAX) {
						dOccupationLikelihood = Numeric::logAddition(dOccupationLikelihood,trellis[(t-1)*iEdges+edgeActiveCurrent[i]->iEdge].dForward);
					}
					// predecessor edges
					for(FBEdgeHMM *edge = edgeActiveCurrent[i]->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
						if (trellis[(t-1)*iEdges+edge->iEdge].dForward != -DBL_MAX) {
							dOccupationLikelihood = Numeric::logAddition(dOccupationLikelihood,trellis[(t-1)*iEdges+edge->iEdge].dForward);	
						}
					}
					assert(dOccupationLikelihood != -DBL_MAX);
				}	
				dOccupationLikelihood += nodeTrellis->dBackward-dLikelihoodUtterance;
				
				double dOccupationProb = exp(nodeTrellis->dForward+nodeTrellis->dBackward-dLikelihoodUtterance);
				StateOcc *stateOcc = Alignment::newStateOcc(edgeActiveCurrent[i]->hmmStateEstimation->getId(),
					dOccupationProb);
				frameAlignment->push_back(stateOcc);	
				
				// compute the occupation for each gaussian in the mixture
				double dGaussianAux = 0.0;
				for(int g = 0 ; g < edgeActiveCurrent[i]->hmmStateEstimation->getMixture().getNumberComponents() ; ++g) {
					
					double dOccLikGaussian = edgeActiveCurrent[i]->hmmStateEstimation->computeEmissionProbabilityGaussian(g,fFeatureVectorAlignment,t) + dOccupationLikelihood + log(edgeActiveCurrent[i]->hmmStateEstimation->getMixture()(g)->weight());	
					
					double dOccProbability = exp(dOccLikGaussian);
					assert(finite(dOccProbability));
					dGaussianAux += dOccProbability;
					assert(!edgeActiveCurrent[i]->vAccumulator.empty());
					// global accumulators
					if (bAccumulatorsLogical) {
						edgeActiveCurrent[i]->vAccumulator[0]->accumulateObservation(fFeatureVectorAccumulation,dOccProbability);
					}
					// local accumulators
					else {
						edgeActiveCurrent[i]->vAccumulator[g]->accumulateObservation(fFeatureVectorAccumulation,dOccProbability);
					}
					dOccupationTotal += dOccProbability;
				}
				//printf("%8.4f\n",dGaussianAux);
			}
			// keep frame alignment
			alignment->addFrameAlignmentBack(frameAlignment);	
			// determine active states for next time frame
			iActiveNext = 0;
			for(int i=0 ; i<iActiveCurrent ; ++i) {
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
		delete [] edgeActiveCurrent;
		delete [] edgeActiveNext;
	}
		
	//double dTimeEnd = TimeUtils::getTimeMilliseconds();	
	
	//double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;	
	//double dTimeSeconds2 = (dTimeMiddle1-dTimeBegin)/1000.0;
	//double dTimeSeconds3 = (dTimeMiddle2-dTimeMiddle1)/1000.0;
	//double dTimeSeconds4 = (dTimeEnd-dTimeMiddle2)/1000.0;
	
	//printf("Xtterance probability: for= %12.4f back= %12.4f diff= %12.8f (%5d frames %5d states) time: %8.4fs %8.4fs\n",dForwardTotal,dBackwardTotal,fabs(dForwardTotal-dBackwardTotal),iFeatures,iEdges,dTimeSeconds,dTimeSeconds2);
	//printf("Xtterance probability: for= %12.4f back= %12.4f diff= %12.8f (%5d frames %5d states) time: %8.4fs %8.4fs %8.4fs %8.4fs\n",dForwardTotal,dBackwardTotal,fabs(dForwardTotal-dBackwardTotal),iFeatures,iEdges,dTimeSeconds,dTimeSeconds2,dTimeSeconds3,dTimeSeconds4);
	
	//print(trellis,iFeatures,iEdges);
	
	//exit(-1);
	
	// clean-up
	deleteTrellis(trellis);
	HMMGraph::destroy(nodes,iNodes);	
	
	*dUtteranceLikelihood = dLikelihoodUtterance;

	*strReturnCode = FB_RETURN_CODE_SUCCESS;
	
	return alignment;
}


// print the trellis
void ForwardBackwardX::print(FBTrellisNode *node, int iRows, int iColumns) {

	/*printf("\n");
	for(int i=0 ; i < iRows ; ++i) {
		for(int j=0 ; j < iColumns ; ++j) {
			printf("%12d  ",j);
		}
	}*/
	printf("\n");
	for(int i=0 ; i < iRows ; ++i) {
		for(int j=0 ; j < iColumns ; ++j) {
			float fScore = node[i*iColumns+j].fScore;
			if (fScore == -FLT_MAX) {
				fScore = -100000.0000;
			}
			double dForward = node[i*iColumns+j].dForward;
			if (dForward == -DBL_MAX) {
				dForward = -100000.0000;
			}
			double dBackward = node[i*iColumns+j].dBackward;
			if (dBackward == -DBL_MAX) {
				dBackward = -100000.0000;
			}	
			//printf("%12.4f  ",node[i*iColumns+j].fScore);
			//printf("%12.4f  ",node[i*iColumns+j].dForward);
			//printf("%12.4f  ",node[i*iColumns+j].dBackward);
			printf("%12.4f %12.4f %12.4f   ",fScore,dForward,dBackward);
			//printf("%12.4f %12.4f   ",fForward,fBackward);
		}
		printf("\n");
	}
}

// forward/backward computation on the trellis
FBTrellisNode *ForwardBackwardX::forwardBackward(int iFeatures, float *fFeatures, int iNodes, FBNodeHMM **nodes, int iEdges, FBNodeHMM *nodeInitial, FBNodeHMM *nodeFinal, float fBeamForward, float fBeamBackward, const char **strReturnCode) {

	//double dTimeBegin = TimeUtils::getTimeMilliseconds();

	// allocate memory for the trellis
	FBTrellisNode *trellis = newTrellis(iFeatures,iEdges,strReturnCode);
	if (trellis == NULL) {
		return NULL;
	}
	
	for(int i=0 ; i < iFeatures*iEdges ; ++i) {
		trellis[i].fScore = -FLT_MAX;
		trellis[i].dForward = -DBL_MAX;
		trellis[i].dBackward = -DBL_MAX;
	}

	// active states
	FBEdgeHMM **edgeActiveCurrent = new FBEdgeHMM*[iEdges];
	FBEdgeHMM **edgeActiveNext = new FBEdgeHMM*[iEdges];
	int iActiveCurrent = 0;
	int iActiveNext = 0;
	
	//double dTimeBeginBwd = TimeUtils::getTimeMilliseconds();	

	// (1) backward-pass
	
	// mark all the edges as inactive except the final edge
	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			edge->iActive = INT_MAX;
		}
	}	
	
	// activate predecessors of terminal edges
	int iTerminalEdges = 0;
	for(FBEdgeHMM *edge = nodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext, ++iTerminalEdges) {
		for(FBEdgeHMM *edge1 = edge->nodePrev->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
			edgeActiveCurrent[iActiveCurrent++] = edge1;
			edge1->iActive = iFeatures-2;
		}
		// self loop
		edgeActiveCurrent[iActiveCurrent++] = edge;
	}
	
	// initialize the backward score for the terminal edges	
	for(FBEdgeHMM *edge = nodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		trellis[((iFeatures-1)*iEdges)+edge->iEdge].dBackward = log(1.0/((float)iTerminalEdges));
	}	
	
	// fill the trellis	
	for(int t = iFeatures-2 ; t >= 0 ; --t) {
		// (a) compute backward scores for the time frame	
		float *fFeatureVector = fFeatures+((t+1)*m_iFeatureDimensionalityAlignment);
		for(int i = 0 ; i < iActiveCurrent ; ++i) {
			int s = edgeActiveCurrent[i]->iEdge;
			FBTrellisNode *nodeTrellis = &trellis[t*iEdges+s];
			// accumulate probability mass from next time-frame
			// (1) same node
			nodeTrellis->dBackward = -DBL_MAX;
			if (iFeatures-2-t >= edgeActiveCurrent[i]->nodeNext->iDistanceEnd) {
				FBTrellisNode *nodeSame = &trellis[(t+1)*iEdges+s];
				if (nodeSame->dBackward != -DBL_MAX) {
					nodeSame->fScore = edgeActiveCurrent[i]->hmmStateEstimation->computeEmissionProbability(fFeatureVector,t+1);
					nodeTrellis->dBackward = nodeSame->dBackward + nodeSame->fScore;
				}
			}
			// (2) next nodes
			for(FBEdgeHMM *edge = edgeActiveCurrent[i]->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
				if (iFeatures-2-t >= edge->nodeNext->iDistanceEnd) {
					FBTrellisNode *nodeNext = &trellis[(t+1)*iEdges+edge->iEdge];
					if (nodeNext->dBackward != -DBL_MAX) {
						nodeNext->fScore = edge->hmmStateEstimation->computeEmissionProbability(fFeatureVector,t+1);
						nodeTrellis->dBackward = Numeric::logAddition(nodeTrellis->dBackward,nodeNext->dBackward+nodeNext->fScore);
					}
				}
			}
		}
		// (b) beam pruning
		// - get the best scoring node
		double dBestScore = -DBL_MAX;
		for(int i = 0 ; i < iActiveCurrent ; ++i) {
			if (trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dBackward > dBestScore) {
				dBestScore = trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dBackward;
			}
		}	
		// - compute the threshold (minimum backward value for a node not to be pruned)
		double dThreshold = dBestScore - fBeamBackward;
		assert(dThreshold > -DBL_MAX);
		// (c) determine active states for previous time frame
		iActiveNext = 0;
		for(int i=0 ; i<iActiveCurrent ; ++i) {
			// actual pruning	
			if (trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dBackward < dThreshold) {
				trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dBackward = -DBL_MAX;
				continue;
			}
			// self-loop
			if (t > edgeActiveCurrent[i]->nodePrev->iDistanceStart) {
				if (edgeActiveCurrent[i]->iActive > t-1) {
					edgeActiveCurrent[i]->iActive = t-1;
					edgeActiveNext[iActiveNext++] = edgeActiveCurrent[i];
				}
			}
			// successors
			for(FBEdgeHMM *edge = edgeActiveCurrent[i]->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
				if (t > edge->nodePrev->iDistanceStart) {
					if (edge->iActive > t-1) {
						edge->iActive = t-1;
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
	
	//double dTimeBeginFwd = TimeUtils::getTimeMilliseconds();
	
	
	//print(trellis,iFeatures,iEdges);
	//exit(-1);
	
	// (2) forward-pass
	
	// compute the only scores that are left to compute 
	// (the backward process does not imply the computation of these scores)
	for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		trellis[edge->iEdge].fScore = edge->hmmStateEstimation->computeEmissionProbability(fFeatures,0);
	}	
	
	// get the utterance likelihood from the initial edges
	// HACK: this is not an exact method to compute the utterance likelihood, although it should be good enough
	double dLikelihoodUtterance = -DBL_MAX;
	for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		if (trellis[edge->iEdge].dBackward != -DBL_MAX) {
			dLikelihoodUtterance = Numeric::logAddition(dLikelihoodUtterance,trellis[edge->iEdge].dBackward+trellis[edge->iEdge].fScore);
		}
	}
	assert(dLikelihoodUtterance != -DBL_MAX);
	double dForwardThreshold = dLikelihoodUtterance + fBeamForward;	
	
	// reset the emission probability computation (to avoid using cached computations that are outdated)
	VHMMState vHMMState;
	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			assert(edge->hmmStateEstimation != NULL);
			vHMMState.push_back(edge->hmmStateEstimation);
		}
	}
	m_hmmManagerAlignment->resetHMMEmissionProbabilityComputation(vHMMState);	
	
	iActiveCurrent = 0;
	iActiveNext = 0;		
	
	//printf("\n");
	
	// mark all the edges as inactive except the initial edge
	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			edge->iActive = -1;
			//printf("%2d (%2d) %2d [%6d]\n",edge->nodePrev->iNode,edge->iEdge,edge->nodeNext->iNode,edge->hmmState->getId());
		}
	}
	
	// compute forward score of initial edges and activate their successors
	int iInitialEdges = 0;
	for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev, ++iInitialEdges) {
		assert(trellis[edge->iEdge].fScore != -FLT_MAX);
		//trellis[edge->iEdge].fScore = edge->hmmState->computeEmissionProbability(fFeatures,0);
		trellis[edge->iEdge].dForward = trellis[edge->iEdge].fScore;
		for(FBEdgeHMM *edge1 = edge->nodeNext->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
			edgeActiveCurrent[iActiveCurrent++] = edge1;
			edge1->iActive = 1;
		}
		// self loop
		edgeActiveCurrent[iActiveCurrent++] = edge;
		edge->iActive = 1;
	}
	
	for(FBEdgeHMM *edge = nodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		trellis[edge->iEdge].dForward += log(1.0/iInitialEdges);
	}	
	
	// fill the trellis	
	for(int t = 1 ; t < iFeatures ; ++t) {
		for(int i = 0 ; i < iActiveCurrent ; ++i) {		
			int s = edgeActiveCurrent[i]->iEdge;
			FBTrellisNode *nodeTrellis = &trellis[t*iEdges+s];
			if (nodeTrellis->dBackward == -DBL_MAX) {
				continue;
			}	
			// accumulate probability mass from previous time-frame
			nodeTrellis->dForward = -DBL_MAX;
			// (1) same node	
			if (t > edgeActiveCurrent[i]->nodePrev->iDistanceStart) {
				if (trellis[(t-1)*iEdges+s].dForward != -DBL_MAX) {
					nodeTrellis->dForward = trellis[(t-1)*iEdges+s].dForward;	
				}
			}
			// (2) previous nodes
			for(FBEdgeHMM *edge = edgeActiveCurrent[i]->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
				if (t > edge->nodePrev->iDistanceStart) {
					if (trellis[(t-1)*iEdges+edge->iEdge].dForward != -DBL_MAX) {
						nodeTrellis->dForward = Numeric::logAddition(trellis[(t-1)*iEdges+edge->iEdge].dForward,nodeTrellis->dForward);
					}
				}
			}
			if (nodeTrellis->dForward != -DBL_MAX) {	
				assert(nodeTrellis->fScore != -FLT_MAX);
				nodeTrellis->dForward += nodeTrellis->fScore;
				if (nodeTrellis->dForward+nodeTrellis->dBackward < dForwardThreshold) {
					nodeTrellis->dForward = -DBL_MAX;
				}
			}
		}
		// determine active states for next time frame
		iActiveNext = 0;
		for(int i=0 ; i<iActiveCurrent ; ++i) {
			// actual pruning	
			if (trellis[t*iEdges+edgeActiveCurrent[i]->iEdge].dForward == -DBL_MAX) {
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
	
	//double dTimeBwd = (dTimeBeginFwd-dTimeBeginBwd)/1000.0;
	//double dTimeFwd = (dTimeEnd-dTimeBeginFwd)/1000.0;
	//double dTimeTotal = (dTimeEnd-dTimeBegin)/1000.0;
	//printf("%8.4f %8.4f total: %8.4f \n",dTimeBwd,dTimeFwd,dTimeTotal);
	
	//print(trellis,iFeatures,iEdges);
	
	// clean-up
	delete [] edgeActiveCurrent;
	delete [] edgeActiveNext;	

	return trellis;
}


// allocate memory for the trellis
FBTrellisNode *ForwardBackwardX::newTrellis(int iFeatures, int iEdges, const char **strReturnCode) {

	// trellis size: # elements
	int iTrellisSize = iFeatures*iEdges;
	
	// check if the memory to be allocated does not exceed the maximum allowed
	int iTrellisSizeBytes = iTrellisSize*sizeof(FBTrellisNode);
	if (iTrellisSizeBytes > m_iTrellisMaxSizeBytes) {	
		*strReturnCode = FB_RETURN_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED;
		return NULL;
	}	
	
	// allocate memory for the trellis (if enabled, try to reuse a cached trellis)
	FBTrellisNode *trellis = NULL;
	if ((m_bTrellisCache == false) || (m_trellisCache == NULL) || (m_iTrellisCacheSize < (iFeatures*iEdges))) {
		// try to allocate memory for the trellis
		try {
			trellis = new FBTrellisNode[iFeatures*iEdges];
		}
		catch (const std::bad_alloc&) {
			*strReturnCode = FB_RETURN_CODE_UTTERANCE_TOO_LONG_UNABLE_TO_CREATE_TRELLIS;
			return NULL;
		}
		// if caching is enabled and the new trellis it not too large, cache it 
		// note: not caching large trellis prevents a more efficient use of the main memory
		if ((m_bTrellisCache) && (iTrellisSizeBytes <= m_iTrellisCacheMaxSizeBytes)) {	
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
	
	*strReturnCode = FB_RETURN_CODE_SUCCESS;
	return trellis;
}

// delete a trellis
void ForwardBackwardX::deleteTrellis(FBTrellisNode *trellis) {

	// clean-up
	if (m_trellisCache != trellis) {
		delete [] trellis;
	}	
}

// compute the number of unused positions in the trellis (used for analysis)
void ForwardBackwardX::countUnusedPositions(FBTrellisNode *trellis, int iFeatures, int iNodes) {

	int iTotal = iFeatures*iNodes;
	int iUnused = 0;

	for(int i=0 ; i < iFeatures ; ++i) {
		for(int j=0 ; j < iNodes ; ++j) {
			if ((trellis[i*iNodes+j].dForward == -DBL_MAX) && (trellis[i*iNodes+j].dBackward == -DBL_MAX)) {
				++iUnused;
			}
		}
	}
	
	double dPercent = 100.0*(((double)iUnused)/((double)iTotal));
	
	printf("total: %10d unused: %10d (%.2f) %10d %10d\n",iTotal,iUnused,dPercent,iFeatures,iNodes);
}

};	// end-of-namespace

