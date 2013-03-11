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
#include "ForwardBackward.h"
#include "Numeric.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
ForwardBackward::ForwardBackward(PhoneSet *phoneSet, HMMManager *hmmManagerEstimation, HMMManager *hmmManagerUpdate, float fForwardPruningBeam, float fBackwardPruningBeam, int iTrellisMaxSize, bool bTrellisCacheEnabled, int iTrellisCacheMaxSize)
{
	m_iFeatureDimensionality = hmmManagerEstimation->getFeatureDimensionality();
	m_iFeatureDimensionalityUpdate = hmmManagerUpdate->getFeatureDimensionality();
	m_phoneSet = phoneSet;
	m_hmmManagerEstimation = hmmManagerEstimation;
	m_hmmManagerUpdate = hmmManagerUpdate;
	m_fForwardPruningBeam = fForwardPruningBeam;
	m_fBackwardPruningBeam = fBackwardPruningBeam;
	m_iTrellisSizeMaxBytes = iTrellisMaxSize*1024*1024;
	m_bTrellisCacheEnabled = bTrellisCacheEnabled;
	// conver from MB to Bytes
	m_iTrellisCacheSizeMaxBytes = iTrellisCacheMaxSize*1024*1024;
	
	// trellis cache
	m_nodeTrellisCache = NULL;
	m_iTrellisCacheSize = 0;
	
	// likelihood cache
	m_mLikelihoodCache = NULL;
	
	// uterance processed counter
	m_iUtterancesProcessed = 0;
	m_iLatticesProcessed = 0;
}

// destructor
ForwardBackward::~ForwardBackward()
{
	if (m_nodeTrellisCache != NULL) {
		delete [] m_nodeTrellisCache;
	}	
	if (m_mLikelihoodCache != NULL) {
		delete [] m_mLikelihoodCache;
	}
}

// create the trellis for an utterance (this trellis will be used to perform the forward-backward algorithm)
Alignment *ForwardBackward::processUtterance(VLexUnit &vLexUnitTranscription, float *fFeaturesAlignment, 
	float *fFeaturesAccumulation, int iFeatures, double *dUtteranceLikelihood, const char **strReturnCode) {

	double dBegin = TimeUtils::getTimeMilliseconds();
	
	assert(fFeaturesAlignment == fFeaturesAccumulation);
	
	bool bSingleGaussian = m_hmmManagerEstimation->isSingleGaussian();
	bool bAccumulatorsLogical = m_hmmManagerEstimation->areAccumulatorsLogical();
	bool bTrellisNotCached = false;
	unsigned char iContextSizeAcc = m_hmmManagerEstimation->getContextSizeAccumulators();
	
	// increment the number of utterances processed
	++m_iUtterancesProcessed;
	
	// (0) make sure there is at least one lexical unit to align to
	if (vLexUnitTranscription.empty()) {
		*strReturnCode = FB_RETURN_CODE_EMPTY_TRANSCRIPTION;
		return NULL;	
	}

	// (2) extract the sequence of HMM-states from the transcription
	VHMMState vHMMStateCompositeEstimation;
	VHMMState vHMMStateCompositeUpdate;
	vector<VAccumulator> vAccumulatorComposite;
	unsigned char iPhoneSilence = m_phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE);
	assert(iPhoneSilence != UCHAR_MAX);
	unsigned char iPosition = UCHAR_MAX;
	unsigned char iPhonesPrev[HMM_CONTEXT_SIZE_MAX]; 
	unsigned char iPhonesNext[HMM_CONTEXT_SIZE_MAX]; 
	for(int i=0 ; i < HMM_CONTEXT_SIZE_MAX ; ++i) {
		iPhonesPrev[i] = iPhoneSilence;
	}
	
	// convert the sequence of lexical units to a sequence of phones
	
	vector<unsigned char> vPhone;				// phones
	vector<unsigned char> vPosition;			// within word position	
	
	// left context padding
	for(int i=0 ; i < iContextSizeAcc ; ++i) {
		vPhone.push_back(iPhoneSilence);
		vPosition.push_back(UCHAR_MAX);
	}	
	// for each lexical unit
	for(VLexUnit::iterator it = vLexUnitTranscription.begin() ; it != vLexUnitTranscription.end() ; ++it) {
		// for each phone	
		assert((*it)->vPhones.empty() == false);
		for(vector<int>::iterator jt = (*it)->vPhones.begin() ; jt != (*it)->vPhones.end() ; ++jt) {
			vPhone.push_back(*jt);
			// monophone
			if ((*it)->vPhones.size() == 1) {
				iPosition = WITHIN_WORD_POSITION_MONOPHONE;	
			}
			// starting position 
			else if (jt == (*it)->vPhones.begin()) {
				iPosition = WITHIN_WORD_POSITION_START;
			} 			
			else {
				vector<int>::iterator mt = jt;
				advance(mt,1);
				// final position
				if (mt == (*it)->vPhones.end()) {
					iPosition = WITHIN_WORD_POSITION_END;
				} 
				// internal position
				else {
					iPosition = WITHIN_WORD_POSITION_INTERNAL;
				}
			}
			vPosition.push_back(iPosition);	
		}	
	}
	// right context padding
	for(int i=0 ; i < iContextSizeAcc ; ++i) {
		vPhone.push_back(iPhoneSilence);
	}
	
	for(int i=iContextSizeAcc ; i < (int)vPhone.size()-iContextSizeAcc ; ++i) {	
		for(int j=1 ; j <= iContextSizeAcc ; ++j) {	
			iPhonesPrev[iContextSizeAcc-j] = vPhone[i-j];
			iPhonesNext[j-1] = vPhone[i+j];
		}
		// create the states corresponding to the phone
		for(int k=0 ; k < NUMBER_HMM_STATES ; ++k) {
			// get the most specific HMM-state available (it can be either a monophone or a clustered triphone)
			HMMState *hmmState = m_hmmManagerEstimation->getHMMState(iPhonesPrev,vPhone[i],iPhonesNext,vPosition[i],k);
			assert(hmmState != NULL);
			if (m_hmmManagerUpdate != m_hmmManagerEstimation) {
				HMMState *hmmState = m_hmmManagerUpdate->getHMMState(iPhonesPrev,vPhone[i],iPhonesNext,vPosition[i],k);
				vHMMStateCompositeUpdate.push_back(hmmState);
			} else {
				vHMMStateCompositeUpdate.push_back(hmmState);
			}	
			vHMMStateCompositeEstimation.push_back(hmmState);
			VAccumulator vAccumulator;
			if (bAccumulatorsLogical) {
				Accumulator *accumulator = m_hmmManagerUpdate->getAccumulator(iPhonesPrev,vPhone[i],iPhonesNext,vPosition[i],k); 
				assert(accumulator);
				vAccumulator.push_back(accumulator);
			} else {
				for(unsigned int g=0 ; g < hmmState->getMixture().getNumberComponents() ; ++g) {
					Accumulator *accumulator = m_hmmManagerUpdate->getAccumulator(hmmState->getId(),g); 
					assert(accumulator);
					vAccumulator.push_back(accumulator);
				}
			}
			vAccumulatorComposite.push_back(vAccumulator);
		}	
	}
	
	// reset the emission probability computation (to avoid using cached computations that are outdated)	
	m_hmmManagerEstimation->resetHMMEmissionProbabilityComputation(vHMMStateCompositeEstimation);
	
	int iHMMStates = (int)vHMMStateCompositeEstimation.size();
	
	// there can't be fewer feature vectors than HMM-states in the composite
	if (iFeatures < iHMMStates) {
		*strReturnCode = FB_RETURN_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS;
		return NULL;
	}

	// trellis size: # elements
	int iTrellisSize = iFeatures*iHMMStates;
	
	// check if the memory to be allocated does not exceed the maximum allowed
	int iTrellisSizeBytes = iHMMStates*iFeatures*sizeof(NodeTrellis);
	if (iTrellisSizeBytes > m_iTrellisSizeMaxBytes) {
		*strReturnCode = FB_RETURN_CODE_UTTERANCE_TOO_LONG_MAXIMUM_TRELLIS_SIZE_EXCEEDED;
		return NULL;
	}
	
	// allocate memory for the trellis (if enabled, try to reuse a cached trellis first)
	NodeTrellis *nodeTrellis = NULL;
	if ((m_bTrellisCacheEnabled == false) || (m_nodeTrellisCache == NULL) || (m_iTrellisCacheSize < (iHMMStates*iFeatures))) {
		// try to allocate memory for the trellis
		try {
			nodeTrellis = new NodeTrellis[iHMMStates*iFeatures];
		}
		catch (const std::bad_alloc&) {
			*strReturnCode = FB_RETURN_CODE_UTTERANCE_TOO_LONG_UNABLE_TO_CREATE_TRELLIS;
			return NULL;
		}
		// if caching is enabled and the new trellis it not too large, cache it 
		// note: not caching large trellis prevents a more efficient use of the main memory
		if ((m_bTrellisCacheEnabled) && (iTrellisSizeBytes <= m_iTrellisCacheSizeMaxBytes)) {	
			// empty the cache if necessary
			if (m_nodeTrellisCache != NULL) {
				delete [] m_nodeTrellisCache;
			}
			// cache the new trellis
			m_nodeTrellisCache = nodeTrellis;
			m_iTrellisCacheSize = iHMMStates*iFeatures;
		} else {
			bTrellisNotCached = true;
		}
	} else {
		nodeTrellis = m_nodeTrellisCache;
	} 
	// initialization
	for(int i=0 ; i < iTrellisSize ; ++i) {
		nodeTrellis[i].dForward = -DBL_MAX;
		nodeTrellis[i].dBackward = -DBL_MAX;
		nodeTrellis[i].fScore = -FLT_MAX;
	}
	
	// (3) compute forward/backward probabilities	
	float fBackwardThreshold = m_fBackwardPruningBeam;
	float f1 = 0.0;
	float f2 = 0.0;
	double dTimeSeconds = 0.0;

	// compute forward-backward scores
	computeFBTrellis(nodeTrellis,vHMMStateCompositeEstimation,iHMMStates,fFeaturesAlignment,iFeatures,fBackwardThreshold);
	//printTrellis(nodeTrellis,iRows,iColumns);
		
	double dEnd = TimeUtils::getTimeMilliseconds();
	dTimeSeconds = (dEnd-dBegin)/1000.0;
	
	// (2.1) sanity checks
	f1 = (float)nodeTrellis[iHMMStates*iFeatures-1].dForward;
	f2 = (float)nodeTrellis[0].dBackward + nodeTrellis[0].fScore;
	if (fabs(f1-f2) > 0.1) {
		*strReturnCode = FB_RETURN_CODE_UTTERANCE_TOO_LONG_NUMERICAL_INACCURACIES;	
	}
	//printf("utterance probability: for= %12.4f back= %12.4f diff= %12.8f (%5d frames %5d states) time: %8.4fs\n",f1,f2,fabs(f1-f2),iFeatureVectors,iHMMStates,dTimeSeconds);	
	
	// (3) compute estimation statistics
	double dProbabilityUtterance = nodeTrellis[iHMMStates*iFeatures-1].dForward;
	
	// update the accumulators of each of the HMM-states in the utterance

	// single gaussian estimation (state-occupation, either physical or logical acumulator)
	double dOccupationTotal = 0.0;
	Accumulator *accumulator = NULL;
	if (bSingleGaussian) {	
		bool bStateData = false;
		for(int i = 0 ; i < iHMMStates ; ++i) {
			if (bAccumulatorsLogical) {
				accumulator = vAccumulatorComposite[i].front();	
			} else {
				accumulator = vAccumulatorComposite[i].front();
			}
			// for each time frame
			for(int t = 0 ; t < iFeatures ; ++t) {
				// only if there is a significant occupation probability
				if ((nodeTrellis[t*iHMMStates+i].dForward != -DBL_MAX) && (nodeTrellis[t*iHMMStates+i].dBackward != -DBL_MAX)) {
					double dProbabilityOccupation = nodeTrellis[t*iHMMStates+i].dForward + nodeTrellis[t*iHMMStates+i].dBackward - dProbabilityUtterance;
					bStateData = true;
					double dOccupation = exp(dProbabilityOccupation);
					accumulator->accumulateObservation(&fFeaturesAccumulation[t*m_iFeatureDimensionalityUpdate],dOccupation);
					dOccupationTotal += dOccupation;
				}
			}
		}
		assert(bStateData);
	} 
	// multiple gaussian estimation (gaussian-occupation, accumulate statistics in the physical HMM-accumulator)
	else {
	
		bool bStateData = false;
		// for each HMM-state
		for(int i = 0 ; i < iHMMStates ; ++i) {
			// for each gaussian in the mixture
			for(unsigned int g = 0 ; g < vHMMStateCompositeUpdate[i]->getMixture().getNumberComponents() ; ++g) {	
				// for each time frame
				for(int t = 0 ; t < iFeatures ; ++t) {
					// only if there is a significant occupation probability
					if ((nodeTrellis[t*iHMMStates+i].dForward != -DBL_MAX) && (nodeTrellis[t*iHMMStates+i].dBackward != -DBL_MAX)) {
						// compute the occupation probability
						double dProbabilityOccupation = -dProbabilityUtterance + nodeTrellis[t*iHMMStates+i].dBackward;
						double dU = 0.0; 
						// t == 0
						if ((t == 0) && (i == 0)) {
							dU = 0.0;
						}
						// t > 0
						else {
							if (nodeTrellis[(t-1)*iHMMStates+i].dForward != -DBL_MAX) {
								dU = nodeTrellis[((t-1)*iHMMStates)+i].dForward;
								if (i > 0) {
									if (nodeTrellis[((t-1)*iHMMStates)+(i-1)].dForward != -DBL_MAX) {	
										dU = Numeric::logAddition(dU,nodeTrellis[((t-1)*iHMMStates)+(i-1)].dForward);
									}
								}
							} else {		
								if (i > 0) {
									if (nodeTrellis[((t-1)*iHMMStates)+(i-1)].dForward != -DBL_MAX) {	
										dU = nodeTrellis[((t-1)*iHMMStates)+(i-1)].dForward;
									}
								}	
							}
						}
						dProbabilityOccupation += dU + log(vHMMStateCompositeEstimation[i]->getMixture()(g)->weight());
						// compute the gaussian score
						dProbabilityOccupation += vHMMStateCompositeEstimation[i]->computeEmissionProbabilityGaussian(g,&(fFeaturesAlignment[t*m_iFeatureDimensionality]),t);
						
						bStateData = true;
						double dOccupation = exp(dProbabilityOccupation);
						assert(finite(dOccupation));	
						vAccumulatorComposite[i][g]->accumulateObservation(&fFeaturesAccumulation[t*m_iFeatureDimensionalityUpdate],dOccupation);
						dOccupationTotal += dOccupation;
					}
				}
			}
		}
		assert(bStateData);	
	}
	
	//printf("frames: %d occupation: %.2f\n",iFeatureVectors,dAux);
		
	*dUtteranceLikelihood = dProbabilityUtterance;
	
	//printf("%d accumulators\n",m_iAccumulatorsAllocated);
	
	/*double dEnd = TimeUtils::getTimeMilliseconds();
	double dMillisecondsInterval1 = dAfterFB - dBegin;	
	double dMillisecondsInterval2 = dEnd - dAfterFB;
	double dTotal = dEnd - dBegin;
	double dRTF = dTotal/(iFeatureVectors*10.0);*/
	//printf("processing time: %.4f %.4f seconds (total= %.4f)(RTF= %.4f)\n",dMillisecondsInterval1/1000.0,dMillisecondsInterval2/1000.0,dTotal/1000.0,dRTF);
	
	//printf("utterance probability: for= %12.4f back= %12.4f diff= %12.8f (%5d frames %5d states) time: %8.4fs\n",f1,f2,fabsf(f1-f2),iFeatureVectors,iHMMStates,dTimeSeconds);	
	
	// clean-up
	if (bTrellisNotCached) {
		delete [] nodeTrellis;
	}	

	*strReturnCode = FB_RETURN_CODE_SUCCESS;

	return NULL;;
} 

// print the given trellis
void ForwardBackward::printTrellis(NodeTrellis *node, int iRows, int iColumns) {

	double dProbabilityUtterance = node[iRows*iColumns-1].dForward;
	
	cout << "utterance probability: " << dProbabilityUtterance << endl << endl;

	// print froward and backward scores
	for(int j=0 ; j < iRows; ++j) {
		//printf("state: %d ",j);
		for(int i=0 ; i < iColumns; ++i) {
			double dForward = node[(iRows*i)+j].dForward;
			double dBackward = node[(iRows*i)+j].dBackward;
			if (dForward == -DBL_MAX) {
				dForward = 1000.0;
			}
			if (dBackward == -DBL_MAX) {
				dBackward = 1000.0;
			}
			/*if (dForward == DBL_MAX) {
				dForward = 0.0;
				dBackward = 0.0;
			}*/
			/*float fProbabilityOccupation = node[i*iRows+j].dForward + node[i*iRows+j].dBackward - fProbabilityUtterance;
			float fOccupation = exp(fProbabilityOccupation);
			if (node[(iRows*i)+j].dForward == DBL_MAX) {
				fOccupation = -1.0;
			}
			//printf("[(t=%5d) (%5.2f) %10.2f,%10.2f (%10.2f)] ",i,fOccupation,dForward,dBackward,node[(iRows*i)+j].fScore);*/
			printf("[t=%5d %12.4f %12.4f] ",i,dForward,dBackward);
			//printf("[(t=%5d) %10.2f,%10.2f (%10.2f)] ",i,dForward,dBackward,node[(iRows*i)+j].fScore);
			continue;
			if (node[(iRows*i)+j].fScore != -FLT_MAX) {
				printf("[(t=%5d) %10.2f %10.2f (%10.2f)] ",i,dForward,dBackward,node[(iRows*i)+j].fScore);
			} else {
				printf("[(t=%5d) %10.2f %10.2f (%10.2f)] ",i,dForward,dBackward,1000.0);
			}
			//printf("[(t=%5d) %5.2f] ",i,fOccupation);
		}
		cout << endl;
	}
}

// aligns a lattice against the given feature vectors using Posterior Probabilities as the edge occupation
// probability (discriminative training)
MOccupation *ForwardBackward::processLattice(HypothesisLattice *lattice, float *fFeatures, int iFeatures, 
	float fScaleAM, float fScaleLM, double &dLikelihood, bool bMMI, float fBoostingFactor, const char **strReturnCode) {
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	bool bAccumulatorsLogical = m_hmmManagerEstimation->areAccumulatorsLogical();
	assert(bAccumulatorsLogical == false);

	int iEdges = -1;
	LEdge **edges = lattice->getEdges(&iEdges);

	MOccupation *mOccupation = new MOccupation();
	
	// allocate a structure to keep the occupation for phonemes within each edge
	double **dOccupationEdges = new double*[iEdges];
	double *dScoreAMEdges = new double[iEdges];
	for(int i=0 ; i < iEdges ; ++i) {
		dScoreAMEdges[i] = 0.0;
		int iFramesEdge = edges[i]->iFrameEnd-edges[i]->iFrameStart+1;
		dOccupationEdges[i] = new double[iFramesEdge*NUMBER_HMM_STATES];
		for(int j=0 ; j < iFramesEdge*NUMBER_HMM_STATES ; ++j) {
			dOccupationEdges[i][j] = 0.0;
		}
	}
	
	// boosted MMI -> compute phone-accuracy for each phone in the lattice
	if (bMMI) {
		// get the time-alignment of the reference
		VLPhoneAlignment *vLPhoneAlignment = lattice->getBestPathAlignment();
		if (vLPhoneAlignment == NULL) {
			*strReturnCode = FB_RETURN_CODE_UNABLE_TO_GET_BEST_PATH_FROM_LATTICE;
			return NULL;
		}
		// compute phone accuracy respect to the reference
		lattice->computePhoneAccuracy(*vLPhoneAlignment,true,false);
		delete vLPhoneAlignment;
	}
	
	// create the likelihood cache
	m_mLikelihoodCache = new MLikelihood;
	
	// backward threshold
	float fBackwardThreshold = m_fBackwardPruningBeam;	
	
	double dOccupationTotalAll = 0.0;
	
	// for each edge (word) in the lattice
	for(int i=0 ; i < iEdges ; ++i) {
	
		LEdge *edge = edges[i];
		
		// compute state occupation for each phoneme within the edge
		int iOffset = 0;
		for(int iPhone = 0 ; iPhone < edge->iPhones ; ++iPhone) {
		
			LPhoneAlignment *phoneAlignment = &edge->phoneAlignment[iPhone];
		
			int iFeaturesPhone = phoneAlignment->iStateEnd[NUMBER_HMM_STATES-1]-phoneAlignment->iStateBegin[0]+1;
			int iFeaturesOffset = phoneAlignment->iStateBegin[0];
			float *fFeaturesPhone = fFeatures+(iFeaturesOffset*m_iFeatureDimensionality);
			
			VHMMState vHMMStateComposite;
			int iHMMStates = NUMBER_HMM_STATES;
			for(int iState = 0; iState < iHMMStates ; ++iState) {
				HMMState *hmmState = m_hmmManagerEstimation->getHMMState(phoneAlignment->iHMMState[iState]);
				if (hmmState == NULL) {
					*strReturnCode = FB_RETURN_CODE_UTTERANCE_UNKNOWN_LATTICE_HMMSTATE;
					return NULL;
				}
				vHMMStateComposite.push_back(hmmState);
			}
		
			NodeTrellis *nodeTrellis = newTrellis(iHMMStates,iFeaturesPhone);
			
			computeFBTrellis(nodeTrellis,vHMMStateComposite,iHMMStates,
				fFeaturesPhone,iFeaturesPhone,fBackwardThreshold,iFeaturesOffset);
				
			// (2.1) sanity checks
			//double d1 = nodeTrellis[iHMMStates*iFeaturesPhone-1].dForward;
			//double d2 = nodeTrellis[0].dBackward + nodeTrellis[0].fScore;
	
			// (3) compute estimation statistics
			double dLikelihoodPhone = nodeTrellis[iHMMStates*iFeaturesPhone-1].dForward;
			dScoreAMEdges[i] += dLikelihoodPhone;
			if (bMMI) {
				dScoreAMEdges[i] -= fBoostingFactor*edge->fPhoneAccuracy[iPhone];
			}
			
			// hmm-state occupation		
			for(int iState = 0 ; iState < iHMMStates ; ++iState) {
				// for each time frame
				for(int t = 0 ; t < iFeaturesPhone ; ++t) {
					// only if there is a significant occupation probability
					NodeTrellis &node = nodeTrellis[t*iHMMStates+iState];
					if ((node.dForward != -DBL_MAX) && (node.dBackward != -DBL_MAX)) {
						double dOccupationLikelihood = node.dForward + node.dBackward - dLikelihoodPhone;
						double dOccupationProb = exp(dOccupationLikelihood);
						dOccupationTotalAll += dOccupationProb;
						dOccupationEdges[i][iOffset+t*NUMBER_HMM_STATES+iState] = dOccupationProb;
					}
				}
			}
			iOffset += iFeaturesPhone*NUMBER_HMM_STATES;
			
			deleteTrellis(nodeTrellis);
		}	
		
		//printf("%12.4f %12.4f\n",edge->fScoreAM,dProbabilityWord);
	}
	
	// update phone-level acoustic scores in lattice edges
	for(int i=0 ; i < iEdges ; ++i) {
		edges[i]->fScoreAM = (float)dScoreAMEdges[i];
	}
	
	// compute edge posterior probabilities (edge occupation probabilities)
	//lattice->storeTextFormat("./lattice.txt");
	lattice->computeForwardBackwardScores(fScaleAM,fScaleLM);
	lattice->computePosteriorProbabilities();
	
	// multiply occupations by the corresponding edge posterior probability and keep them
	for(int i=0 ; i < iEdges ; ++i) {
		int iOffset = 0;
		for(int iPhone=0 ; iPhone < edges[i]->iPhones ; ++iPhone) {
			double dPhoneOcc = 0.0;
			LPhoneAlignment *phoneAlignment = edges[i]->phoneAlignment+iPhone;
			int iFeaturesPhone = phoneAlignment->iStateEnd[NUMBER_HMM_STATES-1]-phoneAlignment->iStateBegin[0]+1;
			for(int iState=0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				for(int t=0 ; t < iFeaturesPhone ; ++t) {	
					pair<int,int> timeState(phoneAlignment->iStateBegin[0]+t,phoneAlignment->iHMMState[iState]);
					double dOccupationProb = dOccupationEdges[i][iOffset+t*NUMBER_HMM_STATES+iState];
					// multiply by the edge occupancy
					dOccupationProb *= edges[i]->fPP;
					MOccupation::iterator it = mOccupation->find(timeState);
					if (it == mOccupation->end()) {
						mOccupation->insert(MOccupation::value_type(timeState,dOccupationProb));
					} else {
						it->second += dOccupationProb;
					}
					dPhoneOcc += dOccupationProb;
				}
			}
			//printf("phone occ: %12.4f (%12d)\n",dPhoneOcc,iFeaturesPhone);
			iOffset += iFeaturesPhone*NUMBER_HMM_STATES;
		}
	}	
	
	// compute the lattice likelihood
	dLikelihood = lattice->getLikelihood();
	
	// clean-up
	for(int i=0 ; i < iEdges ; ++i) {
		delete [] dOccupationEdges[i];
	}	
	delete [] dOccupationEdges;
	delete [] dScoreAMEdges;
	delete m_mLikelihoodCache;
	m_mLikelihoodCache = NULL;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	
	*strReturnCode = FB_RETURN_CODE_SUCCESS;
	
	//printf("%12d %12.6f\n",iFeatures,dOccupationTotalAll);
	
	//printf("# features: %d\n",iFeatures);
	printf("processing time: %.2f seconds (RTF: %6.2f)\n",dTimeSeconds,dTimeSeconds/(((double)iFeatures)/100.0));

	return mOccupation;
}


// given a trellis estimates forward-backward scores
void ForwardBackward::computeFBTrellis(NodeTrellis *nodeTrellis, VHMMState &vHMMStateComposite, int iHMMStates, 
	float *fFeatureVectors, int iFeatureVectors, float fBackwardThreshold, int iOffset) {

	// (3.1) backward probabilities
	NodeTrellis *node;
	NodeTrellis *nodeSuccessor1;
	NodeTrellis *nodeSuccessor2;	
	float *fFeatureVector = NULL;
	double dScoreFrameBest = -DBL_MAX;
	// t = T
	nodeTrellis[(iFeatureVectors*iHMMStates)-1].dBackward = 0.0;
	nodeTrellis[(iFeatureVectors*iHMMStates)-1].fScore = computeLikelihood(vHMMStateComposite[iHMMStates-1],
		&(fFeatureVectors[(iFeatureVectors-1)*m_iFeatureDimensionality]),iFeatureVectors-1+iOffset);
	
	// t < T
	for(int t = iFeatureVectors-1 ; t > 0 ; --t) {
		dScoreFrameBest = -DBL_MAX;	// keep the best frame-level score, paths with a score outside the beam will be pruned-off
		fFeatureVector = fFeatureVectors+(t*m_iFeatureDimensionality);
		// only one successor
		node = &nodeTrellis[(t*iHMMStates)-1];
		nodeSuccessor1 = &nodeTrellis[((t+1)*iHMMStates)-1];
		if ((t > (iHMMStates-1)) && (nodeSuccessor1->dBackward != -DBL_MAX)) {	
			nodeSuccessor1->fScore = computeLikelihood(vHMMStateComposite[iHMMStates-1],fFeatureVector,t+iOffset);		
			node->dBackward = nodeSuccessor1->dBackward + nodeSuccessor1->fScore;
			// keep the best frame-level score
			if (node->dBackward > dScoreFrameBest) {
				dScoreFrameBest = node->dBackward;
			}
		}		
		int iRoof = min(iHMMStates-1,t);
		int iFloor = max(0,(iHMMStates-1)-(iFeatureVectors)+t);
		for(int i = iFloor ; i < iRoof ; ++i) {
			node = &nodeTrellis[((t-1)*iHMMStates)+i];
			nodeSuccessor1 = &nodeTrellis[(t*iHMMStates)+i];
			nodeSuccessor2 = &nodeTrellis[(t*iHMMStates)+(i+1)];
			// two possible successors
			if (nodeSuccessor1->dBackward == -DBL_MAX) {
				if (nodeSuccessor2->dBackward == -DBL_MAX) {
					node->dBackward = -DBL_MAX;
				} else {	
					nodeSuccessor2->fScore = computeLikelihood(vHMMStateComposite[i+1],fFeatureVector,t+iOffset);
					node->dBackward = nodeSuccessor2->dBackward + nodeSuccessor2->fScore;
					// keep the best frame-level score
					if (node->dBackward > dScoreFrameBest) {
						dScoreFrameBest = node->dBackward;
					} 
				}
			} else {
				if (nodeSuccessor2->dBackward == -DBL_MAX) {
					nodeSuccessor1->fScore = computeLikelihood(vHMMStateComposite[i],fFeatureVector,t+iOffset);
					node->dBackward = nodeTrellis[(t*iHMMStates)+i].dBackward + nodeSuccessor1->fScore;
					// keep the best frame-level score
					if (node->dBackward > dScoreFrameBest) {
						dScoreFrameBest = node->dBackward;
					} 
				} else {
					nodeSuccessor1->fScore = computeLikelihood(vHMMStateComposite[i],fFeatureVector,t+iOffset);
					nodeSuccessor2->fScore = computeLikelihood(vHMMStateComposite[i+1],fFeatureVector,t+iOffset);
					node->dBackward = Numeric::logAddition(nodeSuccessor2->dBackward+nodeSuccessor2->fScore,
						nodeSuccessor1->dBackward+nodeSuccessor1->fScore);
					// keep the best frame-level score
					if (node->dBackward > dScoreFrameBest) {
						dScoreFrameBest = node->dBackward;
					} 
				}	
			}
		}
		// apply the backward pruning if it is enabled
		if (fBackwardThreshold != -1) {
			if (t > (iHMMStates-1)) {
				if ((nodeTrellis[(t*iHMMStates)-1].dBackward + fBackwardThreshold) < dScoreFrameBest) {
					nodeTrellis[(t*iHMMStates)-1].dBackward = -DBL_MAX;
				}
			}
			for(int i = 0 ; i < iHMMStates-1 ; ++i) {
				// avoid unnecessary computations
				if ((i >= t) || (((iFeatureVectors-1)-t) < ((iHMMStates-1)-i)-1)) {
					continue;
				}
				if ((nodeTrellis[((t-1)*iHMMStates)+i].dBackward + fBackwardThreshold) < dScoreFrameBest) {
					nodeTrellis[((t-1)*iHMMStates)+i].dBackward = -DBL_MAX;
				}
			}
		}
	}
	
	// compute the one score that is left to compute (the backward process does not imply the computation of this score)
	nodeTrellis[0].fScore = computeLikelihood(vHMMStateComposite[0],fFeatureVectors,0+iOffset);	
	
	// forward threshold	
	double dProbabilityUtterance = nodeTrellis[0].dBackward + nodeTrellis[0].fScore;
	double dForwardThreshold = dProbabilityUtterance + m_fForwardPruningBeam;
	
	// (3.2) forward probabilities
	NodeTrellis *nodePredecessor1;
	NodeTrellis *nodePredecessor2;
	// t = 0
	assert(nodeTrellis[0].fScore != -FLT_MAX);
	nodeTrellis[0].dForward = nodeTrellis[0].fScore;	
	// t > 0
	for(int t = 1 ; t < iFeatureVectors ; ++t) {
		node = &nodeTrellis[t*iHMMStates];
		nodePredecessor2 = &nodeTrellis[((t-1)*iHMMStates)];	
		// only one predecessor
		if (((iFeatureVectors-1)-t)+1 > (iHMMStates-1)) {
			if (node->dBackward == -DBL_MAX) {
				node->dForward = -DBL_MAX;
			} else {
				assert(node->fScore != -FLT_MAX);
				node->dForward = nodePredecessor2->dForward + node->fScore;	
				if ((node->dForward + node->dBackward) < dForwardThreshold) {	
					node->dForward = -DBL_MAX;
				}
			}
		}
		// two possible predecessors
		for(int i = 1 ; i < iHMMStates ; ++i) {
			node = &nodeTrellis[(t*iHMMStates)+i];
			nodePredecessor1 = &nodeTrellis[((t-1)*iHMMStates)+(i-1)];
			nodePredecessor2 = &nodeTrellis[((t-1)*iHMMStates)+i];
			// avoid unnecessary computations	
			if ((i > t) || (((iFeatureVectors-1)-t) < ((iHMMStates-1)-i))) {
				continue;
			}	
			if (node->dBackward == -DBL_MAX) {
				node->dForward = -DBL_MAX;
				continue;
			}
			// two possible predecessors
			if (nodePredecessor2->dForward == -DBL_MAX) {
				if (nodePredecessor1->dForward == -DBL_MAX) {
					node->dForward = -DBL_MAX;
				} else {
					node->dForward = nodePredecessor1->dForward + node->fScore;
					assert(node->fScore != -FLT_MAX);
					if ((node->dForward + node->dBackward) < dForwardThreshold) {	
						node->dForward = -DBL_MAX;
					}
				}
			} else {
				if (nodePredecessor1->dForward == -DBL_MAX) {
					node->dForward = nodePredecessor2->dForward + node->fScore;
					assert(node->fScore != -FLT_MAX);
					if ((node->dForward + node->dBackward) < dForwardThreshold) {	
						node->dForward = -DBL_MAX;
					}	
				} else {
					node->dForward = Numeric::logAddition(nodePredecessor1->dForward,nodePredecessor2->dForward);
					node->dForward += node->fScore;
					assert(node->fScore != -FLT_MAX);
					if ((node->dForward + node->dBackward) < dForwardThreshold) {	
						node->dForward = -DBL_MAX;
					}
				}		
			}
		}	
	}
	
	//printTrellis(nodeTrellis,iHMMStates,iFeatureVectors);
}

// allocate memory for the trellis
NodeTrellis *ForwardBackward::newTrellis(int iHMMStates, int iFeatures) {

	int iTrellisSizeBytes = 0;
	// trellis size: # elements
	int iTrellisSize = iFeatures*iHMMStates;
	
	// allocate memory for the trellis (if enabled, try to reuse a cached trellis)
	NodeTrellis *nodeTrellis = NULL;
	if ((m_bTrellisCacheEnabled == false) || (m_nodeTrellisCache == NULL) || (m_iTrellisCacheSize < (iHMMStates*iFeatures))) {
		// try to allocate memory for the trellis
		try {
			nodeTrellis = new NodeTrellis[iHMMStates*iFeatures];
		}
		catch (const std::bad_alloc&) {
			return NULL;
		}
		// if caching is enabled and the new trellis it not too large, cache it 
		// note: not caching large trellis prevents a more efficient use of the main memory
		if ((m_bTrellisCacheEnabled) && (iTrellisSizeBytes <= m_iTrellisCacheSizeMaxBytes)) {	
			// empty the cache if necessary
			if (m_nodeTrellisCache != NULL) {
				delete [] m_nodeTrellisCache;
			}
			// cache the new trellis
			m_nodeTrellisCache = nodeTrellis;
			m_iTrellisCacheSize = iHMMStates*iFeatures;
		}
	} else {
		nodeTrellis = m_nodeTrellisCache;
	} 
	// initialization
	for(int i=0 ; i < iTrellisSize ; ++i) {
		nodeTrellis[i].dForward = -DBL_MAX;
		nodeTrellis[i].dBackward = -DBL_MAX;
		nodeTrellis[i].fScore = -FLT_MAX;
	}
	
	return nodeTrellis;
}

// delete a trellis
void ForwardBackward::deleteTrellis(NodeTrellis *nodeTrellis) {

	// clean-up
	if (m_nodeTrellisCache != nodeTrellis) {
		delete [] nodeTrellis;
	}	
}

// forward-backward alignment preserving phone-boundaries
Alignment *ForwardBackward::processPhoneAlignment(float *fFeatures, int iFeatures, 
	VLPhoneAlignment *vLPhoneAlignment, double &dLikelihood, const char **strReturnCode) {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	bool bAccumulatorsLogical = m_hmmManagerEstimation->areAccumulatorsLogical();
	assert(bAccumulatorsLogical == false);
	dLikelihood = 0.0;

	Alignment *alignment = new Alignment(ALIGNMENT_TYPE_FORWARD_BACKWARD);
	double dOccupationTotal = 0.0;
	
	// backward threshold
	float fBackwardThreshold = m_fBackwardPruningBeam;
	
	// for each phoneme
	for(VLPhoneAlignment::iterator it = vLPhoneAlignment->begin() ; it != vLPhoneAlignment->end() ; ++it) {
	
		LPhoneAlignment *phoneAlignment = *it;
		
		int iFeaturesPhone = phoneAlignment->iStateEnd[NUMBER_HMM_STATES-1]-phoneAlignment->iStateBegin[0]+1;
		int iFeaturesOffset = phoneAlignment->iStateBegin[0];
		float *fFeaturesPhone = fFeatures+(iFeaturesOffset*m_iFeatureDimensionality);	
		
		VHMMState vHMMStateComposite;
		int iHMMStates = NUMBER_HMM_STATES;
		for(int iState = 0; iState < iHMMStates ; ++iState) {
			HMMState *hmmState = m_hmmManagerEstimation->getHMMState(phoneAlignment->iHMMState[iState]);
			if (hmmState == NULL) {
				*strReturnCode = FB_RETURN_CODE_UTTERANCE_UNKNOWN_LATTICE_HMMSTATE;
				return NULL;
			}
			vHMMStateComposite.push_back(hmmState);
		}
	
		NodeTrellis *nodeTrellis = newTrellis(iHMMStates,iFeaturesPhone);
		
		computeFBTrellis(nodeTrellis,vHMMStateComposite,iHMMStates,
			fFeaturesPhone,iFeaturesPhone,fBackwardThreshold,iFeaturesOffset);
			
		// (2.1) sanity checks
		//double d1 = nodeTrellis[iHMMStates*iFeaturesPhone-1].dForward;
		//double d2 = nodeTrellis[0].dBackward + nodeTrellis[0].fScore;

		// (3) compute estimation statistics
		double dLikelihoodPhone = nodeTrellis[iHMMStates*iFeaturesPhone-1].dForward;
		dLikelihood += dLikelihoodPhone;
		
		// for each time frame
		for(int t = 0 ; t < iFeaturesPhone ; ++t) {
			// hmm-state occupation
			FrameAlignment *frameAlignment = new FrameAlignment;
			for(int iState = 0 ; iState < iHMMStates ; ++iState) {
				// only if there is a significant occupation probability
				NodeTrellis &node = nodeTrellis[t*iHMMStates+iState];
				if ((node.dForward != -DBL_MAX) && (node.dBackward != -DBL_MAX)) {
					double dOccupationLikelihood = node.dForward + node.dBackward - dLikelihoodPhone;
					double dOccupationProb = exp(dOccupationLikelihood);
					StateOcc *stateOcc = new StateOcc;
					stateOcc->iHMMState = phoneAlignment->iHMMState[iState];
					stateOcc->dOccupation = dOccupationProb;
					frameAlignment->push_back(stateOcc);
					dOccupationTotal += dOccupationProb;
				}
			}
			alignment->addFrameAlignmentBack(frameAlignment);
		}
		
		deleteTrellis(nodeTrellis);
	}
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	
	*strReturnCode = FB_RETURN_CODE_SUCCESS;
	
	printf("processing time: %.2f seconds (RTF: %6.2f)\n",dTimeSeconds,dTimeSeconds/(((double)iFeatures)/100.0));

	return alignment;
}

// data conversion
Alignment *ForwardBackward::getAlignment(MOccupation *mOccupation, int iFrames) {

	Alignment *alignment = new Alignment(ALIGNMENT_TYPE_FORWARD_BACKWARD);
	FrameAlignment **frameAlignment = new FrameAlignment*[iFrames];
	for(int i=0 ; i < iFrames ; ++i) {
		frameAlignment[i] = new FrameAlignment;
	}
	int iEvents = 0;
	double dOccupation = 0.0;
	for(MOccupation::iterator it = mOccupation->begin() ; it != mOccupation->end() ; ++it) {
		StateOcc *stateOcc = new StateOcc;
		stateOcc->iHMMState = it->first.second;
		stateOcc->dOccupation = it->second;
		frameAlignment[it->first.first]->push_back(stateOcc);
		dOccupation += it->second;
		++iEvents;
	}
	for(int i=0 ; i < iFrames ; ++i) {
		alignment->addFrameAlignmentBack(frameAlignment[i]);
		//printf("%4d: %d\n",i,frameAlignment[i]->size());
	}
	//exit(-1);
		
	delete [] frameAlignment;

	return alignment;
}

// compute observation prob
float ForwardBackward::computeLikelihood(HMMState *hmmState, float *fFeatures, int iT) {

	// use cache?
	if (m_mLikelihoodCache != NULL) {
		
		pair<int,int> timeState(iT,hmmState->getId());
		MLikelihood::iterator it = m_mLikelihoodCache->find(timeState);
		if (it == m_mLikelihoodCache->end()) {
			float fLikelihood = hmmState->computeEmissionProbability(fFeatures,-1);
			m_mLikelihoodCache->insert(MLikelihood::value_type(timeState,fLikelihood));
			return fLikelihood;
		} else {
			return it->second;
		}	
	}
	// cache disabled: 
	else {
		return hmmState->computeEmissionProbability(fFeatures,-1);
	}
}

};	// end-of-namespace
