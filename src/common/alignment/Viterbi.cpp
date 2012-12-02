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


#include "BestPath.h"
#include "HypothesisLattice.h"
#include "PhoneSet.h"
#include "Viterbi.h"

namespace Bavieca {

// constructor
Viterbi::Viterbi(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager, float fBeamWidth)
{
	m_phoneSet = phoneSet;
	m_hmmManager = hmmManager;
	m_lexiconManager = lexiconManager;
	m_fBeamWidth = fBeamWidth;
	m_nodeTrellisCache = NULL;
	m_iTrellisCacheSize = 0;
	m_fScoreCache = NULL;
	m_iFeatureDimensionality = m_hmmManager->getFeatureDimensionality();
	
	m_iEntriesCacheOld = 0;
	m_fScoreCacheOld = NULL;
}

// destructor
Viterbi::~Viterbi()
{
	if (m_fScoreCacheOld != NULL) {
		delete [] m_fScoreCacheOld;
	}	
	if (m_nodeTrellisCache != NULL) {
		delete [] m_nodeTrellisCache;
	}
}

// align the feature vectors against the text
VPhoneAlignment *Viterbi::align(VLexUnit &vLexUnit, float *fFeatures, int iFeatureVectors, float *fLikelihood) {

	// reset the score cache
	m_fScoreCache = NULL;	

	// build the HMM-state composite from the reference
	VHMMStateDecoding vHMMStateDecodingComposite;
	getHMMStateDecodingComposite(vLexUnit,vHMMStateDecodingComposite);
	unsigned char iErrorCode;
	useCache(false);
	VPhoneAlignment *vPhoneAlignment = alignHMMStates(fFeatures,iFeatureVectors,
		vHMMStateDecodingComposite,vLexUnit,fLikelihood,0,iErrorCode);	

	return vPhoneAlignment;
}



// print a sequence of lexical units
void Viterbi::printLexUnitSequence(VLexUnit &vLexUnit) {

	for(VLexUnit::iterator it = vLexUnit.begin() ; it != vLexUnit.end() ; ++it) {
		if ((*it)->iPronunciation == 0) {
			printf("%s ",m_lexiconManager->getStrLexUnit((*it)->iLexUnit));
		} else {
			printf("%s(%d) ",m_lexiconManager->getStrLexUnit((*it)->iLexUnit),(*it)->iPronunciation+1);
		}
	}
	printf("\n");
}

// return an HMM-composite from a sequence of lexical units
void Viterbi::getHMMStateDecodingComposite(VLexUnit &vLexUnitText, VHMMStateDecoding &vHMMStateDecodingComposite, LexUnit *lexUnitLeft, LexUnit *lexUnitRight) {

	// (1) extract the sequence of HMM-states from the transcription	
	int iPhoneSilence = m_phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE);
	assert(iPhoneSilence != UCHAR_MAX);	
	unsigned char iPhonePrev = iPhoneSilence;
	if (lexUnitLeft != NULL) {
		assert(lexUnitLeft->vPhones.empty() == false);
		iPhonePrev = lexUnitLeft->vPhones.back();
	}
	unsigned char iPhoneNext = UCHAR_MAX;
	unsigned char iPosition = UCHAR_MAX;
	int iPhone = 0;
	// for each lexical unit
	for(VLexUnit::iterator it = vLexUnitText.begin() ; it != vLexUnitText.end() ; ++it, ++iPhone) {
		// for each phone
		assert(iPhonePrev != UCHAR_MAX);
		// iterate through the phones
		for(vector<int>::iterator jt = (*it)->vPhones.begin() ; jt != (*it)->vPhones.end() ; ++jt) {
			vector<int>::iterator kt = jt;
			advance(kt,1);
			// next phone within the current lexical unit
			if (kt != (*it)->vPhones.end()) {
				iPhoneNext = *kt;
			} 
			// look in the next lexical unit
			else {
				VLexUnit::iterator lt = it;
				advance(lt,1);
				// first phone of the next lexical unit
				if (lt != vLexUnitText.end()) {
					iPhoneNext = (*lt)->vPhones.front();	
				} 
				// there are no more phones: assume silence (or first phone from next lexical unit)
				else {
					if (lexUnitRight != NULL) {
						assert(lexUnitRight->vPhones.empty() == false);
						iPhoneNext = lexUnitRight->vPhones.front();
					} else {
						iPhoneNext = iPhoneSilence;
					}
				}
			}
			// monophone position
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
			// create the states corresponding to the phone
			for(int k=0 ; k < NUMBER_HMM_STATES ; ++k) {
				// get the most specific HMM-state available (it can be either a monophone or a clustered triphone)
				HMMStateDecoding *hmmStateDecoding = m_hmmManager->getHMMStateDecoding(&iPhonePrev,*jt,&iPhoneNext,iPosition,k);
				assert(hmmStateDecoding != NULL);
				vHMMStateDecodingComposite.push_back(hmmStateDecoding);
			}
			iPhonePrev = *jt;	
		}
	}
}

// compute emission probability or retrieves it from the cache
float Viterbi::computeEmissionProbability(HMMStateDecoding *hmmStateDecoding, float *fFeatureVector, int iFeatureVector) {

	// cache: check if the score exists in the cache
	if (m_bUseCache == true) {
	
		int iIndex = iFeatureVector*m_iHMMStatesPhysical+hmmStateDecoding->getId();
		if (m_fScoreCache[iIndex] != FLT_MAX) {
			return m_fScoreCache[iIndex];
		}
	
		// compute the score and keep it in the cache
		//m_fScoreCache[iIndex] = hmmState->computeEmissionProbability(fFeatureVector,iFeatureVector);	
		m_fScoreCache[iIndex] = hmmStateDecoding->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,iFeatureVector);	
	
		return m_fScoreCache[iIndex];
	}
	// no cache: compute the score directly
	else {
		return hmmStateDecoding->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,iFeatureVector);;	
	}
}

// align a sequence of HMM-states to the audio
VPhoneAlignment *Viterbi::alignHMMStates(float *fFeatureVectors, int iFeatureVectors, VHMMStateDecoding &vHMMStateDecodingComposite, VLexUnit &vLexUnit, float *fLikelihood, int iOffset, unsigned char &iErrorCode) {

	// check whether there are enough feature frames to perform the alignment
	if ((int)vHMMStateDecodingComposite.size() > iFeatureVectors) {
		iErrorCode = ALIGNER_ERROR_CODE_INSUFFICIENT_NUMBER_FEATURE_VECTORS;
		BVC_WARNING << "not enough feature frames to perform the alignment against the sequence of HMM-states";
		return NULL;
	}

	// reset the emission probability computation (to avoid using cached computations that are outdated)	
	m_hmmManager->resetHMMEmissionProbabilityComputation(vHMMStateDecodingComposite);
	
	// create the cache if necessary
	if (m_bUseCache == true) {
		if (m_fScoreCache == NULL) {
			m_iHMMStatesPhysical = m_hmmManager->getNumberHMMStatesPhysical();
			int iEntriesNeeded = iFeatureVectors*m_iHMMStatesPhysical;
			// (1) try to reuse an old cache
			if (m_iEntriesCacheOld >= iEntriesNeeded) {
				assert(m_fScoreCacheOld != NULL);
				m_fScoreCache = m_fScoreCacheOld;
				for(int i=0 ; i<iEntriesNeeded ; ++i) {
					m_fScoreCache[i] = FLT_MAX;
				}
			}
			// (2) create a brand new cache
			else {
				m_fScoreCache = new float[iEntriesNeeded];
				for(int i=0 ; i<iEntriesNeeded ; ++i) {
					m_fScoreCache[i] = FLT_MAX;
				}
				if (m_fScoreCacheOld != NULL) {
					delete [] m_fScoreCacheOld;
				}	
				m_fScoreCacheOld = m_fScoreCache;
				m_iEntriesCacheOld = iEntriesNeeded;
			}
		} 
	}
	
	// (2) create the trellis
	
	// print the HMM-state composite	
	int iHMMStates = vHMMStateDecodingComposite.size();
	
	// allocate memory for the trellis (try to reuse a cached trellis first)
	VNode *nodeTrellis = NULL;
	if ((m_nodeTrellisCache == NULL) || (m_iTrellisCacheSize < (iHMMStates*iFeatureVectors))) {
	
		// try to allocate memory for the trellis
		try {
			nodeTrellis = new VNode[iHMMStates*iFeatureVectors];
		}
		catch (const std::bad_alloc&) {
			iErrorCode = ALIGNER_ERROR_CODE_UTTERANCE_TOO_LONG_INSUFFICIENT_MEMORY;
			BVC_WARNING << "unable to allocate memory for the trellis, utterance too long?";
			return NULL;
		}	
		delete [] m_nodeTrellisCache;
		m_nodeTrellisCache = nodeTrellis;
		m_iTrellisCacheSize = iHMMStates*iFeatureVectors;
	} else {
		nodeTrellis = m_nodeTrellisCache;
	} 
	// initialization
	for(int t = 0 ; t < iFeatureVectors ; ++t) {
		for(int i = 0 ; i < iHMMStates ; ++i) {
			nodeTrellis[t*iHMMStates+i].dScore = -DBL_MAX;
		}
	}
	
	float fBeam = m_fBeamWidth;
	
	// (3) fill the trellis
	
	// t = 0
	nodeTrellis[0].dScore = computeEmissionProbability(vHMMStateDecodingComposite[0],&(fFeatureVectors[0]),0+iOffset);
	nodeTrellis[0].iHMMStatePrev = -1;
	// t > 0
	for(int t = 1 ; t < iFeatureVectors ; ++t) {
		double dScoreFrameBest = -FLT_MAX;	// keep the best frame-level score (pruning)
		float *fFeatureVector = fFeatureVectors+(t*m_iFeatureDimensionality);
		// only one predecessor
		if (((iFeatureVectors-1)-t)+1 > (iHMMStates-1)) {
			if (nodeTrellis[((t-1)*iHMMStates)].dScore != -DBL_MAX) {
				nodeTrellis[t*iHMMStates].dScore = nodeTrellis[((t-1)*iHMMStates)].dScore + 
				computeEmissionProbability(vHMMStateDecodingComposite[0],fFeatureVector,t+iOffset);
				nodeTrellis[t*iHMMStates].iHMMStatePrev = 0;
				if (nodeTrellis[t*iHMMStates].dScore > dScoreFrameBest) {
					dScoreFrameBest = nodeTrellis[t*iHMMStates].dScore;
				}
			}
		}
		// two possible predecessors
		for(int i = 1 ; i < iHMMStates ; ++i) {
			// avoid unnecessary computations	
			if ((i > t) || (((iFeatureVectors-1)-t) < ((iHMMStates-1)-i))) {
				continue;
			}
			VNode &node = nodeTrellis[(t*iHMMStates)+i];
			VNode &nodePrev1 = nodeTrellis[((t-1)*iHMMStates)+(i-1)];
			VNode &nodePrev2 = nodeTrellis[((t-1)*iHMMStates)+i];
			// two possible predecessors
			if (nodePrev2.dScore == -DBL_MAX) {
				if (nodePrev1.dScore == -DBL_MAX) {
					node.dScore = -DBL_MAX;
				} else {
					node.dScore = nodePrev1.dScore + 
					computeEmissionProbability(vHMMStateDecodingComposite[i],fFeatureVector,t+iOffset);
					node.iHMMStatePrev = i-1;
					if (node.dScore > dScoreFrameBest) {	
						dScoreFrameBest = node.dScore;
					}
				}
			} else {
				if (nodePrev1.dScore == -DBL_MAX) {
					node.dScore = nodePrev2.dScore + 
					computeEmissionProbability(vHMMStateDecodingComposite[i],fFeatureVector,t+iOffset);
					node.iHMMStatePrev = i;
					if (node.dScore > dScoreFrameBest) {	
						dScoreFrameBest = node.dScore;
					}
				} else {
					double dPredecessor1 = nodePrev2.dScore;
					double dPredecessor2 = nodePrev1.dScore;
					if (dPredecessor1 > dPredecessor2) {
						node.dScore = dPredecessor1;
						node.iHMMStatePrev = i;
					} else {
						node.dScore = dPredecessor2;
						node.iHMMStatePrev = i-1;
					}	
					node.dScore += computeEmissionProbability(vHMMStateDecodingComposite[i],fFeatureVector,t+iOffset);
					if (node.dScore > dScoreFrameBest) {	
						dScoreFrameBest = node.dScore;
					}
				}		
			}
		}	
		// apply the pruning if it is enabled
		if (fBeam != -1) {
			// only one predecessor
			if (((iFeatureVectors-1)-t)+1 > (iHMMStates-1)) {
				if (nodeTrellis[t*iHMMStates].dScore < dScoreFrameBest - fBeam) {
					nodeTrellis[t*iHMMStates].dScore = -DBL_MAX;
				} 
			}
			// two possible predecessors
			for(int i = 1 ; i < iHMMStates ; ++i) {
				// avoid unnecessary computations	
				if ((i > t) || (((iFeatureVectors-1)-t) < ((iHMMStates-1)-i))) {
					continue;
				}	
				if (nodeTrellis[(t*iHMMStates)+i].dScore < dScoreFrameBest - fBeam) {
					nodeTrellis[(t*iHMMStates)+i].dScore = -DBL_MAX;
				} 
			}
		}
	}	
		
	// (4) recover the best sequence of states
	VNode *nodeAux = &(nodeTrellis[(iHMMStates*iFeatureVectors)-1]);
	*fLikelihood = (float)nodeAux->dScore;
	//printf("Alignment likelihood: %f\n",*fLikelihood);
	int t = iFeatureVectors-1;
	int iHMMState = iHMMStates-1;
	int iHMMStatePrev = iHMMStates-1;
	bool bNewPhone = false;
	float fLikelihoodPhone = *fLikelihood;
	int iLexUnitIndex = vLexUnit.size()-1; 
	int iPhonesLexUnit = vLexUnit[iLexUnitIndex]->vPhones.size();
	int iPhoneIndexLexUnit = iPhonesLexUnit-1;									// phone index within the current lexical unit
	// create the first alignment unit
	VPhoneAlignment *vPhoneAlignment = new VPhoneAlignment();
	PhoneAlignment *phoneAlignment = new PhoneAlignment;
	phoneAlignment->iPhone = vHMMStateDecodingComposite[iHMMState]->getPhone();
	if (vLexUnit[iLexUnitIndex]->vPhones.size() == 1) {
		phoneAlignment->iPosition = WITHIN_WORD_POSITION_MONOPHONE;
	} else {
		phoneAlignment->iPosition = WITHIN_WORD_POSITION_END;
	}
	phoneAlignment->lexUnit = vLexUnit[iLexUnitIndex];
	phoneAlignment->iStateEnd[vHMMStateDecodingComposite[iHMMState]->getState()] = t+iOffset;
	vPhoneAlignment->push_front(phoneAlignment);
	while(t >= 0) {
		// the phone changes
		if ((vHMMStateDecodingComposite[iHMMState]->getState() == NUMBER_HMM_STATES-1) && (bNewPhone == true)) {
			// increment the phone index within the current lexical unit (or restart it in case we moved to the next lexical unit)
			--iPhoneIndexLexUnit;
			if (iPhoneIndexLexUnit == -1) {
				iLexUnitIndex--;
				if (iLexUnitIndex == -1) {
					assert(0);
				}
				assert(iLexUnitIndex >= 0);
				iPhonesLexUnit = vLexUnit[iLexUnitIndex]->vPhones.size();
				iPhoneIndexLexUnit = iPhonesLexUnit-1;
			}
			// complete the last phone alignment data structure
			assert(vPhoneAlignment->empty() == false);
			vPhoneAlignment->front()->iStateBegin[0] = t+1+iOffset;
			vPhoneAlignment->front()->fLikelihood = fLikelihoodPhone - (float)nodeAux->dScore;
			fLikelihoodPhone = (float)nodeAux->dScore;	
			// create a new phone alignment
			phoneAlignment = new PhoneAlignment;
			phoneAlignment->iPhone = vHMMStateDecodingComposite[iHMMState]->getPhone();
			if (iPhonesLexUnit == 1) {
				phoneAlignment->iPosition = WITHIN_WORD_POSITION_MONOPHONE;
			} else if (iPhoneIndexLexUnit == 0) {
				phoneAlignment->iPosition = WITHIN_WORD_POSITION_START;
			} else if (iPhoneIndexLexUnit == iPhonesLexUnit-1) {
				phoneAlignment->iPosition = WITHIN_WORD_POSITION_END;
			} else {
				assert((iPhoneIndexLexUnit > 0) && (iPhoneIndexLexUnit < iPhonesLexUnit-1));
				phoneAlignment->iPosition = WITHIN_WORD_POSITION_INTERNAL;
			}
			phoneAlignment->lexUnit = vLexUnit[iLexUnitIndex];
			phoneAlignment->iStateEnd[vHMMStateDecodingComposite[iHMMState]->getState()] = t+iOffset;
			vPhoneAlignment->push_front(phoneAlignment);
			bNewPhone = false;	
		} else {
			// next time a final state is observed there will be a new phone
			if (vHMMStateDecodingComposite[iHMMState]->getState() == 0) {
				bNewPhone = true;
			}
			// the state changes
			if (iHMMStatePrev != iHMMState) {
				phoneAlignment->iStateBegin[vHMMStateDecodingComposite[iHMMStatePrev]->getState()] = t+1+iOffset;
				phoneAlignment->iStateEnd[vHMMStateDecodingComposite[iHMMState]->getState()] = t+iOffset;
			}	
		}
		t--;
		iHMMStatePrev = iHMMState;
		iHMMState = nodeAux->iHMMStatePrev;
		assert((iHMMState == iHMMStatePrev-1) || (iHMMState == iHMMStatePrev));
		assert(iHMMState < iHMMStates);
		nodeAux = &(nodeTrellis[(iHMMStates*t)+nodeAux->iHMMStatePrev]);
	}
	assert((vPhoneAlignment->front()->iPosition == WITHIN_WORD_POSITION_START) || (vPhoneAlignment->front()->iPosition == WITHIN_WORD_POSITION_MONOPHONE));	
	vPhoneAlignment->front()->iStateBegin[0] = 0+iOffset;
	vPhoneAlignment->front()->fLikelihood = fLikelihoodPhone;
	
	// sanity checks
	for(VPhoneAlignment::iterator it = vPhoneAlignment->begin() ; it != vPhoneAlignment->end() ; ++it) {
		assert(phoneAlignment->iStateEnd[0]-phoneAlignment->iStateBegin[0] >= 0);
		assert(phoneAlignment->iStateBegin[1]-phoneAlignment->iStateEnd[0] >= 1);
		assert(phoneAlignment->iStateEnd[1]-phoneAlignment->iStateBegin[1] >= 0);
		assert(phoneAlignment->iStateBegin[2]-phoneAlignment->iStateEnd[1] >= 1);
		assert(phoneAlignment->iStateEnd[2]-phoneAlignment->iStateBegin[2] >= 0);
	}	

	return vPhoneAlignment;
}

// return a state-level alignment given the BestPath
VPhoneAlignment *Viterbi::align(float *fFeatures, int iFeatures, BestPath *bestPath) {

	VPhoneAlignment *vPhoneAlignment = new VPhoneAlignment;
	
	// the cache is not neded
	useCache(false);
	
	float fLikelihood = 0.0;
	float fLikelihoodAux = 0.0;
	// get the lexical units in the best path
	LBestPathElement *lBestPathElements = bestPath->getBestPathElements();
	LexUnit *lexUnitPrev = m_lexiconManager->getLexUnitSilence();
	LexUnit *lexUnitNext = NULL;
	// for each lexical unit
	for(LBestPathElement::iterator it = lBestPathElements->begin() ; it != lBestPathElements->end() ; ++it) {
	
		// skip beginning and end of utterance
		if (m_lexiconManager->isSentenceDelimiter((*it)->lexUnit)) {
			continue;
		}
		
		VLexUnit vLexUnit;
		vLexUnit.push_back((*it)->lexUnit);
		// (1) get the state-level alignment of the lexical unit
		// get the context lexical units
		LBestPathElement::iterator jt = it;
		++jt;
		// check if the next element is the end of utterance
		assert(jt != lBestPathElements->end());
		if ((*jt)->lexUnit == m_lexiconManager->m_lexUnitEndSentence) {
			lexUnitNext = m_lexiconManager->getLexUnitSilence();
		} else {
			lexUnitNext = (*jt)->lexUnit;
		}
		// build the HMM-state composite from the reference
		VHMMStateDecoding vHMMStateDecodingComposite;
		getHMMStateDecodingComposite(vLexUnit,vHMMStateDecodingComposite,lexUnitPrev,lexUnitNext);
		// do the force alignment between the HMM-composite and the features
		unsigned char iErrorCode;
		VPhoneAlignment *vPhoneAlignmentLexUnit = alignHMMStates(fFeatures+((*it)->iFrameStart*m_iFeatureDimensionality),(*it)->iFrameEnd-(*it)->iFrameStart+1,vHMMStateDecodingComposite,vLexUnit,&fLikelihoodAux,0,iErrorCode);
		fLikelihood += fLikelihoodAux;
		//printf("likelihood: %10f %10f\n",fLikelihood,fLikelihoodAux);
		if (vPhoneAlignmentLexUnit == NULL) {
			return NULL;
		}
		// update context lexical unit
		lexUnitPrev = (*it)->lexUnit;	
		
		// (2) move the lexical unit alignment to the utterance alignment
		for(VPhoneAlignment::iterator jt = vPhoneAlignmentLexUnit->begin() ; jt != vPhoneAlignmentLexUnit->end() ; ++jt) {
			// fix offsets
			for(int i=0 ; i<NUMBER_HMM_STATES ; ++i) {
				(*jt)->iStateBegin[i] += (*it)->iFrameStart;
				(*jt)->iStateEnd[i] += (*it)->iFrameStart;
			}
			vPhoneAlignment->push_back(*jt);
		}
		delete vPhoneAlignmentLexUnit;
	}
	//AlignmentFile *aux = new AlignmentFile(m_phoneSet,m_lexiconManager);
	//aux->print(*vPhoneAlignment);	
	
	// sanity checks
	assert(vPhoneAlignment->front()->iStateBegin[0] == 0);
	assert(vPhoneAlignment->back()->iStateEnd[NUMBER_HMM_STATES-1] == iFeatures-1);

	return vPhoneAlignment;
}


// align each of the lexical units in the lattice to the given set of feature vectors and store the
// time alignment information into the edges 
bool Viterbi::align(float *fFeatures, int iFeatures, HypothesisLattice *hypothesisLattice) {

	// the lattice must be marked with HMMs
	if (hypothesisLattice->checkProperty(LATTICE_PROPERTY_HMMS,"yes") == false) {
		return false;
	}

	// the cache is needed since we align different lexical units against the same segment of features
	useCache(true);
	//useCache(false);
	allocateCache(iFeatures);

	// reset the score cache
	m_fScoreCache = NULL;	
	
	// get the edges
	int iEdges = -1;
	LEdge **edges = hypothesisLattice->getEdges(&iEdges);
	if (iEdges < 1) {
		return false;	
	}
	for(int i=0 ; i < iEdges ; ++i) {
		assert(edges[i]->iFrameEnd > edges[i]->iFrameStart);
	}
	
	// align each of the edges
	VLexUnit vLexUnit;
	for(int i=0 ; i < iEdges ; ++i) {
	
		LEdge *edge = edges[i];
		LexUnit *lexUnit = edge->lexUnit;
		vLexUnit.push_back(lexUnit);
		
		assert(edge->iPhones > 0);
		
		// build the HMM-state composite from the reference
		assert(edge->phoneAlignment != NULL);	
		VHMMStateDecoding vHMMStateDecodingComposite;
		for(int iPhone = 0 ; iPhone < edge->iPhones ; ++iPhone) {
			for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				HMMStateDecoding *hmmStateDecoding = m_hmmManager->getHMMStateDecoding(edge->phoneAlignment[iPhone].iHMMState[iState]);
				vHMMStateDecodingComposite.push_back(hmmStateDecoding);
			}
		}
		assert(vHMMStateDecodingComposite.empty() == false);
		
		// do the forced alignment between the HMM-composite and the features
		unsigned char iErrorCode;
		int iFeaturesSegment = edge->iFrameEnd-edge->iFrameStart+1;
		int iOffset = edge->iFrameStart;
		float fLikelihood = -FLT_MAX;
		VPhoneAlignment *vPhoneAlignment = alignHMMStates(fFeatures+(m_iFeatureDimensionality*iOffset),iFeaturesSegment,vHMMStateDecodingComposite,
			vLexUnit,&fLikelihood,iOffset,iErrorCode);
		assert(vPhoneAlignment != NULL);
		edge->fScoreAM = fLikelihood;	
		
		assert(edge->phoneAlignment != NULL);
		assert((int)vPhoneAlignment->size() == edge->iPhones);
		int iPhone = 0;
		for(VPhoneAlignment::iterator it = vPhoneAlignment->begin() ; it != vPhoneAlignment->end() ; ++it, ++iPhone) {
			for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				edge->phoneAlignment[iPhone].iStateBegin[iState] = (*it)->iStateBegin[iState];
				edge->phoneAlignment[iPhone].iStateEnd[iState] = (*it)->iStateEnd[iState];
			}
		}	
		
		AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
		vLexUnit.clear();
	}
	
	// update the lattice properties
	hypothesisLattice->setProperty(LATTICE_PROPERTY_AM_PROB,"yes");
	hypothesisLattice->setProperty(LATTICE_PROPERTY_PHONE_ALIGN,"yes");
	
	return true;	
}

// allocate memory for the cache
void Viterbi::allocateCache(int iFeatures) {

	m_iHMMStatesPhysical = m_hmmManager->getNumberHMMStatesPhysical();
	int iEntriesNeeded = iFeatures*m_iHMMStatesPhysical;

	// (1) try to reuse an old cache
	if (m_iEntriesCacheOld >= iEntriesNeeded) {
		assert(m_fScoreCacheOld != NULL);
		m_fScoreCache = m_fScoreCacheOld;
		for(int i=0 ; i<iEntriesNeeded ; ++i) {
			m_fScoreCache[i] = FLT_MAX;
		}
	}
	// (2) create a brand new cache
	else {
		m_fScoreCache = new float[iEntriesNeeded];
		for(int i=0 ; i<iEntriesNeeded ; ++i) {
			m_fScoreCache[i] = FLT_MAX;
		}
		if (m_fScoreCacheOld != NULL) {
			delete [] m_fScoreCacheOld;
		}	
		m_fScoreCacheOld = m_fScoreCache;
		m_iEntriesCacheOld = iEntriesNeeded;
	}
}

};	// end-of-namespace



