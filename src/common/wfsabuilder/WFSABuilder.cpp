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


#include "HMMManager.h"
#include "LMManager.h"
#include "PhoneSet.h"
#include "TimeUtils.h"
#include "WFSABuilder.h"

namespace Bavieca {

// contructor
WFSABuilder::WFSABuilder(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager, LMManager *lmManager, unsigned char iNGram, float fLMScalingFactor) {

	m_phoneSet = phoneSet;
	m_hmmManager = hmmManager;
	m_lexiconManager = lexiconManager;
	m_lmManager = lmManager;
	m_iNGram = iNGram;
	m_fLMScalingFactor = fLMScalingFactor;
	
	m_iLexUnitsTotal = m_lexiconManager->getLexiconSize();
	m_iHMMStatesTotal = 0;
	m_hmmManager->getHMMStates(&m_iHMMStatesTotal);	
}

// destructor
WFSABuilder::~WFSABuilder() {


}

// build the decoding network as a WFSAcceptor
WFSAcceptor *WFSABuilder::build() {

	// build a within-word or cross-word context dependent decoding network depending on the acoustic model context
	// within-word (also for context-independent)
	if (m_hmmManager->getContextModelingOrderHMM() == HMM_CONTEXT_MODELING_MONOPHONES) {
		return buildWordInternal();
	}
	// cross-word 
	else {
		return buildCrossWord();
	}
}

// build a word-internal context dependent decoding network or a context-independent one
WFSAcceptor *WFSABuilder::buildWordInternal2() {

	unsigned char iContextSize = m_hmmManager->getContextSizeHMM();
	assert(iContextSize == 0);
	
	/*for(int i=0 ; i < 21000000 ; ++i) {
		State *stateAux = new State;
		stateAux->vTransition.push_back(new Transition);
	}*/
	/*for(int i=0 ; i < 21000000 ; ++i) {
		State2 *stateAux = new State2;
		stateAux->transitions = new Transition2[1];	
	}*/
	//exit(-1);
	// (1) create the language model graph
	unsigned int iStatesG = 0;
	unsigned int iTransitionsG = 0;
	//VState vStateFinalG;
	State *stateRootG = NULL;
	State *stateFinal = NULL;
	State *states = NULL;
	if (buildG(&states,&iStatesG,&iTransitionsG,&stateRootG,&stateFinal) == false) {
		return NULL;
	}
	// create a list with the states in G sorted by topological order with respect to epsilon transitions

	// (2) build the lexical tree (only one tree needs to be built)
	map<LexUnit*,State*> mLexUnitState;
	unsigned int iTreeDepth = 0;
	State *stateRoot = buildWordInternalLexiconTree(mLexUnitState,&iTreeDepth);
	if (stateRoot == NULL) {
		return NULL;
	}
	//printStateTransitions(stateRoot);
	
	// (3) do the actual incremental composition
	unsigned int iStateIdFinal = iStatesG;
	LLeafToProcess lLeafToProcess;													// leaf queue
	State **stateWaiting = new State*[iTreeDepth];								// waiting states
	VTransition *transitionInsertion = new VTransition[iTreeDepth];		// transitions ready to be inserted
	for(unsigned int i=0 ; i < iTreeDepth ; ++i) {
		stateWaiting[i] = NULL;
	}
	
	// create an array to mark the already processed states
	bool *bStateProcessed = new bool[iStatesG];
	for(unsigned int i=0 ; i < iStatesG ; ++i) {
		bStateProcessed[i] = false;
	}
	
	// traverse G and apply incremental composition to each state
	MStateState mStateState;
	unsigned int iTransitionsProcessed = 0;
	unsigned int iStatesProcessed = 0;
	unsigned int iTransitionsCreated = 0;
	unsigned int iStatesCreated = 0;
	unsigned int iTransitionsRemovedFromG = 0;
	unsigned int iStatesReused = 0;
	LState lState;
	lState.push_back(stateRootG);
	bStateProcessed[stateRootG->iId] = true;
	while(lState.empty() == false) {
		
		// get the next state to process
		State *stateGFrom = lState.front();
		lState.pop_front();
		
		//printf("processing G state:\n");
	
		// process the state	in G by processing all the transitions that come from it
		bool bLexUnitTransitions = false;		// whether actual lexical unit transitions (transitions that are actual words) are seen
		
		//printf("G state to be processed has: %d transitions\n",stateGFrom->vTransition.size());
		
		// (3.1) for each lexical-unit transition coming from the state, find the corresponding leaf in the lexical tree
		VTransition vTransitionSurvive;
		for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
		
			// insert G destination states into the queue (if not already processed)
			if (bStateProcessed[(*it)->state->iId] == false) {
				lState.push_back((*it)->state);
				bStateProcessed[(*it)->state->iId] = true;
			}
			
			// epsilon transition: keep it in the final tree 
			if ((*it)->iSymbol == EPSILON_TRANSITION) {
				vTransitionSurvive.push_back(*it);
				//printf("epsilon processed\n");
			} 
			// beginning/end of sentence: keep it in the final tree
			else if (((int)(*it)->iSymbol == (m_lexiconManager->m_lexUnitBegSentence->iLexUnit|LEX_UNIT_TRANSITION)) || 
				(((int)(*it)->iSymbol == (m_lexiconManager->m_lexUnitEndSentence->iLexUnit|LEX_UNIT_TRANSITION)))) {
				vTransitionSurvive.push_back(*it);
				//printf("%s processed\n",m_lexiconManager->getStrLexUnit((*it)->iSymbol));	
			}	
			// non-epsilon transition: put leaves from each pronuciation of the lexical unit in the queue
			else {
				bLexUnitTransitions = true;
				LexUnitX *lexUnitX = m_lexiconManager->getLexUnit((*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
				for(VLexUnit::iterator jt = lexUnitX->vLexUnitPronunciations.begin() ; jt != lexUnitX->vLexUnitPronunciations.end() ; ++jt) {
					LeafToProcess *leaf = new LeafToProcess;
					leaf->lexUnit = *jt;
					leaf->iHMMStates = (*jt)->vPhones.size()*NUMBER_HMM_STATES;
					leaf->stateLeaf = mLexUnitState[*jt]; 
					leaf->iRightContextGroup = UCHAR_MAX;
					leaf->stateGTo = (*it)->state;
					leaf->fWeight = (*it)->fWeight;
					lLeafToProcess.push_back(leaf);	
					//m_lexiconManager->printLexUnit(*jt);
				}
				delete *it;
			}	
		}
		iTransitionsRemovedFromG += stateGFrom->vTransition.size()-vTransitionSurvive.size();
		// remove old G transitions
		stateGFrom->vTransition.clear();
		// add the G transitions that survive
		for(VTransition::iterator it = vTransitionSurvive.begin() ; it != vTransitionSurvive.end() ; ++it) {
			stateGFrom->vTransition.push_back(*it);
		}
		
		// sort the queue by state number
		lLeafToProcess.sort(WFSABuilder::compareStateNumber);
		
		// process leaves in the queue
		while(lLeafToProcess.empty() == false) {
		
			LeafToProcess *leaf = lLeafToProcess.front();
			lLeafToProcess.pop_front();
			
			// create the lexical unit transition and the state that goes to it
			State *stateLexUnit = newState(iStateIdFinal++,NULL,NULL);
			stateLexUnit->vTransition.push_back(newTransition(((unsigned int)leaf->lexUnit->iLexUnitPron)|LEX_UNIT_TRANSITION,
				leaf->fWeight,leaf->stateGTo));
			checkTransitionSymbolCorrectness(((unsigned int)leaf->lexUnit->iLexUnitPron)|LEX_UNIT_TRANSITION);
			++iStatesCreated;
			++iTransitionsCreated;
			
			// try to reuse a state
			MStateState::iterator ft = mStateState.find(stateLexUnit);
			if (ft != mStateState.end()) {
				//assert(ft->second->vTransition.size() == 1);
				//fWeightRest -= ft->second->front()->fWeight;
				for(VTransition::iterator gt = stateLexUnit->vTransition.begin() ; gt != stateLexUnit->vTransition.end() ; ++gt) {
					delete *gt;
				}
				delete stateLexUnit;
				iStateIdFinal--;
				--iStatesCreated;
				--iTransitionsCreated;
				++iStatesReused;
				stateLexUnit = ft->second;
				if (iStatesReused % 100 == 0) {
					printf("(1) reused: %u total: %u\n",iStatesReused,iStatesCreated);
				}	
			} else {
				//mStateState.insert(MStateState::value_type(stateLexUnit,stateLexUnit));	
			}
			
			int iHMMStates = leaf->iHMMStates;
			State *state = leaf->stateLeaf->stateParent;		// this state goes to the last HMM-transition
			unsigned int iSymbol = leaf->stateLeaf->transitionBW->iSymbol;
			
			// process states from the leaf to the root
			for(int i=iHMMStates-1 ; i >= 0 ; --i) {
			
				//printStateTransitions(state);
			
				// empty position: insert the state into the buffer
				if (stateWaiting[i] == NULL) {
					stateWaiting[i] = state;
					transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
					assert(i%NUMBER_HMM_STATES == (int)iSymbol%NUMBER_HMM_STATES);
					checkTransitionSymbolCorrectness(iSymbol);
					++iTransitionsCreated;
					// if the state goes to a G node we already know the destination state in the final graph
					if (state == leaf->stateLeaf->stateParent) {
						transitionInsertion[i].back()->state = stateLexUnit;
					}
				}
				// different state: move the original state along with its transitions to the final tree
				else if (stateWaiting[i] != state) {
					// (1) the state in the buffer becomes final
					State *stateFinal = newState(iStateIdFinal++,NULL,NULL);
					++iStatesCreated;
					assert(i > 0);		// the root state stays in the buffer until the end
					assert(transitionInsertion[i].empty() == false);
					for(VTransition::iterator jt = transitionInsertion[i].begin() ; jt != transitionInsertion[i].end() ; ++jt) {
						stateFinal->vTransition.push_back(*jt);
						// if no destination state: look for the final state that does have destination state
						if ((*jt)->state == NULL) {
							assert(i == iHMMStates-1);
							bool bFound = false;
							for(int j = iTreeDepth-2 ; j > i ; --j) {
								if (stateWaiting[j] != NULL) {
									while(i != j) {
										State *stateFinal2 = newState(iStateIdFinal++,NULL,NULL);
										++iStatesCreated;
										for(VTransition::iterator kt = transitionInsertion[j].begin() ; kt != transitionInsertion[j].end() ; ++kt) {
											stateFinal2->vTransition.push_back(*kt);
											assert((*kt)->state != NULL);
											// weight pushing
											//applyWeightPushing(*kt);
										}
										stateWaiting[j] = NULL;
										transitionInsertion[j].clear();
										// try to reuse the state
										MStateState::iterator ft = mStateState.find(stateFinal2);
										if (ft != mStateState.end()) {
											for(VTransition::iterator gt = stateFinal2->vTransition.begin() ; gt != stateFinal2->vTransition.end() ; ++gt) {
												delete *gt;
											}
											delete stateFinal2;
											iStateIdFinal--;
											--iStatesCreated;
											--iTransitionsCreated;
											++iStatesReused;
											stateFinal2 = ft->second;
											if (iStatesReused % 100 == 0) {
												printf("(2) reused: %u total: %u\n",iStatesReused,iStatesCreated);
											}
										} else {
											//mStateState.insert(MStateState::value_type(stateFinal2,stateFinal2));
										}
										assert(transitionInsertion[j-1].back()->state == NULL);
										transitionInsertion[j-1].back()->state = stateFinal2;
										--j;
									}
									bFound = true;
								}	
							}
							assert(bFound);
						}
						assert((*jt)->state != NULL);
						// weight pushing
						//applyWeightPushing(*jt);
					}
					transitionInsertion[i].clear();
					// (2) insert the state in the map for later reuse
					MStateState::iterator ft = mStateState.find(stateFinal);
					if (ft != mStateState.end()) {
						for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
							delete *gt;
						}
						delete stateFinal;
						iStateIdFinal--;
						--iStatesCreated;
						--iTransitionsCreated;
						++iStatesReused;
						stateFinal = ft->second;	
						if (iStatesReused % 100 == 0) {
							printf("(3) reused: %u total: %u\n",iStatesReused,iStatesCreated);
						}
					} else {
						//mStateState.insert(MStateState::value_type(stateFinal,stateFinal));	
					}
					assert(transitionInsertion[i-1].back()->state == NULL);
					transitionInsertion[i-1].back()->state = stateFinal;
					// (3) the new state is inserted into the buffer along with its transition
					stateWaiting[i] = state;
					transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
					assert(i%NUMBER_HMM_STATES == (int)iSymbol%NUMBER_HMM_STATES);
					checkTransitionSymbolCorrectness(iSymbol);
					++iTransitionsCreated;
					// if the state goes to a G node we already know the destination state in the final graph
					if (state == leaf->stateLeaf->stateParent) {
						transitionInsertion[i].back()->state = stateLexUnit;						
					} 
				} 
				// equal state: stop the right to left traversal
				else {
					// create a new transition
					transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
					assert(i%NUMBER_HMM_STATES == (int)iSymbol%NUMBER_HMM_STATES);
					checkTransitionSymbolCorrectness(iSymbol);
					++iTransitionsCreated;
					// if the state goes to a G node we already know the destination state in the final graph
					if (state == leaf->stateLeaf->stateParent) {
						transitionInsertion[i].back()->state = stateLexUnit;
					}	
					break;
				}
				
				if (i == 0) {
					assert(state->stateParent == NULL);
					break;	
				} 
				
				// move to the parent state
				iSymbol = state->transitionBW->iSymbol;
				state = state->stateParent;
			}
				
			delete leaf;
		}
		
		// insert all the waiting states in the final graph (if any)
		if (bLexUnitTransitions) {
			assert(stateWaiting[iTreeDepth-1] == NULL);
			assert(stateWaiting[0] != NULL);							// the root always has to be in the buffer	
			for(int i=iTreeDepth-2 ; i >= 0 ; --i) {
				if (stateWaiting[i] != NULL) {	
					while(i >= 0) {
						State *stateFinal = NULL;
						if (i > 0) {
							stateFinal = newState(iStateIdFinal++,NULL,NULL);
							++iStatesCreated;	
						} else {
							// (i == 0) connect the root to the final graph
							stateFinal = stateGFrom;	
						}
						for(VTransition::iterator kt = transitionInsertion[i].begin() ; kt != transitionInsertion[i].end() ; ++kt) {
							stateFinal->vTransition.push_back(*kt);
							assert((*kt)->state != NULL);
						}
						// insert the state in the map for later reuse
						MStateState::iterator ft = mStateState.find(stateFinal);
						if (ft != mStateState.end()) {
							for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
								delete *gt;
							}
							delete stateFinal;
							iStateIdFinal--;
							--iStatesCreated;
							--iTransitionsCreated;
							++iStatesReused;
							stateFinal = ft->second;	
							if (iStatesReused % 100 == 0) {
								printf("(4) reused: %u total: %u\n",iStatesReused,iStatesCreated);
							}
						} else {
							//mStateState.insert(MStateState::value_type(stateFinal,stateFinal));
						}
						if (i > 0) {
							assert(transitionInsertion[i-1].back()->state == NULL);
							transitionInsertion[i-1].back()->state = stateFinal;	
						}
						//printStateTransitions(stateFinal);
						stateWaiting[i] = NULL;
						transitionInsertion[i].clear();
						--i;
					}
					break;
				}
			}
			// sanity check
			for(int i=0 ; i < (int)iTreeDepth-1 ; ++i) {
				assert(stateWaiting[i] == NULL);
				assert(transitionInsertion[i].empty());
			}
		}	
		
		iTransitionsProcessed += stateGFrom->vTransition.size();
		iStatesProcessed++;
			
		//printStateTransitions(stateGFrom);
		//printf("G node processed! (%u states processed, %u transitions processed)\n",iStatesProcessed,iTransitionsProcessed);
		//printf("(%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
	}
	// sanity check
	for(unsigned int i=0 ; i < iStatesG ; ++i) {
		assert(bStateProcessed[i]);
	}
	delete [] bStateProcessed;
	delete [] stateWaiting;
	delete [] transitionInsertion;
	
	// destroy the lexical tree
	destroyTree(stateRoot);
	
	printf("# states reused: %u\n",iStatesReused);
	printf("ending!!! (%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
	//exit(-1);
	
	// sanity checks: 
	// (1) make sure there are not unconnected states
	// (2) make sure transition symbols are correct
	int iHMMStatesTotal = -1;
	m_hmmManager->getHMMStates(&iHMMStatesTotal);
	unsigned int iTransitionsSeen = 0;
	iStatesProcessed = 0;
	bStateProcessed = new bool[iStateIdFinal];
	for(unsigned int i=0 ; i < iStateIdFinal ; ++i) {
		bStateProcessed[i] = false;
	}
	assert(lState.empty());
	lState.push_back(stateRootG);
	bStateProcessed[stateRootG->iId] = true;
	++iStatesProcessed;
	while(lState.empty() == false) {
		
		// get the next state to process
		State *stateGFrom = lState.front();
		lState.pop_front();
		
		// (3.1) for each lexical-unit transition coming from the state, find the corresponding leaf in the lexical tree
		for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
		
			checkTransitionSymbolCorrectness((*it)->iSymbol);
		
			++iTransitionsSeen;
			if (bStateProcessed[(*it)->state->iId] == false) {
				lState.push_back((*it)->state);
				bStateProcessed[(*it)->state->iId] = true;
				++iStatesProcessed;
			}
		}	
	}
	for(unsigned int i=0 ; i < iStateIdFinal ; ++i) {
		assert(bStateProcessed[i]);
	}
	delete [] bStateProcessed;
	printf("transitions: %u seen, %u created, %u removed\n",iTransitionsSeen,iTransitionsCreated,iTransitionsRemovedFromG);
	assert(iTransitionsSeen == (iTransitionsCreated+iTransitionsG-iTransitionsRemovedFromG));
	
	unsigned int iStatesTotal = iStatesG+iStatesCreated;
	unsigned int iTransitionsTotal = iTransitionsG+iTransitionsCreated-iTransitionsRemovedFromG;
	float fNetworkSize = iStatesTotal*4+iTransitionsTotal*12;
	printf("total states: %u, total transitions: %u\n",iStatesTotal,iTransitionsTotal);
	printf("Expected network size: %f MB\n",fNetworkSize/(1024.0*1024.0));
	
	//exit(-1);
	
	// (4) equalize the acceptor input	
	//equalizeInput(stateRootG,iStatesTotal,iTransitionsTotal);
	
	// (5) transform the resulting graph to an acceptor that can be used for decoding
	WFSAcceptor *wfsAcceptor = optimize(stateRootG,iStatesTotal,iTransitionsTotal,stateFinal);

	return wfsAcceptor;
}


// build the grammar (from an n-gram) as an acceptor, it returns the initial state
bool WFSABuilder::buildG(State **statesG, unsigned int *iStatesG, unsigned int *iTransitionsG, State **stateInitialG, State **stateFinalG) {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	unsigned int iStatesCreated = 0;
	unsigned int iTransitionsCreated = 0;
	*stateFinalG = NULL;
	
	// get the list of fillers
	VLexUnit vLexUnitFiller;
	m_lexiconManager->getVLexUnitFiller(vLexUnitFiller);

	// uniform: a single state with one transition per word, same probability
	if (m_iNGram == LM_NGRAM_ZEROGRAM) {
	
		unsigned int iStateID = 0;	
		
		State *state = newState(iStateID++,NULL,NULL);
		++iStatesCreated;
		
		// get the unigrams
		int iUnigrams = -1;
		m_lmManager->getUnigrams(iUnigrams);
		
		// count the unigrams that are relevant in order to estimate the uniform probability
		unsigned int iUnigramsRelevant = 0;
		for(int i=0 ; i < iUnigrams ; ++i) {
		
			// skip unigrams for <UNK> and <s> lexical units
			if ((i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) ||
				(i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
				continue;
			}
			
			++iUnigramsRelevant;
		}	
		float fProbability = log10f(1.0/((float)iUnigramsRelevant));
		
		for(int i=0 ; i < iUnigrams ; ++i) {
		
			// note that all the pronunciation variants have the same insertion penalty associated
			LexUnit *lexUnit = m_lexiconManager->getLexUnit(i,0);
		
			// skip unigrams for <UNK> and <s> lexical units
			if ((i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) ||
				(i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
				continue;
			}
			Transition *transition = NULL;
			// </s> does not have insertion penalty
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				transition = newTransition(i|LEX_UNIT_TRANSITION,fProbability*m_fLMScalingFactor,state);
			} else {
				transition = newTransition(i|LEX_UNIT_TRANSITION,fProbability*m_fLMScalingFactor+lexUnit->fInsertionPenalty,state);
			}	
			state->vTransition.push_back(transition);
			++iTransitionsCreated;
		}
		
		// create a self-transition for each filler
		for(VLexUnit::iterator it = vLexUnitFiller.begin() ; it != vLexUnitFiller.end() ; ++it) {
			Transition *transitionFiller = newTransition((*it)->iLexUnit|LEX_UNIT_TRANSITION,(*it)->fInsertionPenalty,state);
			state->vTransition.push_back(transitionFiller);
			++iTransitionsCreated;
		}		
		
		*iStatesG = iStatesCreated;
		*iTransitionsG = iTransitionsCreated;
		*stateInitialG = state;
		*stateFinalG = state;
	} 
	// unigram: a single state and one transition per word (no backoffs)
	else if (m_iNGram == LM_NGRAM_UNIGRAM) {	
	
		unsigned int iStateID = 0;	
		
		State *state = newState(iStateID++,NULL,NULL);
		++iStatesCreated;
		
		// get the unigrams
		int iUnigrams = -1;
		Unigram *unigrams = m_lmManager->getUnigrams(iUnigrams);
		for(int i=0 ; i < iUnigrams ; ++i) {
		
			// skip unigrams for <UNK>, <s> and </s> lexical units
			if ((i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) ||
				(i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
				continue;
			}
			
			// note that all the pronunciation variants have the same insertion penalty associated
			LexUnit *lexUnit = m_lexiconManager->getLexUnit(i,0);
			
			Transition *transition = NULL;
			// </s> does not have insertion penalty
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				transition = newTransition(i|LEX_UNIT_TRANSITION,unigrams[i].fProbability*m_fLMScalingFactor,state);
			} else {
				transition = newTransition(i|LEX_UNIT_TRANSITION,unigrams[i].fProbability*m_fLMScalingFactor+lexUnit->fInsertionPenalty,state);
			}	
			state->vTransition.push_back(transition);
			++iTransitionsCreated;
		}
		
		// create a self-transition for each filler
		for(VLexUnit::iterator it = vLexUnitFiller.begin() ; it != vLexUnitFiller.end() ; ++it) {
			Transition *transitionFiller = newTransition((*it)->iLexUnit|LEX_UNIT_TRANSITION,(*it)->fInsertionPenalty,state);
			state->vTransition.push_back(transitionFiller);
			++iTransitionsCreated;
		}		
		
		*iStatesG = iStatesCreated;
		*iTransitionsG = iTransitionsCreated;
		*statesG = state;
		*stateInitialG = state;
		*stateFinalG = state;
	}
	// bigram: one state per unigram (+ unigram backoff state) and one transition per existing bigram + backoff transitions
	else if (m_iNGram == LM_NGRAM_BIGRAM) {
		
		unsigned int iStateID = 0;
	
		// get the unigrams
		int iUnigrams = -1;
		int iUnigramsDiscarded = 0;
		Unigram *unigrams = m_lmManager->getUnigrams(iUnigrams);
		
		// allocate memory to store the states (one state per unigram plus the bigram backoff)
		*statesG = new State[iUnigrams+1];	
		for(int i=0 ; i < iUnigrams+1 ; ++i) {
			(*statesG)[i].transitionBW = NULL;
			(*statesG)[i].stateParent = NULL;
		}
		
		// create a map from unigram index (lexical unit index) to state
		State **stateMap = new State*[iUnigrams];
		for(int i=0 ; i < iUnigrams ; ++i) {
			stateMap[i] = NULL;
		}		
		
		// create the back-off state (to handle all the bigram back-offs)
		State *stateBackoff = &(*statesG)[iStateID];
		stateBackoff->iId = iStateID++;
		++iStatesCreated;	
		
		// create one state per unigram
		for(int i=0 ; i < iUnigrams ; ++i) {
			// skip the unknown lexical unit
			if (i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
				++iUnigramsDiscarded;
				continue;
			}	
			
			// create the state
			State *state = &(*statesG)[iStateID];
			state->iId = iStateID++;
			++iStatesCreated;
			stateMap[i] = state;
		
			// keep the final state (end of sentence)
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				*stateFinalG = state;
			}
		}
		
		// create a transition from the back-off state to all the unigram states
		for(int i=0 ; i < iUnigrams ; ++i) {
			// skip lexical units for which no state was created
			if (stateMap[i] == NULL) {
				continue;
			}
			// skip the <s>
			if (i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) {
				continue;
			}
			
			// note that all the pronunciation variants have the same insertion penalty associated
			LexUnit *lexUnit = m_lexiconManager->getLexUnit(i,0);
						
			Transition *transition = NULL;
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				// no insertion penalty
				transition = newTransition(i|LEX_UNIT_TRANSITION,unigrams[i].fProbability*m_fLMScalingFactor,stateMap[i]);
			} else {
				// insertion penalty
				transition = newTransition(i|LEX_UNIT_TRANSITION,unigrams[i].fProbability*m_fLMScalingFactor+lexUnit->fInsertionPenalty,stateMap[i]);
			}	
			++iTransitionsCreated;
			stateBackoff->vTransition.push_back(transition);
			assert(stateBackoff != transition->state);
		}
				
		// create the transitions from each unigram state (unigram states are first in the vector)
		for(int i=0 ; i < iUnigrams ; ++i) {
			// skip lexical units for which no state was created
			if (stateMap[i] == NULL) {
				continue;
			}				
			// skip the </s> lexical unit (this is the final state)
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				continue;
			}
			
			// create a self-transition for each filler
			for(VLexUnit::iterator it = vLexUnitFiller.begin() ; it != vLexUnitFiller.end() ; ++it) {
				Transition *transitionFiller = newTransition((*it)->iLexUnit|LEX_UNIT_TRANSITION,(*it)->fInsertionPenalty,stateMap[i]);
				stateMap[i]->vTransition.push_back(transitionFiller);
				++iTransitionsCreated;
			}
			
			// get the existing bigrams for this unigram
			int iBigrams = unigrams[i].iBigrams;
			Bigram *bigrams = unigrams[i].bigrams;
			int iBigramsDiscarded = 0;
			for(int j=0; j<iBigrams ; ++j) {
				// skip lexical units for which no state was created
				if (stateMap[bigrams[j].iLexUnit] == NULL) {
					++iBigramsDiscarded;
					continue;
				}
				// skip useless transition <s> <s>
				if ((i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) && (i == bigrams[j].iLexUnit)) {
					++iBigramsDiscarded;	
					continue;
				}
				
				// note that all the pronunciation variants have the same insertion penalty associated
				LexUnit *lexUnit = m_lexiconManager->getLexUnit(bigrams[j].iLexUnit,0);
							
				Transition *transition = NULL;
				if ((bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) || 	
					(bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) {	
					// no insertion penalty
					transition = newTransition(bigrams[j].iLexUnit|LEX_UNIT_TRANSITION,
						bigrams[j].fProbability*m_fLMScalingFactor,stateMap[bigrams[j].iLexUnit]);
				} else {
					// insertion penalty
					transition = newTransition(bigrams[j].iLexUnit|LEX_UNIT_TRANSITION,
						bigrams[j].fProbability*m_fLMScalingFactor+lexUnit->fInsertionPenalty,stateMap[bigrams[j].iLexUnit]);
				}
				++iTransitionsCreated;
				stateMap[i]->vTransition.push_back(transition);
			}	
			// create a transition to the bigram backoff state (if necessary)
			if (iBigrams-iBigramsDiscarded < iUnigrams-iUnigramsDiscarded) {  // only if there are back-offs, this is the usual case
				// create a transition to the bigram back-off state
				Transition *transition = newTransition(EPSILON_TRANSITION,unigrams[i].fProbabilityBackoff*m_fLMScalingFactor,stateBackoff);
				stateMap[i]->vTransition.push_back(transition);
				++iTransitionsCreated;
				assert(stateMap[i] != transition->state);
			}
		}
		
		// create the initial state and connect it to the state (<s>)
		State *stateInitial = stateMap[m_lexiconManager->m_lexUnitBegSentence->iLexUnit];
		assert(stateInitial != NULL);
		
		*iStatesG = iStatesCreated;
		*iTransitionsG = iTransitionsCreated;
		*stateInitialG = stateInitial;
		
		delete [] stateMap;	
	}	
	// trigram: one state per existing bigram, one state per trigram backoff, one state for all the bigram backoffs
	else if (m_iNGram == LM_NGRAM_TRIGRAM) {
	
		// make sure there is one and only one initial state: either (</s>,<s>) or (<s>,<s>)
		Bigram *bigramAux1 = m_lmManager->getBigram(m_lexiconManager->getLexUnitId(LEX_UNIT_END_SENTENCE),
			m_lexiconManager->getLexUnitId(LEX_UNIT_BEGINNING_SENTENCE));
		Bigram *bigramAux2 = m_lmManager->getBigram(m_lexiconManager->getLexUnitId(LEX_UNIT_BEGINNING_SENTENCE),
			m_lexiconManager->getLexUnitId(LEX_UNIT_BEGINNING_SENTENCE));
		// report error
		if (((bigramAux1 == NULL) && (bigramAux2 == NULL)) || ((bigramAux1 != NULL) && (bigramAux2 != NULL))) {
			printf("unable to find the initial state of the WFST while building G, the initial state must be either (</s>,<s>) or (</s>,<s>)\n");	
			return false;
		}
	
		// get the unigrams (important: the array of unigrams has iVocabularySize index 
		// (so lexUnit Ids can be used to access it))
		int iUnigrams = -1;
		Unigram *unigrams = m_lmManager->getUnigrams(iUnigrams);
		// create a map from unigram index (lexical unit index) to state
		State **stateMapUnigram = new State*[iUnigrams];
		for(int i=0 ; i < iUnigrams ; ++i) {
			stateMapUnigram[i] = NULL;
		}
		
		// get the bigrams
		int iBigrams = -1;
		Bigram *bigrams = m_lmManager->getBigrams(iBigrams);
		// create a map from bigram index (lexical unit index) to bigram state
		State **stateMapBigram = new State*[iBigrams];
		for(int i=0 ; i < iBigrams ; ++i) {
			stateMapBigram[i] = NULL;
		}		
	
		// allocate memory to store the states (one state per unigram and bigram plus the bigram backoff)
		*statesG = new State[iUnigrams+iBigrams+1];	
		for(int i=0 ; i < iUnigrams+iBigrams+1 ; ++i) {
			(*statesG)[i].transitionBW = NULL;
			(*statesG)[i].stateParent = NULL;
		}	
		
		// map used to keep final states, no transitions will be created from these states
		//MStateBool mStateFinal;
	
		// TODO: there cannot be transitions from sentence marker states, for example </s> <s> cannot go to <s> <s> but end, the decoder is supposed to
		// be utterance based so the end of sentence will be the end of utterance, there cannot be end of sentence markers within the utterance
	
		unsigned int iStateID = 0;
		
		// create the bigram back-off state (to handle all the bigram back-offs)
		State *stateBigramBackoff = &(*statesG)[iStateID];
		stateBigramBackoff->iId = iStateID++;
		++iStatesCreated;
		
		// create the unigram states (one per unigram)
		for(int i=0 ; i < iUnigrams ; ++i) {	
		
			// -> skip the unknown lexical unit
			// -> the beginning of sentence bigram (</s>,<s>) or (<s>,<s>) is the initial state 
			// case 1) (</s>,<s>) -> (<s>) -> (w1) when the bigram (<s>,w1) does not exist
			// case 2) (</s>,<s>) -> (<s>,w1) when the bigram (<s>,w1) exists
			// -> the end of sentence unigram (</s>) is the final state
			// case 1) (w1,w2) -> (</s>) when (w2,</s>) does not exist
			// case 2) (w1,w2) -> (w2,</s>) -> (</s>) then (w2,</s>) exists
			if (i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
				continue;
			}
		
			State *state = &(*statesG)[iStateID];
			state->iId = iStateID++;
			++iStatesCreated;
			stateMapUnigram[i] = state;
			
			// </s> is the final state 
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				*stateFinalG = state;
				// no transitions from the final state
				continue;
			}			
			
			// create a transition to the bigram back-off (if necessary)
			if (unigrams[i].iBigrams < iUnigrams) {

				Transition *transition = newTransition(EPSILON_TRANSITION,unigrams[i].fProbabilityBackoff*m_fLMScalingFactor,stateBigramBackoff);
				++iTransitionsCreated;
				state->vTransition.push_back(transition);
			}			
			
			// create a self-transition for each filler
			for(VLexUnit::iterator it = vLexUnitFiller.begin() ; it != vLexUnitFiller.end() ; ++it) {
				Transition *transitionFiller = newTransition((*it)->iLexUnit|LEX_UNIT_TRANSITION,(*it)->fInsertionPenalty,state);
				state->vTransition.push_back(transitionFiller);
				++iTransitionsCreated;
			}
		}
			
		// create the bigram states (one per bigram)
		*stateInitialG = NULL;
		for(int i=0 ; i < iUnigrams ;  ++i) {
		
			// skip bigrams that come from the unknown lexical unit
			if (i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
				continue;
			}
		
			for(int j=0 ; j < unigrams[i].iBigrams ; ++j) {
		
				// skip bigrams that go to the unknown lexical unit 
				if (unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
					continue;
				}
				
				// skip bigrams that go to the start of sentence except the initial state (</s>,<s>) or (<s>,<s>)
				if ((unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) && 
					(i != m_lexiconManager->m_lexUnitEndSentence->iLexUnit) && 
					(i != m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					continue;
				}
				// skip bigrams that go to the end of sentence (there is only one final state </s>)
				if (unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
					continue;
				}
				// skip bigrams that come from the end of sentence except the initial state: (</s> <s>)
				if ((i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) && 
					(unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					continue;
				}
				// allow bigrams that come from the start of sentence except: (<s>,</s>) and (<s>,<unk>)
				if ((i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) && 
					((unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) ||
					(unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit))) {
					continue;
				}
			
				// get the position of the bigram in the global array of bigrams
				int iIndex = &(unigrams[i].bigrams[j])-bigrams;
				assert((iIndex >= 0) && (iIndex < iBigrams));	
				
				State *state = &(*statesG)[iStateID];
				state->iId = iStateID++;
				++iStatesCreated;
				stateMapBigram[iIndex] = state;
				
				// keep the initial state
				if (((i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) || 
					(i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) && 
					(unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					if (*stateInitialG != NULL) {
						printf("multiple initial state!!\n");
						return false;
					}
					*stateInitialG = state; 
				}
				
				// create a self-transition for each filler
				// (only if the bigram state is reached outputting an actual lexical unit)
				// IMP: do it for the first state too!!
				if ((state == *stateInitialG) || ((unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitBegSentence->iLexUnit) &&
					(unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitEndSentence->iLexUnit) &&
					(unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitUnknown->iLexUnit))) {	
					for(VLexUnit::iterator it = vLexUnitFiller.begin() ; it != vLexUnitFiller.end() ; ++it) {
						Transition *transitionFiller = newTransition((*it)->iLexUnit|LEX_UNIT_TRANSITION,(*it)->fInsertionPenalty,state);
						state->vTransition.push_back(transitionFiller);
						++iTransitionsCreated;
					}
				}
			}
		}
		if (*stateInitialG == NULL) {
			printf("Error: no initial state defined\n");
			return false;
		}
		
		// connect unigram states to brigram states
		for(int i=0 ; i < iUnigrams ; ++i) {	
		
			// skip non existent unigrams
			if (stateMapUnigram[i] == NULL) {
				continue;
			}
			
			// (</s>) is the final state
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				continue;
			}
			
			// create a transition to the bigram states that come from this unigram	
			for(int j=0 ; j < unigrams[i].iBigrams ; ++j) {
			
				// connect to the final state </s> ?
				if (unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
					
					// no insertion penalty
					float fProbability = unigrams[i].bigrams[j].fProbability*m_fLMScalingFactor;
					Transition *transition = newTransition(unigrams[i].bigrams[j].iLexUnit|LEX_UNIT_TRANSITION,fProbability,*stateFinalG);
					++iTransitionsCreated;
					stateMapUnigram[i]->vTransition.push_back(transition);
					continue;
				}	
			
				// get the position of the bigram in the global array of bigrams
				int iIndex = &(unigrams[i].bigrams[j])-bigrams;
				assert((iIndex >= 0) && (iIndex < iBigrams));	
				// skip non existent bigrams
				if (stateMapBigram[iIndex] == NULL) {
					continue;
				}	
				
				float fProbability = unigrams[i].bigrams[j].fProbability*m_fLMScalingFactor;
				// insertion penalty?
				if ((unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitBegSentence->iLexUnit) &&
					(unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) {	
					// note that all the pronunciation variants have the same insertion penalty associated
					LexUnit *lexUnit = m_lexiconManager->getLexUnit(unigrams[i].bigrams[j].iLexUnit,0);
					fProbability += lexUnit->fInsertionPenalty;
				} 
				
				Transition *transition = newTransition(unigrams[i].bigrams[j].iLexUnit|LEX_UNIT_TRANSITION,fProbability,stateMapBigram[iIndex]);
				++iTransitionsCreated;
				stateMapUnigram[i]->vTransition.push_back(transition);
			}
		}
		
		// create a transition from the bigram back-off state to all the unigram states
		for(int i=0 ; i < iUnigrams ; ++i) {
		
			// skip non existent unigrams
			if (stateMapUnigram[i] == NULL) {
				continue;
			}
			
			// do not create a back-off to the beginning of sentence
			if (i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) {
				continue;
			}
			
			float fProbability = unigrams[i].fProbability*m_fLMScalingFactor;
			// insertion penalty?
			if ((i != m_lexiconManager->m_lexUnitBegSentence->iLexUnit) &&
				(i != m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) {	
				// note that all the pronunciation variants have the same insertion penalty associated
				LexUnit *lexUnit = m_lexiconManager->getLexUnit(i,0);
				fProbability += lexUnit->fInsertionPenalty;
			} 
			
			Transition *transition = newTransition(i|LEX_UNIT_TRANSITION,fProbability,stateMapUnigram[i]);
			++iTransitionsCreated;
			stateBigramBackoff->vTransition.push_back(transition);
		}
		
		// create the transitions from each of the bigram states (bigram states are first in the vector)
		for(int i=0 ; i < iBigrams ; ++i) {
		
			// skip non-existent bigram states
			if (stateMapBigram[i] == NULL) {
				continue;
			}
			
			// skip bigram states that are final (word,</s>)
			//if (mStateFinal.find(stateMapBigram[i]) != mStateFinal.end()) {
			//	continue;
			//}
		
			// get the existing trigrams for this bigram
			int iTrigrams = bigrams[i].iTrigrams;
			Trigram *trigrams = bigrams[i].trigrams;
			for(int j=0; j<iTrigrams ; ++j) {
				// get the position of the bigram in the global array of bigrams
				Bigram *bigramAux = m_lmManager->getBigram(bigrams[i].iLexUnit,trigrams[j].iLexUnit);
				int iIndex = (bigramAux-bigrams);
				assert((iIndex >= 0) && (iIndex < iBigrams));
				
				// skip non existent bigrams states
				if (stateMapBigram[iIndex] == NULL) {
					continue;
				}
				
				float fProbability = trigrams[j].fProbability*m_fLMScalingFactor;
				// insertion penalty?
				if ((trigrams[j].iLexUnit != m_lexiconManager->m_lexUnitBegSentence->iLexUnit) &&
					(trigrams[j].iLexUnit != m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) {
					// note that all the pronunciation variants have the same insertion penalty associated
					LexUnit *lexUnit = m_lexiconManager->getLexUnit(trigrams[j].iLexUnit,0);		
					fProbability += lexUnit->fInsertionPenalty;
				} 				
				
				Transition *transition = newTransition(trigrams[j].iLexUnit|LEX_UNIT_TRANSITION,fProbability,stateMapBigram[iIndex]);
				++iTransitionsCreated;
				stateMapBigram[i]->vTransition.push_back(transition);
			}
			// create a transition to the right unigram state (with the backoff link)
			if (iTrigrams < iUnigrams) {
			
				if (stateMapUnigram[bigrams[i].iLexUnit] != NULL) {
					Transition *transition = newTransition(EPSILON_TRANSITION,bigrams[i].fProbabilityBackoff*m_fLMScalingFactor,
						stateMapUnigram[bigrams[i].iLexUnit]);
					stateMapBigram[i]->vTransition.push_back(transition);	
					++iTransitionsCreated;
				}	
			}
		}
		printf("states created: %d\n",iStatesCreated);
		for(unsigned int i=0 ; i < iStatesCreated ; ++i) {
			printStateTransitions(&((*statesG)[i]));
		}
		
		// (1) get the initial state, which is the bigram state corresponding to (</s> <s>) or (<s> <s>)
		// get the bigram (</s> <s>)
		/*Bigram *bigramAux = m_lmManager->getBigram(m_lexiconManager->getLexUnitId(LEX_UNIT_END_SENTENCE),m_lexiconManager->getLexUnitId(LEX_UNIT_BEGINNING_SENTENCE));
		// alternatively get the bigram (<s> <s>)
		if (bigramAux == NULL) {
			bigramAux = m_lmManager->getBigram(m_lexiconManager->getLexUnitId(LEX_UNIT_BEGINNING_SENTENCE),m_lexiconManager->getLexUnitId(LEX_UNIT_BEGINNING_SENTENCE));
		} 
		// report error
		if (bigramAux == NULL) {
			printf("unable to find the initial state of the WFST while building G\n");	
			return false;
		}
		// get the position of the bigram in the global array of bigrams
		int iIndex = bigramAux-bigrams;
		assert((iIndex >= 0) && (iIndex < iBigrams));
		// create the initial state
		State *stateInitial = stateMapBigram[iIndex];
		assert(stateInitial != NULL);*/
		
		// perform sanity checks to make sure that all the states created are connected (seen) and all the transitions too,
		bool *bStateProcessed = new bool[iStatesCreated];
		for(unsigned int i=0 ; i<iStatesCreated ; ++i) {
			bStateProcessed[i] = false;
		}
		unsigned int iStatesSeen = 0;
		unsigned int iTransitionSeen = 0;
		LState lState;
		lState.push_back(*stateInitialG);
		bStateProcessed[(*stateInitialG)->iId] = true;
		while(lState.empty() == false) {
			
			// get the next state to process
			State *stateGFrom = lState.front();
			lState.pop_front();
			
			for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
			
				// insert G destination states into the queue (if not already processed)
				if (bStateProcessed[(*it)->state->iId] == false) {
					lState.push_back((*it)->state);
					bStateProcessed[(*it)->state->iId] = true;
				}
				++iTransitionSeen;
			}
			++iStatesSeen;	
		}
		for(unsigned int i=0 ; i<iStatesCreated ; ++i) {
			if (bStateProcessed[i] == false) {
				printf("not seen: %u\n",i);
			}
			//assert(bStateProcessed[i]);
		}	
		assert(iStatesSeen == iStatesCreated);
		assert(iTransitionSeen == iTransitionsCreated);
		printf("building G: states created: %d states seen: %d\n",iStatesCreated,iStatesSeen);
		printf("building G: trans created: %d tran seen: %d\n",iTransitionsCreated,iTransitionSeen);
		delete [] bStateProcessed;
		
		delete [] stateMapUnigram;	
		delete [] stateMapBigram;	
		
		*iStatesG = iStatesCreated;
		*iTransitionsG = iTransitionsCreated;
		//*stateInitialG = stateInitial;
	}	
	// fourgram: one state per existing trigram + one state per existing bigram + one state per existing unigram
	else {
		assert(m_iNGram == LM_NGRAM_FOURGRAM);
		
		// not supported	
		return false;
	}
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	
	printf("-- G graph -----------------------------\n");
	printf(" n-gram: %s\n",m_lmManager->getStrNGram());
	printf(" # states:      %12u\n",iStatesCreated);
	printf(" # transitions: %12u\n",iTransitionsCreated);
	printf(" building time: %.2f seconds\n",(dTimeEnd-dTimeBegin)/1000.0);
	if (*stateFinalG == NULL) {
		printf(" warning: no final state was found!!\n");
	} else {
		printf(" final state found\n");
	}
	printf("----------------------------------------\n");	
	
	return true;
}

// return a list containing the left context of filler lexical units (including silence)
bool *WFSABuilder::getFillerFirstPhones(unsigned char *iElements) {

	// allocate and initialize memory for the right context
	bool *bPhone = new bool[m_phoneSet->size()];
	for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
		bPhone[i] = false;
	}
	
	VLexUnitX *lexiconX = m_lexiconManager->getLexiconXReference();
   for(VLexUnitX::iterator it = lexiconX->begin() ; it != lexiconX->end() ; ++it) {
   
		// these lexical units are not inserted in the lexicon WFST
		if (((*it)->iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit) || 
			((*it)->iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) ||
			((*it)->iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit))
		{
			continue;
		}
   	
   	// for each alternative pronunciation	
   	for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {	
			if (m_lexiconManager->isFiller(*jt) == false) {
				continue;
			}	
   		bPhone[(*jt)->vPhones.front()] = true;
			++(*iElements);
   	}
   } 	

	return bPhone;
}


// return a data structure containing the left context of G-nodes
// it represents the set of 
unsigned char **WFSABuilder::getLeftContextGNodes(State *states, unsigned int iStatesG, int iStateGFinal) {
	
	// determine left context for each state
	
	// (1) allocate memory to keep the left contexts, each node in G 
	bool **bLeftContext = new bool*[iStatesG];
	for(unsigned int i = 0 ; i < iStatesG ; ++i) {
		bLeftContext[i] = new bool[m_phoneSet->size()];
		for(unsigned int j = 0 ; j < m_phoneSet->size() ; ++j) {
			bLeftContext[i][j] = false;
		}
	}
	
	// keep track of the # of epsilon transitions that arrive at every state
	int *iEpsilonTransitions = new int[iStatesG];
	for(unsigned int i = 0 ; i < iStatesG ; ++i) {
		iEpsilonTransitions[i] = 0;
	}
	
	// traverse the states in G by topological order of epsilon transitions (states without an epsilon transition are processed first)
	// for each state in G find the lexical units that can come from it by looking at its transitions, each 
	// transition represents a lexical unit and all the pronunciations for the lexical unit need to be considered
	// for the left context of destination G-nodes
	list<pair<int,Transition*> > vTransitionEpsilon;
	for(int i=iStatesG-1 ; i >= 0 ; --i) {
	
		Transition *transitionEpsilon = NULL;
		for(VTransition::iterator it = states[i].vTransition.begin() ; it != states[i].vTransition.end() ; ++it) {
		
			// lex-unit transition (it can't be <s> or </s>)
			if (((*it)->iSymbol & LEX_UNIT_TRANSITION) && 
				((int)(*it)->iSymbol != (m_lexiconManager->m_lexUnitBegSentence->iLexUnit|LEX_UNIT_TRANSITION)) && 
				((int)(*it)->iSymbol != (m_lexiconManager->m_lexUnitEndSentence->iLexUnit|LEX_UNIT_TRANSITION))) {
				// get the lexical unit and for each pronunciation keep the last phone, which is a left context of
				// the destination G-node
				LexUnitX *lexUnitX = m_lexiconManager->getLexUnit((*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
				for(VLexUnit::iterator jt = lexUnitX->vLexUnitPronunciations.begin() ; jt != lexUnitX->vLexUnitPronunciations.end() ; ++jt) {
					bLeftContext[(*it)->state->iId][(*jt)->vPhones.back()] = true;
				}
			} 
			// epsilon transition (transfer all the contexts)
			else if ((*it)->iSymbol & EPSILON_TRANSITION) {
				iEpsilonTransitions[(*it)->state->iId]++;
				assert(transitionEpsilon == NULL);
				transitionEpsilon = *it;
				vTransitionEpsilon.push_back(pair<int,Transition*>(i,*it));
			}
		}
	}
	// TODO: check whether this is correct for n-grams higher than trigram, if two consecutive epsilon
	// transitions occur then this process is incorrect and a true topological order needs to be done
	// epsilon transitions need to be processed at the end, when every state has its final set of left contexts
	
	while(vTransitionEpsilon.empty() == false) {
		
		pair<int,Transition*> aux = vTransitionEpsilon.front();
		vTransitionEpsilon.pop_front();
		
		// only transfer the context when all the incoming epsilon transitions are processed
		if (iEpsilonTransitions[aux.first] == 0) {
			for(unsigned int j = 0 ; j < m_phoneSet->size() ; ++j) {
				if (bLeftContext[aux.first][j]) {
					bLeftContext[aux.second->state->iId][j] = true;
				}
			}
			--iEpsilonTransitions[aux.second->state->iId];
		} else {
			vTransitionEpsilon.push_back(aux);
		}		
	}	
	
	// sanity check, all the states except the final state must have at least one left-context 
	// (the first state will have silence+fillers)
	for(unsigned int i = 0 ; i < iStatesG ; ++i) {
		if (i == (unsigned int)iStateGFinal) {
			continue;
		}
		int iElements = 0;
		for(unsigned int j = 0 ; j < m_phoneSet->size() ; ++j) {
			if (bLeftContext[i][j]) {
				++iElements;
			}
		}
		assert(iElements > 0);
		if (iElements == 0) {
			printf("problem with state: %d\n",i);
		}
	}
	
	unsigned int iMemoryUsed1 = iStatesG*m_phoneSet->size()*sizeof(bool)+iStatesG*sizeof(bool*);
	unsigned int iMemoryUsed2 = iStatesG*sizeof(unsigned char*);
	
	// compact the data structure to save memory
	unsigned char **iLeftContext = new unsigned char*[iStatesG];
	for(unsigned int i = 0 ; i < iStatesG ; ++i) {
		// get the array size
		unsigned char iSize = 0;
		for(unsigned int j = 0 ; j < m_phoneSet->size() ; ++j) {
			if (bLeftContext[i][j]) {
				++iSize;
			}
		}
		// allocate memory for the array
		iLeftContext[i] = new unsigned char[iSize+1];
		iMemoryUsed2 += (iSize+1)*sizeof(unsigned char);
		unsigned char iOffset = 0;
		for(unsigned int j = 0 ; j < m_phoneSet->size() ; ++j) {
			if (bLeftContext[i][j]) {
				iLeftContext[i][iOffset++] = j;
			}
		}
		iLeftContext[i][iOffset] = UCHAR_MAX;	
		// release the old array
		delete [] bLeftContext[i];
	}	
	delete [] bLeftContext;
	
	printf("memory used: %u -> %u\n",iMemoryUsed1,iMemoryUsed2);
		
	return iLeftContext;
}




// create the word-internal/context-independent lexical tree
State *WFSABuilder::buildWordInternalLexiconTree(map<LexUnit*,State*> &mLexUnitLeafState, unsigned int *iTreeDepth) {

	unsigned int iStatesTotal = 0;
	unsigned int iTransitionsTotal = 0;
	unsigned char iPhoneLeftContext = m_phoneSet->getPhoneIndexSilence();
	State *stateRoot = newState(0,NULL,NULL);
	unsigned int iStateId = 1;	
	
	// for each lexical unit
	VLexUnitX *lexiconX = m_lexiconManager->getLexiconXReference();
	for(VLexUnitX::iterator it = lexiconX->begin() ; it != lexiconX->end() ; ++it) {
	
		// skip these lexical units
		if (((*it)->iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit) ||
			((*it)->iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) ||
			((*it)->iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit))
		{
			continue;
		}
		
		// for each alternative pronunciation	
		for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
		
			State *state = stateRoot;
			unsigned char iPhonePrev = iPhoneLeftContext;
		
			// for each phone
			unsigned char iPhoneNumber = 0;
			for(vector<int>::iterator kt = (*jt)->vPhones.begin() ; kt != (*jt)->vPhones.end() ; ++kt, ++iPhoneNumber) {	
			
				// get the next phone
				vector<int>::iterator lt = kt;
				advance(lt,1);	
				unsigned char iPhoneNext = UCHAR_MAX;
				
				// not the last phoneme within the lexical unit
				if (lt != (*jt)->vPhones.end()) {	
					iPhoneNext = *lt;
				} 
				// last phoneme within the lexical unit
				else {	
					iPhoneNext = m_phoneSet->getPhoneIndexSilence();
				}
				assert(iPhoneNext <= m_phoneSet->size());
				
				for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				
					// get the HMM-state
					unsigned char iPosition;
					if ((*jt)->vPhones.size() == 1) {
						iPosition = WITHIN_WORD_POSITION_MONOPHONE;
					} else if (iPhoneNumber == 0) {
						iPosition = WITHIN_WORD_POSITION_START;
					} else if (iPhoneNumber > 0) {
						iPosition = WITHIN_WORD_POSITION_INTERNAL;	
					} else if (iPhoneNumber == (*jt)->vPhones.size()-1) {
						iPosition = WITHIN_WORD_POSITION_END;	
					} else {
						assert(0);
					}
					
					unsigned int iHMMState = m_hmmManager->getHMMStateIndex(&iPhonePrev,*kt,&iPhoneNext,iPosition,iState);
					// try to reuse an already existing transition (same HMM-state)
					bool bFound = false;
					for(VTransition::iterator mt = state->vTransition.begin() ; mt != state->vTransition.end() ; ++mt) {
						if ((*mt)->iSymbol == iHMMState) {
							state = (*mt)->state;
							bFound = true;
							break;
						}
					}
					// cannot reuse, create a new transition and a new state
					if (bFound == false) {
						// create the next state
						Transition *transitionAux = newTransition(iHMMState,0.0,NULL);
						transitionAux->state = newState(iStateId++,transitionAux,state);
						++iStatesTotal;
						state->vTransition.push_back(transitionAux);
						++iTransitionsTotal;
						state = transitionAux->state;
					}
						
					// insert the leaf state in the map
					if ((iState == (NUMBER_HMM_STATES-1)) && (iPhoneNumber == ((*jt)->vPhones.size()-1))) {
						mLexUnitLeafState.insert(map<LexUnit*,State*>::value_type(*jt,state));
					}
				}
				iPhonePrev = *kt;	
			}
		}
	}
	
	// number the nodes of the tree in a post-order fashion
	iStateId = 0;
	preOrderNumbering(stateRoot,&iStateId);
	
	// get the maximum tree depth (all the trees have the same depth)
	*iTreeDepth = getTreeDepth(stateRoot);
	
	// count the number of leaves (#leaves <= #lexical units, this is because of homophony)
	LState lState;
	for(map<LexUnit*,State*>::iterator it = mLexUnitLeafState.begin() ; it != mLexUnitLeafState.end() ; ++it) {
		lState.push_back(it->second);
	}
	lState.sort();
	lState.unique();
	unsigned int iLeaves = lState.size();
	lState.clear();
	//unsigned int iLeaves2 = getLeaves(stateRoot);
	
	cout << "# states:        " << iStatesTotal << endl;
	cout << "# transitions:   " << iTransitionsTotal << endl;
	cout << "# lexical units: " << mLexUnitLeafState.size() << endl;
	cout << "# leaves:        " << iLeaves << endl;
	cout << "tree depth:      " << *iTreeDepth << endl;	

	return stateRoot;
}


// create the left-context-dependent map of trees
//LexiconTree **WFSABuilder::buildCrossWordLexiconTrees(unsigned int *iTreeDepth, unsigned char ***iRightContextGroups, unsigned int *iGroups) {
LexiconTree **WFSABuilder::buildCrossWordLexiconTrees(unsigned int *iTreeDepth, vector<unsigned char*> &vRightContextGroups) {

	unsigned int iTreesTotal = 0;
	unsigned int iStatesTotal = 0;
	unsigned int iTransitionsTotal = 0;
	unsigned int iPhonesLongerLexUnit = 0;
	unsigned char iRightContextSubset[m_phoneSet->size()+1];

	// (1) determine the set of left and right contexts for the given lexicon and context modeling order
	unsigned char iContextSize = m_hmmManager->getContextSizeHMM();
	assert(iContextSize == 1);
	
	// allocate and initialize memory for the right context
	bool *bPhoneLeft = new bool[m_phoneSet->size()];
	bool *bPhoneRight = new bool[m_phoneSet->size()];
	for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
		bPhoneLeft[i] = false;
		bPhoneRight[i] = false;
	}
	
	VLexUnitX *lexiconX = m_lexiconManager->getLexiconXReference();
   for(VLexUnitX::iterator it = lexiconX->begin() ; it != lexiconX->end() ; ++it) {
   
		// these lexical units are not inserted in the lexicon WFST
		if (((*it)->iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit) || 
			((*it)->iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) ||
			((*it)->iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit))
		{
			continue;
		}
   	
   	// for each alternative pronunciation	
   	for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {	
   		bPhoneLeft[(*jt)->vPhones.front()] = true;
   		bPhoneRight[(*jt)->vPhones.back()] = true;
   	}
   } 
 
   // count different left and right contexts
   unsigned char iLeftContexts = 0;
   unsigned char iRightContexts = 0;
   for(unsigned char i = 0 ; i< m_phoneSet->size() ; ++i) {
  		if (bPhoneRight[i]) {
  			iLeftContexts++;
  		}
  		if (bPhoneLeft[i]) {
  			iRightContexts++;
  		}
  	}
	
	// (2) create a lexicon prefix tree for each possible left context
	
	State **stateLeaf = new State*[m_phoneSet->size()];
	
	
	//*iRightContextGroups = new unsigned char*[MAX_RIGHT_CONTEXT_GROUPS];
	//for(unsigned int i=0 ; i < MAX_RIGHT_CONTEXT_GROUPS ; ++i) {
	//	(*iRightContextGroups)[i] = NULL;
	//}
	//*iGroups = 0;
	
   // allocate memory for the array of trees
	LexiconTree **lexiconTree = new LexiconTree*[m_phoneSet->size()]; 
	
	// create a different tree for each possible left context
	for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
	
		//double dTimeBeginTree = TimeUtils::getTimeMilliseconds();
		double dTimeSubsets = 0.0;
		
		if (bPhoneRight[i] == false) {
			lexiconTree[i] = NULL;
			continue;
		}
		
		++iTreesTotal;
		unsigned int iStateId = 0;
		lexiconTree[i] = new LexiconTree;
		lexiconTree[i]->stateRoot = newState(iStateId++,NULL,NULL);
		unsigned int iLexiconSize = m_lexiconManager->getLexiconSize();			// # lexical units including alternative pronunciations
		lexiconTree[i]->iLeafRightContext = new unsigned int[iLexiconSize];
		lexiconTree[i]->leafRightContext = new LeafRightContext*[iLexiconSize];
	
		unsigned char iPhoneLeftContext = i;
		State *stateRoot = lexiconTree[i]->stateRoot;	
		
		VLeafRightContext vLeafRightContext;
		LeafRightContext leafRightContext;	
		
		// for each lexical unit
		for(VLexUnitX::iterator it = lexiconX->begin() ; it != lexiconX->end() ; ++it) {
		
			// skip these lexical units
			if (((*it)->iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit) ||
				((*it)->iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) ||
				((*it)->iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit))
			{
				lexiconTree[i]->leafRightContext[(*it)->iLexUnit] = NULL;
				continue;
			}
			
			// for each alternative pronunciation	
			for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
			
				State *state = stateRoot;
				
				if ((*jt)->vPhones.size() > iPhonesLongerLexUnit) {
					iPhonesLongerLexUnit = (*jt)->vPhones.size();
				}
				
				unsigned char iPhonePrev = iPhoneLeftContext;
				
				//m_lexiconManager->printLexUnit(*jt);
				
				// for each phone
				unsigned char iPhoneNumber = 0;
				for(vector<int>::iterator kt = (*jt)->vPhones.begin() ; kt != (*jt)->vPhones.end() ; ++kt, ++iPhoneNumber) {
					// get the next phone
					vector<int>::iterator lt = kt;
					advance(lt,1);
					
					// not the last phoneme within the lexical unit
					if (lt != (*jt)->vPhones.end()) {
					
						unsigned char iPhoneNext = *lt;
						assert(iPhoneNext <= m_phoneSet->size());
							
						for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
						
							// get the HMM-state
							unsigned char iPosition;
							if ((*jt)->vPhones.size() == 1) {
								iPosition = WITHIN_WORD_POSITION_MONOPHONE;
							} else if (iPhoneNumber == 0) {
								iPosition = WITHIN_WORD_POSITION_START;	
							} if (iPhoneNumber > 0) {
								iPosition = WITHIN_WORD_POSITION_INTERNAL;	
							}
							unsigned int iHMMState = m_hmmManager->getHMMStateIndex(&iPhonePrev,*kt,&iPhoneNext,iPosition,iState);
							// try to reuse an already existing transition (same HMM-state)
							bool bFound = false;
							for(VTransition::iterator jt = state->vTransition.begin() ; jt != state->vTransition.end() ; ++jt) {
								if ((*jt)->iSymbol == iHMMState) {
									state = (*jt)->state;
									bFound = true;
									break;
								}
							}
							// cannot reuse, create a new transition and a new state
							if (bFound == false) {
								// create the next state
								Transition *transitionAux = newTransition(iHMMState,0.0,NULL);
								transitionAux->state = newState(iStateId++,transitionAux,state);
								++iStatesTotal;
								state->vTransition.push_back(transitionAux);
								++iTransitionsTotal;
								state = transitionAux->state;
							}	
						}
					}
					// last phoneme within the lexical unit (use all the possible right contexts)
					else {
					
						State *stateCommon = state;
						
						// get the within word position (it can be "start" in the case of monophone lexical units)
						unsigned char iPosition = WITHIN_WORD_POSITION_END;
						if (iPhoneNumber == 0) {
							iPosition = WITHIN_WORD_POSITION_START;	
						}
					
						for(unsigned char iPhoneNext = 0 ; iPhoneNext < m_phoneSet->size() ; ++iPhoneNext) {
						
							if (bPhoneLeft[iPhoneNext] == false) {
								stateLeaf[iPhoneNext] = NULL;
								continue;
							}
							
							/*if ((*jt == m_lexiconManager->getLexUnitPronunciation("A")) && (i == m_phoneSet->getPhoneIndex("AY"))) {
								printf("next phone: %s counter %d\n",m_phoneSet->getStrPhone(iPhoneNext),++iCounter);
							}*/
							
							state = stateCommon;
							assert(iPhoneNext <= m_phoneSet->size());
								
							for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
							
								// get the HMM-state
								unsigned int iHMMState = m_hmmManager->getHMMStateIndex(&iPhonePrev,*kt,&iPhoneNext,iPosition,iState);
								// try to reuse an already existing transition (same HMM-state)
								bool bFound = false;
								for(VTransition::iterator mt = state->vTransition.begin() ; mt != state->vTransition.end() ; ++mt) {
									if ((*mt)->iSymbol == iHMMState) {
										state = (*mt)->state;
										bFound = true;
										break;
									}
								}
								// cannot reuse, create a new transition and a new state
								if (bFound == false) {
									// create the next state
									Transition *transitionAux = newTransition(iHMMState,0.0,NULL);
									transitionAux->state = newState(iStateId++,transitionAux,state);
									++iStatesTotal;
									state->vTransition.push_back(transitionAux);
									++iTransitionsTotal;
									state = transitionAux->state;
								}
									
								// insert the leaf state in the map
								if ((iState == (NUMBER_HMM_STATES-1)) && (iPhoneNumber == ((*jt)->vPhones.size()-1))) {
									stateLeaf[iPhoneNext] = state;
								}
							}
						}
					}
					iPhonePrev = *kt;
				}
				
				double dTimeBeginSubsets = TimeUtils::getTimeMilliseconds();
				
				// attach the right context group to each lexical unit leaf	
				for(unsigned char j = 0 ; j < m_phoneSet->size() ; ++j) {
					if (stateLeaf[j] != NULL) {
					
						//if ((*jt == m_lexiconManager->getLexUnitPronunciation("A")) && (i == m_phoneSet->getPhoneIndex("AY"))) {
						//	printf("looking group for: %s\n",m_phoneSet->getStrPhone(j));
						//}
					
						// create the subset to check
						unsigned char iOffset = 0;
						iRightContextSubset[iOffset++] = j; 
						for(unsigned char k = j+1 ; k < m_phoneSet->size() ; ++k) {
							if (stateLeaf[j] == stateLeaf[k]) {
								stateLeaf[k] = NULL;
								iRightContextSubset[iOffset++] = k; 
								//if ((*jt == m_lexiconManager->getLexUnitPronunciation("A")) && (i == m_phoneSet->getPhoneIndex("AY"))) {
								//	printf("%s same state\n",m_phoneSet->getStrPhone(k));
								//}
							}
						}
						unsigned char iSubsetElements = iOffset;
						iRightContextSubset[iOffset++] = UCHAR_MAX; 
						// does the subset exist already?
						bool bEqual = false;
						unsigned int k=0;
						for(vector<unsigned char*>::iterator it = vRightContextGroups.begin() ; it != vRightContextGroups.end() ; ++it, ++k) {
							bEqual = true;
							unsigned char l = 0;
							for( ; l < m_phoneSet->size() ; ++l) {
								if ((*it)[l] != iRightContextSubset[l]) {
									bEqual = false;
									break;
								}	
								if ((*it)[l] == UCHAR_MAX) {
									assert(iRightContextSubset[l] == UCHAR_MAX);
									assert(l > 0);
									break;
								}
							}
							if (bEqual) {
								break;
							}	
						}
						
						
						/*for( ; k < MAX_RIGHT_CONTEXT_GROUPS ; ++k) {	
							if ((*iRightContextGroups)[k] == NULL) {
								break;
							}
							bFound = true;
							unsigned char l = 0;
							for( ; l < m_phoneSet->size() ; ++l) {
								if ((*iRightContextGroups)[k][l] != iRightContextSubset[l]) {
									bFound = false;
									break;
								}	
								if ((*iRightContextGroups)[k][l] == UCHAR_MAX) {
									assert(iRightContextSubset[l] == UCHAR_MAX);
									assert(l > 0);
									break;
								}
							}
							if (bFound) {
								break;
							}	
						}*/
						// create it
						if (bEqual == false) {
							//assert(k < MAX_RIGHT_CONTEXT_GROUPS);
							//(*iRightContextGroups)[k] = new unsigned char[iSubsetElements+1];
							vRightContextGroups.push_back(new unsigned char[iSubsetElements+1]);
							//printf("new subset(%d): %d elements:",*iGroups,iSubsetElements);
							for(unsigned char l = 0 ; l < iSubsetElements+1 ; ++l) {	// copy the UCHAR_MAX too
								//(*iRightContextGroups)[k][l] = iRightContextSubset[l];
								vRightContextGroups[k][l] = iRightContextSubset[l];
								//if (l < iSubsetElements) {
									//printf(" %s",m_phoneSet->getStrPhone((*iRightContextGroups)[k][l]));
								//}
							}
							//++*iGroups;
							//printf("\n");
						} else {
							//printf("reused subset: %d elements\n",iSubsetElements);
						}
						// keep the lexical unit along with the group index
						leafRightContext.stateLeaf = stateLeaf[j];
						leafRightContext.iRightContextGroup = k;
						vLeafRightContext.push_back(leafRightContext);
						
						//if ((*jt == m_lexiconManager->getLexUnitPronunciation("A")) && (i == m_phoneSet->getPhoneIndex("AY"))) {
						//	printf("using group: %u\n",k);
						//}
					}
				}	
				//printf("done\n");
				
				double dTimeEndSubsets = TimeUtils::getTimeMilliseconds();
				dTimeSubsets += (dTimeEndSubsets-dTimeBeginSubsets);
				
				// allocate memory for the (leaf state) <-> (right context group) pairs
				assert(vLeafRightContext.size() < UCHAR_MAX);
				lexiconTree[i]->iLeafRightContext[(*jt)->iLexUnitPron] = vLeafRightContext.size();
				lexiconTree[i]->leafRightContext[(*jt)->iLexUnitPron] = new LeafRightContext[vLeafRightContext.size()];
				
				unsigned char j = 0;
				for(VLeafRightContext::iterator kt = vLeafRightContext.begin() ; kt != vLeafRightContext.end() ; ++kt, ++j) {
					lexiconTree[i]->leafRightContext[(*jt)->iLexUnitPron][j].stateLeaf = kt->stateLeaf;
					lexiconTree[i]->leafRightContext[(*jt)->iLexUnitPron][j].iRightContextGroup = kt->iRightContextGroup;
						//if ((*jt == m_lexiconManager->getLexUnitPronunciation("A")) && (i == m_phoneSet->getPhoneIndex("AY"))) {
						//	printf("final group: %d\n",kt->iRightContextGroup);
						//}
				}
				vLeafRightContext.clear();
				
				// sanity check
				int iHMMStates = (*jt)->vPhones.size()*3;
				for(unsigned char j = 0 ; j < lexiconTree[i]->iLeafRightContext[(*jt)->iLexUnitPron] ; ++j) {
					State *stateAux = lexiconTree[i]->leafRightContext[(*jt)->iLexUnitPron][j].stateLeaf;
					for(unsigned char k = 0 ; k < iHMMStates ; ++k) {
						assert(stateAux->stateParent != NULL);
						stateAux = stateAux->stateParent;
					}
					assert(stateAux->stateParent == NULL);
				}	
			}
		}
	
		// number the nodes of the tree in a post-order fashion
		iStateId = 0;
		preOrderNumbering(stateRoot,&iStateId);
		
		//double dTimeEndTree = TimeUtils::getTimeMilliseconds();
		//printf("%s tree building time: %.2f seconds, subsets: %.2f\n",m_phoneSet->getStrPhone(i),(dTimeEndTree-dTimeBeginTree)/1000.0,(dTimeSubsets/1000.0));	
		
	}
	
	// get the maximum tree depth (all the trees have the same depth so use the first one)
	for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
		if (bPhoneRight[i] == true) {
			*iTreeDepth = getTreeDepth(lexiconTree[i]->stateRoot);
			break;
		}
	}	
	
	// get the number of HMM-states
	int iHMMStates = 0;
	if (m_hmmManager->getHMMStatesDecoding(&iHMMStates) == NULL) {
		return NULL;
	}
	
	unsigned char iContextModelingOrder = m_hmmManager->getContextModelingOrderHMM();
	const char *strContextModelingOrder = Accumulator::getContextModelingOrder(iContextModelingOrder);
	
	cout << "# HMM-states:    " << iHMMStates << " (context modeling order: " 
		<< strContextModelingOrder << ")" << endl;
	cout << "# left contexts: " << iLeftContexts << ", # right contexts: " << iRightContexts << endl;
	cout << "# states:        " << iStatesTotal << endl;
   cout << "# transitions:   " << iTransitionsTotal << endl;
	cout << "max tree depth:  " << *iTreeDepth << endl;
	cout << "# right context subsets: " << vRightContextGroups.size() << endl;
   
   delete [] stateLeaf;
   delete [] bPhoneLeft;
   delete [] bPhoneRight;
   
   return lexiconTree;
}

// pre-order state numbering
void WFSABuilder::preOrderNumbering(State *state, unsigned int *iId) {

	state->iId = (*iId)++;
	
	for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		preOrderNumbering((*it)->state,iId);
	}
}

// return the maximum tree depth
int WFSABuilder::getTreeDepth(State *state, int iDepth) {

	int iDepthMax = iDepth;
	for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		int iDepthAux = getTreeDepth((*it)->state,iDepth+1);
		if (iDepthAux > iDepthMax) {
			iDepthMax = iDepthAux; 
		}
	}
	
	return iDepthMax;
}


// make each state in the acceptor (WFSA) accessible only by transitions with the same input symbol (for more efficient decoding)
void WFSABuilder::equalizeInput(State *stateInitial, unsigned int iStates, unsigned int iTransitions, unsigned int *iStatesFinal, unsigned int *iTransitionsFinal) {

	//double dTimeBegin = TimeUtils::getTimeMilliseconds();

	VState vStates;							// states yet to be processed
	int iStatesCreated = 0;					// keeps the count of created states
	int iTransitionsCreated = 0;			// keeps the count of created transitions
	unsigned int iStateId = iStates;
	
	// data structures needed to keep track of seen states and pairs (symbol,state), each symbol can only be reached by one kind of state	
	bool *bStateSeen = new bool[iStates*2];
	unsigned int *iStateSymbol = new unsigned int[iStates*2];
	for(unsigned int i=0 ; i<iStates*2 ; ++i) {
		bStateSeen[i] = false;
		iStateSymbol[i] = UINT_MAX;
	}
	unsigned int iStatesSeen = 0;
	unsigned int iTransitionsSeen = 0;
	
	// (1) each of the transitions should go either to a state with the same incoming input symbol associated or to a state that is not yet in the map
	
	// keep the original number of states and transitions
	unsigned int iStatesOriginal = iStates;
	unsigned int iTransitionsOriginal = iTransitions;
	
	// insert the initial state
	vStates.push_back(stateInitial);		
	bStateSeen[stateInitial->iId] = true;
	++iStatesSeen;
	// process all the states in the machine
	while(vStates.empty() == false) {
		
		State *state = vStates.back();
		vStates.pop_back();
		
		for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		
			if (((*it)->iSymbol & LEX_UNIT_TRANSITION) || ((*it)->iSymbol & EPSILON_TRANSITION)) {
				if (bStateSeen[(*it)->state->iId] == false) {
					vStates.push_back((*it)->state);
					bStateSeen[(*it)->state->iId] = true;
					++iStatesSeen;
				}
				continue;
			}
		
			// check if the destination state is in the map
			// destination state not seen yet: attach the symbol to it
			if (iStateSymbol[(*it)->state->iId] == UINT_MAX) {
				iStateSymbol[(*it)->state->iId] = (*it)->iSymbol;
			} 
			else {
				// if the symbol does not match create a new state (unless it was just created, the transducer might not be deterministic)
				if ((*it)->iSymbol != iStateSymbol[(*it)->state->iId]) {	
					State *stateNext = newState(iStateId++,NULL,NULL);
					//iStateSymbol[(*it)->state->iId] = (*it)->iSymbol;
					iStateSymbol[iStateId] = (*it)->iSymbol;
					++iStatesCreated;
					// replicate outgoing transitions in the original state
					for(VTransition::iterator kt = (*it)->state->vTransition.begin() ; kt != (*it)->state->vTransition.end() ; ++kt) {
						stateNext->vTransition.push_back(newTransition((*kt)->iSymbol,(*kt)->fWeight,(*kt)->state));
						++iTransitionsCreated;
					}
					// move the transitions from the original state to the new state
					bool bChanged = false;
					for(VTransition::iterator kt = it ; kt != state->vTransition.end() ; ++kt) {
						if (((*kt)->iSymbol == (*it)->iSymbol) && ((*kt)->state == (*it)->state)) {
							(*kt)->state = stateNext;
							bChanged = true;
						}
					}
					assert(bChanged);
					// insert the created state in the queue
					vStates.push_back(stateNext);
					bStateSeen[stateNext->iId] = true;
					++iStatesSeen;					
				}
			}
			// if the destination state is not processed	yet add it to the queue
			if (bStateSeen[(*it)->state->iId] == false) {
				vStates.push_back((*it)->state);
				// mark the state as seen
				bStateSeen[(*it)->state->iId] = true;
				++iStatesSeen;				
			}
		}
	}
	
	delete [] bStateSeen;
	delete [] iStateSymbol;
	
	assert(iStatesSeen == (iStates+iStatesCreated));
	
	//printf("# seen states: %d\n",iStatesSeen);
	//printf("# added states: %d\n",iStatesCreated);
	//printf("# added transitions: %d\n",iTransitionsCreated);
	
	assert(vStates.empty());
		
	float fStatesIncrease = 100*((float)(iStatesCreated)/((float)iStatesOriginal));
	float fTransitionsIncrease = 100*((float)(iTransitionsCreated)/((float)iTransitionsOriginal));
	//printf("Increase: #states: %.2f%%  #transitions: %.2f%%\n",fStatesIncrease,fTransitionsIncrease);
	
	// sanity check: make sure all states are reachable by only one kind of input symbol
	
	// data structures needed to keep track of seen states and pairs (symbol,state), each symbol can only be reached by one kind of state	
	bStateSeen = new bool[iStates*2];
	iStateSymbol = new unsigned int[iStates*2];
	for(unsigned int i=0 ; i<iStates*2 ; ++i) {
		bStateSeen[i] = false;
		iStateSymbol[i] = UINT_MAX;
	}
	iStatesSeen = 0;	
	iTransitionsSeen = 0;
	
	vStates.push_back(stateInitial);	
	bStateSeen[stateInitial->iId] = true;
	iTransitionsSeen += stateInitial->vTransition.size();
	
	while(vStates.empty() == false) {
		
		State *state = vStates.back();
		vStates.pop_back();
		
		for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		
			// look for two lexical unit symbols in a row
			if ((*it)->iSymbol & LEX_UNIT_TRANSITION) {
				for(VTransition::iterator jt = (*it)->state->vTransition.begin() ; jt != (*it)->state->vTransition.end() ; ++jt) {
					if ((*jt)->iSymbol & LEX_UNIT_TRANSITION) {
						if ((int)((*jt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) != m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
							int iLexUnit1 = ((*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
							int iLexUnit2 = ((*jt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
							printf("two lex units in a row!! %s %s\n",m_lexiconManager->getStrLexUnitPron(iLexUnit1),m_lexiconManager->getStrLexUnitPron(iLexUnit2));
						}
					}
				}			
			}
		
			// check if the destination state is in the map
			// destination state not seen yet: attach the symbol to it
			if (iStateSymbol[(*it)->state->iId] == UINT_MAX) {	
				iStateSymbol[(*it)->state->iId] = (*it)->iSymbol;
			} else {
				// no state can be reached by different leaf-transitions
				if (((*it)->iSymbol < LEX_UNIT_TRANSITION) && ((*it)->iSymbol >= 0) && 
					(iStateSymbol[(*it)->state->iId] < LEX_UNIT_TRANSITION) && (iStateSymbol[(*it)->state->iId] >= 0)) {
					if ((*it)->iSymbol != iStateSymbol[(*it)->state->iId]) {
						printf("%u != %u\n",(*it)->iSymbol,iStateSymbol[(*it)->state->iId]);
					}
					assert((*it)->iSymbol == iStateSymbol[(*it)->state->iId]);
				}
			}
			// if the destination state is not processed	yet add it to the queue
			if (bStateSeen[(*it)->state->iId] == false) {
				vStates.push_back((*it)->state);
				// mark the state as seen
				bStateSeen[(*it)->state->iId] = true;
				iTransitionsSeen += (*it)->state->vTransition.size();	
			}
		}	
	}	
	
	delete [] bStateSeen;
	delete [] iStateSymbol;	
	
	//printf("# seen states: %d\n",iStatesSeen);
	//printf("# seen transitions: %u\n",iTransitionsSeen);
	
	*iStatesFinal = iStates+iStatesCreated;
	*iTransitionsFinal = iTransitions+iTransitionsCreated;
	
	assert(vStates.empty());
	
	//double dTimeEnd = TimeUtils::getTimeMilliseconds();
	
	printf("-- input equalization ----------------------\n");
	printf(" # states original:      %9u\n",iStatesOriginal);
	printf(" # transitions original: %9u\n",iTransitionsOriginal);
	printf(" # states added:         %9u (%.2f%%)\n",iStatesCreated,fStatesIncrease);
	printf(" # transitions added:    %9u (%.2f%%)\n",iTransitionsCreated,fTransitionsIncrease);
	printf("--------------------------------------------\n");
}

// optimize the acceptor for decoding
WFSAcceptor *WFSABuilder::optimize(State *stateInitial, unsigned int iStates, unsigned int iTransitions, 
	State *stateFinalG) {
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	// allocate memory for the states and transitions
	// # states = # states in original machine + one extra state to keep the end of the array of transitions of the previous state
	StateX *states = new StateX[iStates+1];
	StateX *stateNextAvailable = states;
	
	// # transitions = # transition in original machine + one fake transition for each final state (to keep its weight)
	TransitionX *transitions = new TransitionX[iTransitions+1];
	TransitionX *transitionNextAvailable = transitions;
	
	VStatePair vStatePair;
	//MStateStateX mStateSeen;
	StateX **stateX = new StateX*[iStates+1];
	for(unsigned int i=0 ; i<iStates+1 ; ++i) {
		stateX[i] = NULL;
	}
	
	// get the number of lexical units
	float fMax = -10000000000.0;
	float fMin = 10000000000.0;
	
	int iTransitionsSeen = 0;
	
	// insert the first state in the queue
	StatePair statePairInitial;
	statePairInitial.state = stateInitial;
	statePairInitial.stateX = stateNextAvailable;
	vStatePair.push_back(statePairInitial);
	//mStateSeen.insert(MStateStateX::value_type(stateInitial,stateNextAvailable));
	stateX[stateInitial->iId] = stateNextAvailable;
	*stateNextAvailable = transitionNextAvailable;
	transitionNextAvailable += stateInitial->vTransition.size();
	++stateNextAvailable;
	iTransitionsSeen += stateInitial->vTransition.size();
	
	// process all the states in the transducer
	while(vStatePair.empty() == false) {
		
		StatePair statePair = vStatePair.back();
		vStatePair.pop_back();
		
		// is a final state? (only if destination state == NULL)
		if (statePair.state->vTransition.empty()) {
			// create a transition to keep the weight
			(*statePair.stateX)->iSymbol = FAKE_TRANSITION;
			//(*statePair.stateX)->fWeight = (- m_mStateWeightFinal[statePair.state]);
			(*statePair.stateX)->fWeight = 0.0;
			(*statePair.stateX)->state = NULL;
			continue;
		} 
		
		TransitionX *transitionNextAvailableBase = *statePair.stateX;
		
		// create the transitions in the acceptor
		for(VTransition::iterator it = statePair.state->vTransition.begin() ; it != statePair.state->vTransition.end() ; ++it) {
			// epsilon transition	
			if ((*it)->iSymbol == EPSILON_TRANSITION) {
				transitionNextAvailableBase->iSymbol = EPSILON_TRANSITION;	
			} 
			// lex-unit transition (lexical unit)
			else if ((*it)->iSymbol & LEX_UNIT_TRANSITION) {
				//printf("weight: %.2f\n",(*it)->fWeight);
				if ((*it)->fWeight > fMax) {
					fMax = (*it)->fWeight;
				} else if ((*it)->fWeight < fMin) {
					fMin = (*it)->fWeight;
				}
				int iLexUnit = (int)(*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				assert(iLexUnit < m_iLexUnitsTotal);
				// replace the start and end of sentence lexical unit by an epsilon 
				/*if ((iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) || (iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) {
					transitionNextAvailableBase->iSymbol = EPSILON_TRANSITION;
				} else {*/
					transitionNextAvailableBase->iSymbol = (*it)->iSymbol;
				//}
			}
			// leaf transition (HMM-state)
			else {
				assert((int)(*it)->iSymbol < m_iHMMStatesTotal);
				transitionNextAvailableBase->iSymbol = (*it)->iSymbol;	
			}			
			//transitionNextAvailableBase->fWeight = (-(*it)->fWeight);			// negate the weight
			transitionNextAvailableBase->fWeight = (*it)->fWeight;
			
			// is the destination state already seen?
			//MStateStateX::iterator jt = mStateSeen.find((*it)->state);
			//if (jt == mStateSeen.end()) {
			if (stateX[(*it)->state->iId] == NULL) {
			
				// add the state to the queue
				StatePair statePairAux;
				statePairAux.state = (*it)->state;
				statePairAux.stateX = stateNextAvailable;
				//(*(statePairAux.stateX))->iSymbol = 12345;
				vStatePair.push_back(statePairAux);	
				// destination state
				transitionNextAvailableBase->state = stateNextAvailable;
				*stateNextAvailable = transitionNextAvailable;
				//mStateSeen.insert(MStateStateX::value_type((*it)->state,stateNextAvailable));
				stateX[(*it)->state->iId] = stateNextAvailable;
				++stateNextAvailable;
				transitionNextAvailable += max((int)((*it)->state->vTransition.size()),1);		// we will create a fake transition for each final state
			} else {
				// destination state
				//transitionNextAvailableBase->state = jt->second;	
				transitionNextAvailableBase->state = stateX[(*it)->state->iId];
			}
			++transitionNextAvailableBase;
			++iTransitionsSeen;
		}	
	}
	
	// keep the end of the array of transitions for the last state
	states[iStates] = transitionNextAvailable;
	
	// sanity checks
	//printf("transitions seen: %d\n",iTransitionsSeen);
	//int iAux2 = mStateSeen.size();
	//printf("%u %u\n",mStateSeen.size(),iStates);
	//assert(mStateSeen.size() == iStates);
	//printf("%u %u\n",stateNextAvailable-states,iStates);
	assert(stateNextAvailable-states == (int)iStates);
	//printf("%u %u\n",transitionNextAvailable-transitions,iTransitions+1);
	assert(transitionNextAvailable-transitions == (int)iTransitions+1);	
	
	for(unsigned int i=0; i < iTransitions+1 ; ++i) {
		assert(((int)transitions[i].iSymbol < m_iHMMStatesTotal) || 
				(transitions[i].iSymbol == EPSILON_TRANSITION) || 
				(transitions[i].iSymbol == FAKE_TRANSITION) || 
				(transitions[i].iSymbol & LEX_UNIT_TRANSITION));
	}
	/*
	// traverse the whole machine
	VStateX vStateX;
	MStateX mStateX;
	
	vStateX.push_back(&states[0]);
	mStateX.insert(MStateX::value_type(&states[0],false));
	unsigned int iFinalStatesFound = 0;
	unsigned int iStatesSeen = 0;
	
	while(vStateX.empty() == false) {
	
		StateX *stateX = vStateX.back();
		vStateX.pop_back();
		++iStatesSeen;
	
		TransitionX *transition =  *stateX;
		TransitionX *transitionEnd = *(stateX+1);
	
		while(transition != transitionEnd) {
		
			if (transition->state == NULL) {
				++iFinalStatesFound;
				++transition;
				continue;
			}
		
			if (mStateX.find(transition->state) == mStateX.end()) {
				vStateX.push_back(transition->state);
				bool bLexUnitInput = false;
				if (transition->iSymbol & LEX_UNIT_TRANSITION) {
					bLexUnitInput = true;
					assert(mStateX[stateX] == false);
				} else {
					//printf("%u\n",transition->iInput);
				}
				mStateX.insert(MStateX::value_type(transition->state,bLexUnitInput));
			} 
			
			++transition;
		}
	}
	
	printf("# states visited in optimized acceptor: %d (%u final states found)\n",mStateX.size(),iFinalStatesFound);
	*/
	// the initial state is the first state
	StateX *stateInitialOptimized = &states[0];	
	
	delete [] stateX;	
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	
	printf("-- acceptor optimization (reducing memory use) ----\n");
	printf(" max weight: %8.2f\n",fMax);
	printf(" min weight: %8.2f\n",fMin);
	printf(" time:       %8.2f seconds\n",(dTimeEnd-dTimeBegin)/1000.0);
	printf("--------------------------------------------------\n");	

	return new WFSAcceptor(states,iStates,stateInitialOptimized,
		transitions,iTransitions+1,m_iHMMStatesTotal,m_iLexUnitsTotal);
}


// detroy a lexicon tree
void WFSABuilder::destroyTree(State *state) {

}





// build a word-internal context dependent decoding network or a context-independent one
WFSAcceptor *WFSABuilder::buildWordInternal() {

	unsigned char iContextSize = m_hmmManager->getContextSizeHMM();
	assert(iContextSize == 0);
	
	// (1) create the language model graph
	unsigned int iStatesG = 0;
	unsigned int iTransitionsG = 0;
	//VState vStateFinalG;
	State *stateFinalG = NULL;
	State *states = NULL;
	State *stateRootG = NULL;
	if (buildG(&states,&iStatesG,&iTransitionsG,&stateRootG,&stateFinalG) == false) {
		return NULL;
	}
	// create a list with the states in G sorted by topological order with respect to epsilon transitions

	// (2) build the lexical tree (only one tree needs to be built)
	map<LexUnit*,State*> mLexUnitState;
	unsigned int iTreeDepth = 0;
	State *stateRoot = buildWordInternalLexiconTree(mLexUnitState,&iTreeDepth);
	if (stateRoot == NULL) {
		return NULL;
	}
	//printStateTransitions(stateRoot);
	
	// (3) do the actual incremental composition
	unsigned int iStateIdFinal = iStatesG;
	LLeafToProcess lLeafToProcess;													// leaf queue
	State **stateWaiting = new State*[iTreeDepth];								// waiting states
	VTransition *transitionInsertion = new VTransition[iTreeDepth];		// transitions ready to be inserted
	for(unsigned int i=0 ; i < iTreeDepth ; ++i) {
		stateWaiting[i] = NULL;
	}
	
	// create an array to mark the already processed states
	bool *bStateProcessed = new bool[iStatesG];
	for(unsigned int i=0 ; i < iStatesG ; ++i) {
		bStateProcessed[i] = false;
	}
	
	// traverse G and apply incremental composition to each state
	MStateState mStateState;
	unsigned int iTransitionsProcessed = 0;
	unsigned int iStatesProcessed = 0;
	unsigned int iTransitionsCreated = 0;
	unsigned int iStatesCreated = 0;
	unsigned int iTransitionsRemovedFromG = 0;
	unsigned int iStatesReused = 0;
	LState lState;
	lState.push_back(stateRootG);
	bStateProcessed[stateRootG->iId] = true;
	while(lState.empty() == false) {
		
		// get the next state to process
		State *stateGFrom = lState.front();
		lState.pop_front();
		
		//printf("processing G state:\n");
	
		// process the state	in G by processing all the transitions that come from it
		bool bLexUnitTransitions = false;		// whether actual lexical unit transitions (transitions that are actual words) are seen
		
		//printf("G state to be processed has: %d transitions\n",stateGFrom->vTransition.size());
		
		// (3.1) for each lexical-unit transition coming from the state, find the corresponding leaf in the lexical tree
		VTransition vTransitionSurvive;
		for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
		
			// insert G destination states into the queue (if not already processed)
			if (bStateProcessed[(*it)->state->iId] == false) {
				lState.push_back((*it)->state);
				bStateProcessed[(*it)->state->iId] = true;
			}	
			// epsilon transition: keep it in the final tree 
			if ((*it)->iSymbol == EPSILON_TRANSITION) {
				vTransitionSurvive.push_back(*it);
				//printf("epsilon processed\n");
			} 
			// beginning/end of sentence: keep it in the final tree
			else if (((int)(*it)->iSymbol == (m_lexiconManager->m_lexUnitBegSentence->iLexUnit|LEX_UNIT_TRANSITION)) || 
				(((int)(*it)->iSymbol == (m_lexiconManager->m_lexUnitEndSentence->iLexUnit|LEX_UNIT_TRANSITION)))) {
				vTransitionSurvive.push_back(*it);
				//printf("%s processed\n",m_lexiconManager->getStrLexUnit((*it)->iSymbol));	
			}	
			// non-epsilon transition: put leaves from each pronuciation of the lexical unit in the queue
			else {
				bLexUnitTransitions = true;
				LexUnitX *lexUnitX = m_lexiconManager->getLexUnit((*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
				for(VLexUnit::iterator jt = lexUnitX->vLexUnitPronunciations.begin() ; jt != lexUnitX->vLexUnitPronunciations.end() ; ++jt) {
					LeafToProcess *leaf = new LeafToProcess;
					leaf->lexUnit = *jt;
					leaf->iHMMStates = (*jt)->vPhones.size()*NUMBER_HMM_STATES;
					leaf->stateLeaf = mLexUnitState[*jt]; 
					leaf->iRightContextGroup = UCHAR_MAX;
					leaf->stateGTo = (*it)->state;
					leaf->fWeight = (*it)->fWeight;
					lLeafToProcess.push_back(leaf);	
					//m_lexiconManager->printLexUnit(*jt);
				}
				delete *it;
			}	
		}
		
		iTransitionsRemovedFromG += stateGFrom->vTransition.size()-vTransitionSurvive.size();
		// remove old G transitions
		stateGFrom->vTransition.clear();
		// add the G transitions that survive
		for(VTransition::iterator it = vTransitionSurvive.begin() ; it != vTransitionSurvive.end() ; ++it) {
			stateGFrom->vTransition.push_back(*it);
		}
		
		// sort the queue by state number
		lLeafToProcess.sort(WFSABuilder::compareStateNumber);
		
		// process leaves in the queue
		while(lLeafToProcess.empty() == false) {
		
			LeafToProcess *leaf = lLeafToProcess.front();
			lLeafToProcess.pop_front();
			
			int iHMMStates = leaf->iHMMStates;				
			State *state = leaf->stateLeaf;					// this state goes to the lexical unit leaf
			unsigned int iSymbol = ((unsigned int)leaf->lexUnit->iLexUnitPron)|LEX_UNIT_TRANSITION;
			
			// process states from the state that goes to the leaf (lexical unit) to the root (that goes to the first HMM-state)
			assert(iHMMStates < (int)iTreeDepth);
			for(int i=iHMMStates ; i >= 0 ; --i) {
			
				//printStateTransitions(state);
			
				// empty position: insert the state into the buffer
				if (stateWaiting[i] == NULL) {
					stateWaiting[i] = state;
					transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
					assert((i%NUMBER_HMM_STATES == (int)iSymbol%NUMBER_HMM_STATES) || (iSymbol >= LEX_UNIT_TRANSITION));
					checkTransitionSymbolCorrectness(iSymbol);
					++iTransitionsCreated;
					// if the state goes to a G node we already know the destination state in the final graph
					if (state == leaf->stateLeaf) {
						transitionInsertion[i].back()->state = leaf->stateGTo;
						transitionInsertion[i].back()->fWeight = leaf->fWeight;
					}
				}
				// different state: move the original state along with its transitions to the final tree
				else if (stateWaiting[i] != state) {
					// (1) the state in the buffer becomes final
					State *stateFinal = newState(iStateIdFinal++,NULL,NULL);
					++iStatesCreated;
					assert(i > 0);		// the root state stays in the buffer until the end
					assert(transitionInsertion[i].empty() == false);
					for(VTransition::iterator jt = transitionInsertion[i].begin() ; jt != transitionInsertion[i].end() ; ++jt) {
						stateFinal->vTransition.push_back(*jt);
						// if no destination state: look for the final state that does have destination state
						if ((*jt)->state == NULL) {	
							assert(i == iHMMStates);
							bool bFound = false;
							for(int j = iTreeDepth-1 ; j > i ; --j) {
								if (stateWaiting[j] != NULL) {
									while(i != j) {
										State *stateFinal2 = newState(iStateIdFinal++,NULL,NULL);
										++iStatesCreated;
										for(VTransition::iterator kt = transitionInsertion[j].begin() ; kt != transitionInsertion[j].end() ; ++kt) {
											stateFinal2->vTransition.push_back(*kt);
											assert((*kt)->state != NULL);
										}
										stateWaiting[j] = NULL;
										transitionInsertion[j].clear();
										// apply weight pushing
										float fWeightRest2 = applyWeightPushing(stateFinal2);	
										// try to reuse the state
										MStateState::iterator ft = mStateState.find(stateFinal2);
										if (ft != mStateState.end()) {
											for(VTransition::iterator gt = stateFinal2->vTransition.begin() ; gt != stateFinal2->vTransition.end() ; ++gt) {
												delete *gt;
												--iTransitionsCreated;
											}
											delete stateFinal2;
											iStateIdFinal--;
											--iStatesCreated;		
											++iStatesReused;
											stateFinal2 = ft->second;
											if (iStatesReused % 10000 == 0) {
												printf("(2) reused: %u total: %u (trans: %u)\n",iStatesReused,iStatesCreated,iTransitionsCreated);
											}
										} else {
											mStateState.insert(MStateState::value_type(stateFinal2,stateFinal2));
										}
										assert(transitionInsertion[j-1].back()->state == NULL);
										transitionInsertion[j-1].back()->state = stateFinal2;
										transitionInsertion[j-1].back()->fWeight = fWeightRest2;
										--j;
									}
									bFound = true;
								}	
							}
							assert(bFound);
						}
						assert((*jt)->state != NULL);
					}
					transitionInsertion[i].clear();
					// apply weight pushing
					float fWeightRest = applyWeightPushing(stateFinal);	
					// (2) insert the state in the map for later reuse
					MStateState::iterator ft = mStateState.find(stateFinal);
					if (ft != mStateState.end()) {
						for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
							delete *gt;
							--iTransitionsCreated;
						}
						delete stateFinal;
						iStateIdFinal--;
						--iStatesCreated;	
						++iStatesReused;
						stateFinal = ft->second;	
						if (iStatesReused % 10000 == 0) {
							printf("(3) reused: %u total: %u\n",iStatesReused,iStatesCreated);
						}
					} else {
						mStateState.insert(MStateState::value_type(stateFinal,stateFinal));	
					}
					assert(transitionInsertion[i-1].back()->state == NULL);
					transitionInsertion[i-1].back()->state = stateFinal;
					transitionInsertion[i-1].back()->fWeight = fWeightRest;
					// (3) the new state is inserted into the buffer along with its transition
					stateWaiting[i] = state;
					transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
					assert((i%NUMBER_HMM_STATES == (int)iSymbol%NUMBER_HMM_STATES) || (iSymbol >= LEX_UNIT_TRANSITION));
					checkTransitionSymbolCorrectness(iSymbol);
					++iTransitionsCreated;
					// if the state goes to a G node we already know the destination state in the final graph
					if (state == leaf->stateLeaf) {
						transitionInsertion[i].back()->fWeight = leaf->fWeight;
						transitionInsertion[i].back()->state = leaf->stateGTo;
					} 
				} 
				// equal state: stop the right to left traversal
				else {
					// create a new transition
					transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
					assert((i%NUMBER_HMM_STATES == (int)iSymbol%NUMBER_HMM_STATES) || (iSymbol >= LEX_UNIT_TRANSITION));
					checkTransitionSymbolCorrectness(iSymbol);
					++iTransitionsCreated;
					// if the state goes to a G node we already know the destination state in the final graph
					if (state == leaf->stateLeaf) {
						transitionInsertion[i].back()->fWeight = leaf->fWeight;
						transitionInsertion[i].back()->state = leaf->stateGTo;
					}	
					break;
				}
				
				if (i == 0) {
					assert(state->stateParent == NULL);
					break;	
				} 
				
				// move to the parent state
				iSymbol = state->transitionBW->iSymbol;
				state = state->stateParent;
			}
				
			delete leaf;
		}
		
		// insert all the waiting states in the final graph (if any)
		if (bLexUnitTransitions) {
			//assert(stateWaiting[iTreeDepth-1] == NULL);
			assert(stateWaiting[0] != NULL);							// the root always has to be in the buffer	
			for(int i=iTreeDepth-1 ; i >= 0 ; --i) {
				if (stateWaiting[i] != NULL) {	
					while(i >= 0) {
						State *stateFinal = NULL;
						if (i > 0) {
							stateFinal = newState(iStateIdFinal++,NULL,NULL);
							++iStatesCreated;	
						} else {
							// (i == 0) connect the root to the final graph
							stateFinal = stateGFrom;	
						}
						for(VTransition::iterator kt = transitionInsertion[i].begin() ; kt != transitionInsertion[i].end() ; ++kt) {
							stateFinal->vTransition.push_back(*kt);
							assert((*kt)->state != NULL);
						}
						float fWeightRest = 0.0;
						if (i > 0) {
							// apply weight pushing
							fWeightRest = applyWeightPushing(stateFinal);
						}
						
						// insert the state in the map for later reuse
						MStateState::iterator ft = mStateState.find(stateFinal);
						if (ft != mStateState.end()) {
							for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
								delete *gt;
								--iTransitionsCreated;
							}
							delete stateFinal;
							iStateIdFinal--;
							--iStatesCreated;	
							++iStatesReused;
							stateFinal = ft->second;	
							if (iStatesReused % 10000 == 0) {
								printf("(4) reused: %u total: %u\n",iStatesReused,iStatesCreated);
							}
						} else {
							mStateState.insert(MStateState::value_type(stateFinal,stateFinal));
						}
						if (i > 0) {
							assert(transitionInsertion[i-1].back()->state == NULL);
							transitionInsertion[i-1].back()->state = stateFinal;	
							transitionInsertion[i-1].back()->fWeight = fWeightRest;	
						}
						//printStateTransitions(stateFinal);
						stateWaiting[i] = NULL;
						transitionInsertion[i].clear();
						--i;
					}
					break;
				}
			}
			// sanity check
			for(unsigned int i=0 ; i < iTreeDepth-1 ; ++i) {
				assert(stateWaiting[i] == NULL);
				assert(transitionInsertion[i].empty());
			}
		}	
		
		iTransitionsProcessed += stateGFrom->vTransition.size();
		iStatesProcessed++;
			
		//printStateTransitions(stateGFrom);
		if (iStatesProcessed % 10000 == 0) {
			printf("G node processed! (%u states processed, %u transitions processed)\n",iStatesProcessed,iTransitionsProcessed);
		}
		//printf("(%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
	}
	// sanity check
	for(unsigned int i=0 ; i < iStatesG ; ++i) {
		assert(bStateProcessed[i]);
	}
	delete [] bStateProcessed;
	delete [] stateWaiting;
	delete [] transitionInsertion;
	
	// destroy the lexical tree
	destroyTree(stateRoot);
	
	printf("# states reused: %u\n",iStatesReused);
	printf("ending!!! (%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
	
	// sanity checks: 
	// (1) make sure there are not unconnected states
	// (2) make sure transition symbols are correct
	int iHMMStatesTotal = -1;
	m_hmmManager->getHMMStates(&iHMMStatesTotal);
	unsigned int iTransitionsSeen = 0;
	iStatesProcessed = 0;
	bStateProcessed = new bool[iStateIdFinal];
	for(unsigned int i=0 ; i < iStateIdFinal ; ++i) {
		bStateProcessed[i] = false;
	}
	assert(lState.empty());
	lState.push_back(stateRootG);
	bStateProcessed[stateRootG->iId] = true;
	++iStatesProcessed;
	while(lState.empty() == false) {
		
		// get the next state to process
		State *stateGFrom = lState.front();
		lState.pop_front();
		
		// (3.1) for each lexical-unit transition coming from the state, find the corresponding leaf in the lexical tree
		for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
		
			checkTransitionSymbolCorrectness((*it)->iSymbol);
		
			++iTransitionsSeen;
			if (bStateProcessed[(*it)->state->iId] == false) {
				lState.push_back((*it)->state);
				bStateProcessed[(*it)->state->iId] = true;
				++iStatesProcessed;
			}
		}	
	}
	for(unsigned int i=0 ; i < iStateIdFinal ; ++i) {
		assert(bStateProcessed[i]);
	}
	delete [] bStateProcessed;
	printf("transitions: %u seen, %u created, %u removed\n",iTransitionsSeen,iTransitionsCreated,iTransitionsRemovedFromG);
	assert(iTransitionsSeen == (iTransitionsCreated+iTransitionsG-iTransitionsRemovedFromG));
	
	unsigned int iStatesTotal = iStatesG+iStatesCreated;
	unsigned int iTransitionsTotal = iTransitionsG+iTransitionsCreated-iTransitionsRemovedFromG;
	float fNetworkSize = iStatesTotal*4+iTransitionsTotal*12;
	printf("total states: %u, total transitions: %u\n",iStatesTotal,iTransitionsTotal);
	printf("Expected network size: %f MB\n",fNetworkSize/(1024.0*1024.0));
	
	//exit(-1);
	
	// (4) equalize the acceptor input
	unsigned int iStatesTotalFinal = 0;
	unsigned int iTransitionsTotalFinal = 0;
	equalizeInput(stateRootG,iStatesTotal,iTransitionsTotal,&iStatesTotalFinal,&iTransitionsTotalFinal);
	
	// (5) transform the resulting graph to an acceptor that can be used for decoding
	WFSAcceptor *wfsAcceptor = optimize(stateRootG,iStatesTotalFinal,iTransitionsTotalFinal,stateFinalG);

	return wfsAcceptor;
}


// build a cross-word context dependent decoding network 
WFSAcceptor *WFSABuilder::buildCrossWord() {

	unsigned char iContextSize = m_hmmManager->getContextSizeHMM();
	assert(iContextSize == 1);
	
	// (1) create the language model graph
		
	unsigned int iStatesG = 0;
	unsigned int iTransitionsG = 0;
	State *stateFinalG = NULL;
	State *statesG = NULL;			
	// states in G will be topologically ordered by reverse epsilon transitions
	State *stateRootG = NULL;
	if (buildG(&statesG,&iStatesG,&iTransitionsG,&stateRootG,&stateFinalG) == false) {
		return NULL;
	}	
	
	// (2) build the lexical trees (one tree per possible left context)
	double dTimeBeginL = TimeUtils::getTimeMilliseconds();
	unsigned int iTreeDepth = 0;	
	//unsigned char **iRightContextGroups = NULL;
	//unsigned int iGroups = 0;
	vector<unsigned char*> vRightContextGroups;
	LexiconTree **lexiconTree = buildCrossWordLexiconTrees(&iTreeDepth,vRightContextGroups);
	if (lexiconTree == NULL) {
		return NULL;
	}
	double dTimeEndL = TimeUtils::getTimeMilliseconds();
	printf("L building time: %.2f seconds\n",(dTimeEndL-dTimeBeginL)/1000.0);
	printf("-----------------------------------------------\n");
	//printStateTransitions(stateRoot); 
	
	// (3) do the actual incremental composition
	unsigned int iStateIdFinal = 0;		// G states wont be part of the final graph
	LLeafToProcess lLeafToProcess;													// leaf queue
	State **stateWaiting = new State*[iTreeDepth];								// waiting states
	VTransition *transitionInsertion = new VTransition[iTreeDepth];		// transitions ready to be inserted
	for(unsigned int i=0 ; i < iTreeDepth ; ++i) {
		stateWaiting[i] = NULL;
	}

	// maps a pair (state in G, left context) to a root state in L
	MGStateLRoot mGStateLRoot;

	// maps a leaf in the unconnected graph to its right context
	//MLeafTransitionRightContext mLeafRightContextTemp;
	//MLeafTransitionVRightContext mLeafRightContext;
	VTransition vTransitionLexUnit;
	
	// (1) traverse G and apply incremental composition to each state
	MStateState mStateState;
	unsigned int iTransitionsProcessed = 0;
	unsigned int iStatesProcessed = 0;
	unsigned int iLeavesProcessed = 0;
	unsigned int iTransitionsCreated = 0;
	unsigned int iStatesCreated = 0;
	unsigned int iStatesCreatedTemporal = 0;		// temporal states: are not part of the final graph and do not have ID
	unsigned int iRootsCreatedFinalGraph = 0;
	unsigned int iStatesReused = 0;	
	unsigned int iTransitionsAuxiliarEpsilon = 0;
	// create the first state in the final graph
	State *stateRoot = newState(iStateIdFinal++,NULL,NULL);
	++iStatesCreated;	
	
	// create a definitive state that replaces the final state in G (all states in G are removed and are not part of the final graph)
	State *stateFinalGDefinitive = NULL;
	if (stateFinalG != NULL) {
		stateFinalGDefinitive = newState(iStateIdFinal++,NULL,NULL);
		++iStatesCreated;
	}	
	
	// get the left context information
	unsigned char **iLeftContext = getLeftContextGNodes(statesG,iStatesG,stateFinalG->iId);
	
	for(unsigned int g = 0 ; g < iStatesG ; ++g) {
		
		// get the next state to process
		State *stateGFrom = &statesG[g];
		Transition *transitionEpsilon = NULL;		// at most one epsilon transition coming from each state
		Transition *transitionFinal = NULL;
		
		//printf("processing G state (id= %u, transitions= %u): (%u states processed already)\n",stateGFrom->iId,iStatesProcessed,stateGFrom->vTransition.size());	
		
		// for each posible left context	
		for(unsigned char l=0 ; iLeftContext[stateGFrom->iId][l] != UCHAR_MAX ; ++l) {
		
			// get the left context
			unsigned char iLeftPhone = iLeftContext[stateGFrom->iId][l];
			assert(lexiconTree[iLeftPhone] != NULL);
			
			State *stateRootL = newState(UINT_MAX,NULL,NULL);			// this state wont be part of the final graph
			mGStateLRoot.insert(MGStateLRoot::value_type(pair<State*,unsigned char>(stateGFrom,iLeftPhone),stateRootL));	
			iStatesCreatedTemporal++;
			++iRootsCreatedFinalGraph;
			
			// process the state	in G by processing all the transitions that come from it
			bool bLexUnitTransitions = false;		// whether actual lexical unit transitions (transitions that are actual words) are seen
			
			//printf("G state to be processed (%x) has: %d transitions\n",stateGFrom,stateGFrom->vTransition.size());
			//printf("processing left context: %s\n",m_phoneSet->getStrPhone(l));
			
			// (3.1) for each lexical-unit transition coming from the state, find the corresponding leaves in the lexical tree
			for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
			
				// epsilon transition: keep it in the final tree
				if ((*it)->iSymbol == EPSILON_TRANSITION) {
				
					// there can be at most one epsilon transition coming from any node
					if (transitionEpsilon != NULL) {
						assert(transitionEpsilon == *it);
					}
					transitionEpsilon = *it;	
					//printf("epsilon processed\n");
				} 
				// beginning/end of sentence: keep it in the final tree
				else if ((int)(*it)->iSymbol == (m_lexiconManager->m_lexUnitEndSentence->iLexUnit|LEX_UNIT_TRANSITION)) {
					//printf("end of sentence transition!\n");
					//printf("error: unexpected %s was processed\n",m_lexiconManager->getStrLexUnit((*it)->iSymbol));
					//exit(-1);
					transitionFinal = *it;
				}	
				// non-epsilon transition: put leaves from each pronunciation of the lexical unit in the queue
				else {	
				
					bLexUnitTransitions = true;
					LexUnitX *lexUnitX = m_lexiconManager->getLexUnit((*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
					assert(m_lexiconManager->isStandard(lexUnitX->iLexUnit) || m_lexiconManager->isFiller(lexUnitX->iLexUnit));
					for(VLexUnit::iterator jt = lexUnitX->vLexUnitPronunciations.begin() ; jt != lexUnitX->vLexUnitPronunciations.end() ; ++jt) {
					
						// get the lex unit leaves in the left-context lexical tree
						unsigned char iLeaves = lexiconTree[iLeftPhone]->iLeafRightContext[(*jt)->iLexUnitPron];
						LeafRightContext *leafRightContext = lexiconTree[iLeftPhone]->leafRightContext[(*jt)->iLexUnitPron];
						/*if (*jt == m_lexiconManager->getLexUnit("<UHM>",0)) {
							bool bStop = 0;
						}	*/
						
						for(unsigned int iLeaf = 0 ; iLeaf < iLeaves ; ++iLeaf) {
							LeafToProcess *leaf = new LeafToProcess;
							leaf->lexUnit = *jt;
							leaf->iHMMStates = (*jt)->vPhones.size()*NUMBER_HMM_STATES;
							leaf->stateLeaf = leafRightContext[iLeaf].stateLeaf;
							leaf->iRightContextGroup = leafRightContext[iLeaf].iRightContextGroup;
							leaf->stateGTo = (*it)->state;	
							leaf->fWeight = (*it)->fWeight;
							lLeafToProcess.push_back(leaf);	
						}
						/*++iLeavesQueued;
						if (iLeavesQueued % 1000 == 0) {
							printf("leaves queued: %d\n",iLeavesQueued);
						}*/
						
						// sanity check
						/*bool *bAux = new bool[m_phoneSet->size()];
						for(int i=0 ; i < m_phoneSet->size() ; ++i) {
							bAux[i] = false;
						}
						for(LLeafToProcess::iterator kt = lLeafToProcess.begin() ; kt != lLeafToProcess.end() ; ++kt) {
							unsigned char *iRightContextGroup = iRightContextGroups[(*kt)->iRightContextGroup];
							for(unsigned int i=0 ; iRightContextGroup[i] != UCHAR_MAX ; ++i) {
								bAux[iRightContextGroup[i]] = true;
							}
						}*/
						/*int iCounter = 0;
						for(int i=0 ; i < m_phoneSet->size() ; ++i) {
							if (bAux[i]) {
								++iCounter;
							}
						}
						if (iCounter != 43) {
							m_lexiconManager->printLexUnit(*jt);						
							printf("error: counter: %d expected: 43\n",iCounter);
							printf("left phone: %s\n",m_phoneSet->getStrPhone(iLeftPhone));
							for(LLeafToProcess::iterator kt = lLeafToProcess.begin() ; kt != lLeafToProcess.end() ; ++kt) {
								printf("leaf group: %d\n",(*kt)->iRightContextGroup);
							}
						} else {
							//m_lexiconManager->printLexUnit(*jt);
						}*/
						
						
						
						//m_lexiconManager->printLexUnit(*jt);
					}
					//delete *it;
				}	
			}
			
			// sort leaves in the queue by state number
			lLeafToProcess.sort(WFSABuilder::compareStateNumber);
			
			//printf("processing leaves...\n");
			
			//printf("leaves: %d (G transitions: %d, alternative pron: xxx)\n",lLeafToProcess.size(),stateGFrom->vTransition.size());
		
			// process leaves in the queue
			while(lLeafToProcess.empty() == false) {
			
				LeafToProcess *leaf = lLeafToProcess.front();
				lLeafToProcess.pop_front();
				
				//m_lexiconManager->printLexUnit(leaf->lexUnit);
				//printf("%d leaf: %x\n",leaf->lexUnit->iLexUnitPron,leaf->stateLeaf);
				
				int iHMMStates = leaf->iHMMStates;
				State *state = leaf->stateLeaf;		// this state goes to the lexical unit transition
				unsigned int iSymbol = ((unsigned int)leaf->lexUnit->iLexUnitPron)|LEX_UNIT_TRANSITION;
				
				// process states from the one that goes to the lexical unit transition to the one that goes to the first HMM
				for(int i=iHMMStates; i >= 0 ; --i) {
				
					//printf(" %x\n",state);
					//printStateTransitions(state);
				
					// empty position: insert the state into the buffer
					if (stateWaiting[i] == NULL) {
						stateWaiting[i] = state;
						transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
						assert((iSymbol >= LEX_UNIT_TRANSITION) || (m_hmmManager->getHMMStateDecoding(iSymbol)->getState() == i%NUMBER_HMM_STATES));
						checkTransitionSymbolCorrectness(iSymbol);
						++iTransitionsCreated;
						// if the state goes to a G node we already know the destination state in the final graph
						if (state == leaf->stateLeaf) {
							transitionInsertion[i].back()->fWeight = leaf->fWeight;
							transitionInsertion[i].back()->state = leaf->stateGTo;
							assert(transitionInsertion[i].back()->iSymbol & LEX_UNIT_TRANSITION);
							transitionInsertion[i].back()->iRightContextGroup = leaf->iRightContextGroup;
							//mLeafRightContextTemp.insert(MLeafTransitionRightContext::value_type(transitionInsertion[i].back(),leaf->iRightContextGroup));
						}
					}
					// different state: move the original state along with its transitions to the final tree
					else if (stateWaiting[i] != state) {
						// (1) the state in the buffer becomes final
						State *stateFinal = newState(iStateIdFinal++,NULL,NULL);
						++iStatesCreated;
						assert(i > 0);				// the root state stays in the buffer until the end
						assert(transitionInsertion[i].empty() == false);
						for(VTransition::iterator jt = transitionInsertion[i].begin() ; jt != transitionInsertion[i].end() ; ++jt) {
							stateFinal->vTransition.push_back(*jt);
							// if no destination state: look for the final state that does have destination state
							if ((*jt)->state == NULL) {
								assert(i == iHMMStates);
								bool bFound = false;
								for(int j = iTreeDepth-1 ; j > i ; --j) {
									if (stateWaiting[j] != NULL) {
										while(i != j) {
											State *stateFinal2 = newState(iStateIdFinal++,NULL,NULL);
											++iStatesCreated;
											for(VTransition::iterator kt = transitionInsertion[j].begin() ; kt != transitionInsertion[j].end() ; ++kt) {
												stateFinal2->vTransition.push_back(*kt);
												assert((*kt)->state != NULL);
											}
											stateWaiting[j] = NULL;
											transitionInsertion[j].clear();
											// apply weight pushing
											float fWeightRest2 = applyWeightPushing(stateFinal2);	
											// try to reuse a state
											MStateState::iterator ft = mStateState.find(stateFinal2);
											if (ft != mStateState.end()) {
												// remove the redundant transitions
												for(VTransition::iterator gt = stateFinal2->vTransition.begin() ; gt != stateFinal2->vTransition.end() ; ++gt) {
													delete *gt;
													--iTransitionsCreated;
												}
												delete stateFinal2;
												--iStateIdFinal;
												--iStatesCreated;	
												++iStatesReused;
												stateFinal2 = ft->second;	
												if (iStatesReused % 10000 == 0) {
													//printf("(1) reused: %u total: %u (transitions: %u)\n",iStatesReused,iStatesCreated,iTransitionsCreated);
												}
											} else {
												for(VTransition::iterator gt = stateFinal2->vTransition.begin() ; gt != stateFinal2->vTransition.end() ; ++gt) {
													if ((*gt)->iSymbol & LEX_UNIT_TRANSITION) {	
														vTransitionLexUnit.push_back(*gt);
													}
												}	
												mStateState.insert(MStateState::value_type(stateFinal2,stateFinal2));
											}
											// now the state is final
											Transition *transitionNULL2 = transitionInsertion[j-1].back();
											if (transitionInsertion[j-1].back()->state != NULL) {
												unsigned int iFound2 = 0;
												for(VTransition::iterator kt = transitionInsertion[j-1].begin() ; kt != transitionInsertion[j-1].end() ; ++kt) {
													if ((*kt)->state == NULL) {
														transitionNULL2 = *kt;
														++iFound2;
													}
												}
												assert(iFound2 == 1);
											}
											assert(transitionNULL2->state == NULL);
											transitionNULL2->state = stateFinal2;
											transitionNULL2->fWeight = fWeightRest2;
											--j;
										}
										bFound = true;
									}	
								}
								assert(bFound);
							}
							assert((*jt)->state != NULL);
						}
						transitionInsertion[i].clear();
						// apply weight pushing
						float fWeightRest = applyWeightPushing(stateFinal);	
						// try to reuse a state
						MStateState::iterator ft = mStateState.find(stateFinal);
						if (ft != mStateState.end()) {				
							// remove the redundant transitions
							for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
								delete *gt;
								--iTransitionsCreated;
							}
							delete stateFinal;
							iStateIdFinal--;
							--iStatesCreated;	
							++iStatesReused;
							stateFinal = ft->second;	
							if (iStatesReused % 10000 == 0) {
								//printf("(2) reused: %u total: %u (transitions: %u)\n",iStatesReused,iStatesCreated,iTransitionsCreated);
							}
						} else {
							for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
								if ((*gt)->iSymbol & LEX_UNIT_TRANSITION) {	
									vTransitionLexUnit.push_back(*gt);
								}
							}	
							mStateState.insert(MStateState::value_type(stateFinal,stateFinal));
						}
						//
						Transition *transitionNULL = transitionInsertion[i-1].back();
						if (transitionInsertion[i-1].back()->state != NULL) {
							unsigned int iFound = 0;
							for(VTransition::iterator jt = transitionInsertion[i-1].begin() ; jt != transitionInsertion[i-1].end() ; ++jt) {
								if ((*jt)->state == NULL) {
									transitionNULL = *jt;
									++iFound;
								}
							}
							assert(iFound == 1);
						}
						assert(transitionNULL->state == NULL);
						transitionNULL->state = stateFinal;
						transitionNULL->fWeight = fWeightRest;
						// (2) the new state is inserted into the buffer along with its transition
						stateWaiting[i] = state;
						transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
						assert((iSymbol >= LEX_UNIT_TRANSITION) || (m_hmmManager->getHMMStateDecoding(iSymbol)->getState() == i%NUMBER_HMM_STATES));
						checkTransitionSymbolCorrectness(iSymbol);
						++iTransitionsCreated;
						// if the state goes to a G node we already know the destination state in the final graph
						if (state == leaf->stateLeaf) {
							transitionInsertion[i].back()->fWeight = leaf->fWeight;
							transitionInsertion[i].back()->state = leaf->stateGTo;
							assert(transitionInsertion[i].back()->iSymbol & LEX_UNIT_TRANSITION);
							transitionInsertion[i].back()->iRightContextGroup = leaf->iRightContextGroup;
							//mLeafRightContextTemp.insert(MLeafTransitionRightContext::value_type(transitionInsertion[i].back(),leaf->iRightContextGroup));
						}
					} 
					// equal state: stop the right to left traversal
					else {
						// create a new transition
						transitionInsertion[i].push_back(newTransition(iSymbol,0.0,NULL));
						assert((iSymbol >= LEX_UNIT_TRANSITION) || (m_hmmManager->getHMMStateDecoding(iSymbol)->getState() == i%NUMBER_HMM_STATES));
						checkTransitionSymbolCorrectness(iSymbol);
						++iTransitionsCreated;
						// if the state goes to a G node we already know the destination state in the final graph
						if (state == leaf->stateLeaf) {
							transitionInsertion[i].back()->fWeight = leaf->fWeight;
							transitionInsertion[i].back()->state = leaf->stateGTo;
							assert(transitionInsertion[i].back()->iSymbol & LEX_UNIT_TRANSITION);
							transitionInsertion[i].back()->iRightContextGroup = leaf->iRightContextGroup;
							//mLeafRightContextTemp.insert(MLeafTransitionRightContext::value_type(transitionInsertion[i].back(),leaf->iRightContextGroup));
						}
						break;
					}
					
					if (i == 0) {
						assert(state->stateParent == NULL);
						break;	
					} 
					
					// move to the parent state
					iSymbol = state->transitionBW->iSymbol;
					state = state->stateParent;
				}
				
				++iLeavesProcessed;
				delete leaf;
			}
			
			// insert all the waiting states in the final graph (if any)
			if (bLexUnitTransitions) {
				//assert(stateWaiting[iTreeDepth-1] == NULL);
				assert(stateWaiting[0] != NULL);							// the root always has to be in the buffer	
				for(int i=iTreeDepth-1 ; i >= 0 ; --i) {
					if (stateWaiting[i] != NULL) {	
						while(i >= 0) {
							State *stateFinal = NULL;
							if (i > 0) {
								stateFinal = newState(iStateIdFinal++,NULL,NULL);
								++iStatesCreated;
							} else {
								// (i == 0) connect the left-context specific root to the final graph
								stateFinal = stateRootL;
							}
							for(VTransition::iterator kt = transitionInsertion[i].begin() ; kt != transitionInsertion[i].end() ; ++kt) {
								stateFinal->vTransition.push_back(*kt);
								assert((*kt)->state != NULL);
							}
							float fWeightRest = 0.0;
							if (i > 0) {
								// apply weight pushing
								fWeightRest = applyWeightPushing(stateFinal);
							}
							// try to reuse
							MStateState::iterator ft = mStateState.find(stateFinal);
							if (ft != mStateState.end()) {
								// remove the redundant transitions
								for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
									delete *gt;
									--iTransitionsCreated;
								}
								// i == 0
								if (stateFinal == stateRootL) {
									MGStateLRoot::iterator gt = mGStateLRoot.find(pair<State*,unsigned char>(stateGFrom,iLeftPhone));
									gt->second = ft->second;
								}	
								delete stateFinal;
								iStateIdFinal--;
								--iStatesCreated;	
								++iStatesReused;
								stateFinal = ft->second;	
								if (iStatesReused % 10000 == 0) {
									//printf("(3) reused: %u total: %u (transitions: %u)\n",iStatesReused,iStatesCreated,iTransitionsCreated);
								}
							} else {
								for(VTransition::iterator gt = stateFinal->vTransition.begin() ; gt != stateFinal->vTransition.end() ; ++gt) {
									if ((*gt)->iSymbol & LEX_UNIT_TRANSITION) {	
										vTransitionLexUnit.push_back(*gt);
									}
								}
								//mStateState.insert(MStateState::value_type(stateFinal,stateFinal));
							}	
							if (i > 0) {
								assert(transitionInsertion[i-1].back()->state == NULL);
								transitionInsertion[i-1].back()->state = stateFinal;
								transitionInsertion[i-1].back()->fWeight = fWeightRest;
							}
							//printStateTransitions(stateFinal);
							stateWaiting[i] = NULL;
							transitionInsertion[i].clear();
							--i;
						}
						break;
					}
				}
				// sanity check
				for(unsigned int i=0 ; i < iTreeDepth-1 ; ++i) {
					assert(stateWaiting[i] == NULL);
					assert(transitionInsertion[i].empty());
				}
			}	
			
			//printStateTransitions(stateRootL);
			//int iStop = 0;
		}
		
		// process the epsilon transition (if any)
		if (transitionEpsilon != NULL) {
				
			// it goes to a node that was already processed so the auxiliar left-context dependent roots must exist
			State *stateGTo = transitionEpsilon->state;
			// for each left context create an epsilon transition that transfers that context backwards
			for(unsigned char i =0 ; iLeftContext[stateGTo->iId][i] != UCHAR_MAX; ++i) {	
				unsigned char lLeft = iLeftContext[stateGTo->iId][i];
				MGStateLRoot::iterator jt = mGStateLRoot.find(pair<State*,unsigned char>(stateGTo,lLeft));
				assert(jt != mGStateLRoot.end());
				// sanity check, the L root must go to </s>
				bool bFound = false;
				for(VTransition::iterator nt = jt->second->vTransition.begin() ; nt != jt->second->vTransition.end() ; ++nt) {
					if ((int)((*nt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
						bFound = true;
						break;
					}
					if ((*nt)->iSymbol & EPSILON_TRANSITION) {		// for n-gram > 2
						bFound = true;
						break;
					}
				}
				assert(bFound);
				// epsilon transition: from an auxiliar L root to an auxiliar L root
				Transition *transitionAux = newTransition(EPSILON_TRANSITION,transitionEpsilon->fWeight,jt->second);
				++iTransitionsAuxiliarEpsilon;
				// try to reuse an auxiliar L root
				MGStateLRoot::iterator kt = mGStateLRoot.find(pair<State*,unsigned char>(stateGFrom,lLeft));
				if (kt != mGStateLRoot.end()) {
					kt->second->vTransition.push_back(transitionAux);
				} else {
					State *stateRootLAux = newState(UINT_MAX,NULL,NULL);			// this state wont be part of the final graph
					iStatesCreatedTemporal++;
					stateRootLAux->vTransition.push_back(transitionAux);
					mGStateLRoot.insert(MGStateLRoot::value_type(pair<State*,unsigned char>(stateGFrom,lLeft),stateRootLAux));
				}	
			}	
		}
		// process the transition to the final state </s> if any
		if (transitionFinal != NULL) {
		
			// for each posible left context	
			for(unsigned char l=0 ; iLeftContext[stateGFrom->iId][l] != UCHAR_MAX ; ++l) {
				MGStateLRoot::iterator kt = mGStateLRoot.find(pair<State*,unsigned char>(stateGFrom,iLeftContext[stateGFrom->iId][l]));
				assert(kt != mGStateLRoot.end());
				Transition *transitionFinalAux = newTransition(transitionFinal->iSymbol,transitionFinal->fWeight,stateFinalGDefinitive);
				kt->second->vTransition.push_back(transitionFinalAux);
			}
		}
		
		iTransitionsProcessed += stateGFrom->vTransition.size();
		iStatesProcessed++;
		
		if (iStatesProcessed % 100 == 0) {
			//printf("reused: %u (states: %u transitions: %u)\n",iStatesReused,iStatesCreated,iTransitionsCreated);
			printf("G states: %u (states: %u transitions: %u)\n",iStatesProcessed,iStatesCreated,iTransitionsCreated);
		}
		
		//printStateTransitions(stateGFrom);
		//printf("G node processed! (%u states processed, %u transitions processed)\n",iStatesProcessed,iTransitionsProcessed);
		//printf("(%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
		
		// remove the original transitions
		for(VTransition::iterator it = stateGFrom->vTransition.begin() ; it != stateGFrom->vTransition.end() ; ++it) {
			delete *it;
		}
		stateGFrom->vTransition.clear();
	}
	delete [] stateWaiting;
	delete [] transitionInsertion;
	
	// TODO TODO TODO TODO TODO TODO remove lexical units that come from states in G except <s>,</s> and epsilon ones
	
	// destroy all the lexical trees
	for(unsigned char iLeftContext = 0 ; iLeftContext < m_phoneSet->size() ; ++iLeftContext) {	
		if (lexiconTree[iLeftContext] != NULL) {
			destroyTree(lexiconTree[iLeftContext]->stateRoot);	
			unsigned int iLexiconSize = m_lexiconManager->getLexiconSize();
			for(unsigned int i=0 ; i < iLexiconSize ; ++i) {
				if (lexiconTree[iLeftContext]->leafRightContext[i] != NULL) {
					delete [] lexiconTree[iLeftContext]->leafRightContext[i];
				}
			}
			delete [] lexiconTree[iLeftContext]->iLeafRightContext;
			delete [] lexiconTree[iLeftContext]->leafRightContext;
			delete lexiconTree[iLeftContext];
		}
	}
	delete [] lexiconTree;
	
	for(unsigned int i=0 ; i < iStatesG ; ++i) {
		delete [] iLeftContext[i];
	}
	delete [] iLeftContext;
	
	printf("composition ends!!! (%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
	
	// (2) make connections
	map<State*,bool> mCheck;
	bool bPhones[100];
	bool bPhones2[100];
	for(int i=0 ; i<100 ; ++i) {
		bPhones[i] = false;
		bPhones2[i] = false;
	}
	
	// sanity check
	map<unsigned int,bool*> mLexUnitGroupMembers;
	for(VTransition::iterator it = vTransitionLexUnit.begin() ; it != vTransitionLexUnit.end() ; ++it) {
		
		//unsigned char *iRightContextGroup = iRightContextGroups[(*it)->iRightContextGroup];
 		unsigned char *iRightContextGroup = vRightContextGroups[(*it)->iRightContextGroup];;
		LexUnit *lexUnit = m_lexiconManager->getLexUnitPron((*it)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
		map<unsigned int,bool*>::iterator jt = mLexUnitGroupMembers.find(lexUnit->iLexUnitPron);
		bool *aux = NULL;
		if (jt == mLexUnitGroupMembers.end()) {
			aux = new bool[m_phoneSet->size()];
			for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
				aux[i] = false;
			}
			mLexUnitGroupMembers.insert(map<unsigned int,bool*>::value_type(lexUnit->iLexUnitPron,aux));
		} else {
			aux = jt->second;
		}
		for(unsigned int i=0 ; iRightContextGroup[i] != UCHAR_MAX ; ++i) {
			aux[iRightContextGroup[i]] = true;
		}	
	}
	unsigned int iCounterReference = UINT_MAX;
	cout << "lexical units with right context groups: " << mLexUnitGroupMembers.size() << endl;
	for(map<unsigned int,bool*>::iterator it = mLexUnitGroupMembers.begin() ; it != mLexUnitGroupMembers.end() ; ++it) {
		unsigned int iCounter = 0;
		for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
			if (it->second[i]) {
				++iCounter;
			}
		}		
		if (iCounterReference == UINT_MAX) {
			iCounterReference = iCounter;
		} 
		if (iCounter != iCounterReference) {
			printf("counter: %d expected: %u\n",iCounter,iCounterReference);
			m_lexiconManager->print(m_lexiconManager->getLexUnitPron(it->first));		
		}
		delete [] it->second;
		//assert(iCounter == iCounterReference);
	}
	
	// remove auxiliar L root nodes from the map of final nodes
	for(MGStateLRoot::iterator it = mGStateLRoot.begin() ; it != mGStateLRoot.end() ; ++it) {
		MStateState::iterator jt = mStateState.find(it->second);
		//assert(jt != mStateState.end());
		// note that there are duplicated values in the map and they can only be removed once
		if (jt != mStateState.end()) {
			assert(0);
			mStateState.erase(jt);
		}
	}	
	
	// get the first phonetic symbol of silence/filler lexical units
	unsigned char iFirstPhonesFillerSil = 0;
	bool *bFirstPhoneFillerSil = getFillerFirstPhones(&iFirstPhonesFillerSil);
	
	// for each leaf transition that needs to be connected	
	unsigned int iLeavesConnected = 0;
	unsigned int iLinks = 0;
	bool *bPhoneGroup = new bool[m_phoneSet->size()];
	VState vStateStack;
	vector<float> vWeight;
	for(VTransition::iterator it = vTransitionLexUnit.begin() ; it != vTransitionLexUnit.end() ; ++it) {
	
		Transition *transitionLexUnit = *it;
		LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(transitionLexUnit->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
		unsigned char iLeftContextAux = lexUnit->vPhones.back();
		State *stateGTo = transitionLexUnit->state;
		bPhones2[iLeftContextAux] = true;
		
		//printf("%-20s %u\n",m_lexiconManager->getStrLexUnitPron(lexUnit->iLexUnitPron),(*it)->iRightContextGroup);
			
		// get the right contexts of the leaf HMMs
		
		unsigned char *iRightContextGroup = vRightContextGroups[(*it)->iRightContextGroup];
		unsigned char iGroupElements = 0;
		for(unsigned char i=0; i < m_phoneSet->size() ; ++i) {
			bPhoneGroup[i] = false;
		}
		for(unsigned char i=0; iRightContextGroup[i] != UCHAR_MAX;++i) {
			bPhoneGroup[iRightContextGroup[i]] = true;
			++iGroupElements;
		}
		
		// if the leaf's right context is exclusively composed of filler phones (or silence), ignore the epsilon
		bool bIgnoreEpsilon = true;	
		for(unsigned char i=0; iRightContextGroup[i] != UCHAR_MAX;++i) {
			if (bFirstPhoneFillerSil[iRightContextGroup[i]] == false) {
				bIgnoreEpsilon = false;
				break;
			}
		}
		bIgnoreEpsilon = false;
		
		// get the auxiliar left-context specific root state from the state in G plus the
		MGStateLRoot::iterator jt = mGStateLRoot.find(pair<State*,unsigned char>(stateGTo,iLeftContextAux));	
		assert(jt != mGStateLRoot.end());
		State *stateRootL = jt->second;
		
		// put the auxiliar root state in the stack
		vStateStack.push_back(stateRootL);
		//printf("L inserted: %x\n",stateRootL);
		
		// add to the stack states coming from an epsilon transition
		if (bIgnoreEpsilon == false) {	
			//printf("looking for epsilon in: %x\n",stateRootL);
			bool bInserted;
			do {
				bInserted = false;
				State *stateHead = vStateStack.back();
				unsigned char iFound = 0;
				for(VTransition::iterator jt = stateHead->vTransition.begin() ; jt != stateHead->vTransition.end() ; ++jt) {
					if ((*jt)->iSymbol & EPSILON_TRANSITION) {
					
						// sanity check, the L root must go to </s> or epsilon
						bool bFound = false;
						for(VTransition::iterator nt = (*jt)->state->vTransition.begin() ; nt != (*jt)->state->vTransition.end() ; ++nt) {
							if ((int)((*nt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
								bFound = true;
								break;
							}
							if ((*nt)->iSymbol & EPSILON_TRANSITION) {
								bFound = true;
								break;
							}
						}
						assert(bFound);	
						//printf("epsilon and then </s>: %x %x\n",stateHead,(*jt)->state);
					
						vStateStack.push_back((*jt)->state);
						vWeight.push_back((*jt)->fWeight);
						bInserted = true;
						++iFound;
						//break;	// there cannot be more than one epsilon state leaving a state
					}
				}
				assert(iFound <= 1);
			} while(bInserted);
			
			// the first state in the stack (the one inserted the latest) has no epsilon transitions
		}
		
		State *stateTo = NULL;
		bool bStateToReused = false;
		float fWeightTo = 0.0;	
		int iTransitionsToFinal = 0;
		int iTransitionsToEpsilon = 0;
		
		// connect one by one states in the stack
		//printf("(%d) processing elements in the queue, %d elements\n",iProcessing++,vStateStack.size());
		while(vStateStack.empty() == false) {
		
			State *stateRootLAux = vStateStack.back();
			vStateStack.pop_back();
			//printf("rootL: %x\n",stateRootLAux);
		
			// replace the destination state in G by a new state in the final graph
			// the lexical unit transition goes to this state and this state goes to the first HMM-state leaves of the next subtree 
			State *stateLink = newState(iStateIdFinal++,NULL,NULL);
			++iStatesCreated;	
			
			// create a new transition for each of the right contexts	
			bool bConnected = false;
			bool bEpsilon = false;
			bool bFinal = false;
			for(VTransition::iterator kt = stateRootLAux->vTransition.begin() ; kt != stateRootLAux->vTransition.end() ; ++kt) {
				// epsilon transition
				if ((*kt)->iSymbol & EPSILON_TRANSITION) {
					bEpsilon = true;
					++iTransitionsToEpsilon;
				} 
				// </s> transition
				else if ((*kt)->iSymbol & LEX_UNIT_TRANSITION) {
					
					assert((int)((*kt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit);
					assert((*kt)->state == stateFinalGDefinitive);
					Transition *transitionFinal = newTransition((*kt)->iSymbol,(*kt)->fWeight,(*kt)->state);
					stateLink->vTransition.push_back(transitionFinal);
					++iTransitionsCreated;
					++iLinks;
					bFinal = true;	
					++iTransitionsToFinal;
				} 
				// HMM-state transition
				else {	
					assert((*kt)->iSymbol < LEX_UNIT_TRANSITION);
					// get the phone associated to the HMM-state
					unsigned char iPhone = m_hmmManager->getPhoneticSymbolFromHMMStateDecoding((*kt)->iSymbol);
					// check if the phone falls within the group
					if (bPhoneGroup[iPhone]) {
						// replicate the transition (the original cannot be reused since needs to remain unchanged for further processing)
						Transition *transitionLink = newTransition((*kt)->iSymbol,(*kt)->fWeight,(*kt)->state);
						stateLink->vTransition.push_back(transitionLink);
						++iTransitionsCreated;
						++iLinks;
						bConnected = true;
					}
				}
			}
			assert(bConnected || bEpsilon || bFinal);
			assert(bFinal || bEpsilon);
			
			// add the epsilon transition if needed
			if (stateTo != NULL) {
				assert(vWeight.empty() == false);
				float fWeightEpsilon = vWeight.back();
				vWeight.pop_back();
				Transition *transitionEpsilon = newTransition(EPSILON_TRANSITION,fWeightTo+fWeightEpsilon,stateTo);
				stateLink->vTransition.push_back(transitionEpsilon);
				++iTransitionsCreated;
				++iLinks;
			} else {
				assert(fWeightTo == 0.0);
			}
			
			// try to reuse a state
			
			// apply weight pushing
			float fWeightRest = applyWeightPushing(stateLink);	
			// try to reuse a state
			MStateState::iterator ft = mStateState.find(stateLink);
			if (ft != mStateState.end()) {
				// remove the redundant transitions
				for(VTransition::iterator gt = stateLink->vTransition.begin() ; gt != stateLink->vTransition.end() ; ++gt) {
					delete *gt;
					--iTransitionsCreated;
					--iLinks;	
				}
				delete stateLink;
				--iStateIdFinal;
				--iStatesCreated;	
				++iStatesReused;
				stateLink = ft->second;	
				bStateToReused = true;
			} else {
				/*for(VTransition::iterator gt = stateLink->vTransition.begin() ; gt != stateLink->vTransition.end() ; ++gt) {
					assert(((*gt)->iSymbol & LEX_UNIT_TRANSITION) == 0);
				}	*/
				assert(stateLink->iId < iStateIdFinal);
				mStateState.insert(MStateState::value_type(stateLink,stateLink));
				bStateToReused = false;
			}
			assert(stateLink->iId < iStateIdFinal);
				
			stateTo = stateLink;	
			fWeightTo = fWeightRest;
			iLeavesConnected++;	
		}
		assert(vWeight.empty());
		if (iTransitionsToFinal < 1) {
			printf("bad: %d %d\n",iTransitionsToFinal,iTransitionsToEpsilon);
		} else {
			//printf("good\n");
		}
		assert(iTransitionsToFinal >= 1);
		
		transitionLexUnit->state = stateTo;
		transitionLexUnit->fWeight += fWeightTo;			// this is what makes global minimization needed
		
		// connect the state that replaces the G state to the final state if not already connected
		/*if (stateFinalG != NULL) {
			if (bStateToReused == false) {	
				// connect to the final state in G if necessary
				bool bGoesToFinal = false;
				bool bConnected = false;
				for(VTransition::iterator jt = stateGTo->vTransition.begin() ; jt != stateGTo->vTransition.end() ; ++jt) {
					if ((*jt)->state == stateFinalG) {	
						// sanity check
						for(VTransition::iterator kt = stateTo->vTransition.begin() ; kt != stateTo->vTransition.end() ; ++kt) {
							assert((*kt)->state != stateFinalGDefinitive);
						}
						// connect the state to the final state	
						assert((*jt)->iSymbol & LEX_UNIT_TRANSITION);
						int iLexUnitPron = (*jt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
						assert(iLexUnitPron == m_lexiconManager->m_lexUnitEndSentence->iLexUnit);
						Transition *transitionFinal = newTransition((*jt)->iSymbol,(*jt)->fWeight,stateFinalGDefinitive);
						stateTo->vTransition.push_back(transitionFinal);
						++iTransitionsCreated;
						++iLinks;
						bConnected = true;
						break;
					}
				}
				assert(bConnected == true);
			}
		}*/
	}
	delete [] bPhoneGroup;
	delete [] bFirstPhoneFillerSil;
	
	assert(stateFinalGDefinitive->vTransition.empty());
	
	// destroy the right context groups
	for(vector<unsigned char*>::iterator it = vRightContextGroups.begin() ; it != vRightContextGroups.end() ; ++it) {
		delete [] *it;	
	}	
	
	// (2.1) connect the root node of the final graph to all the starting transitions
	MStateBool mStateBoolSeen;		// multiple entries in mGStateLRoot may have the same value
	assert(vStateStack.empty());
	for(unsigned char i=0 ; i<m_phoneSet->size() ; ++i) {
		MGStateLRoot::iterator it = mGStateLRoot.find(pair<State*,unsigned char>(stateRootG,i));
		if (it != mGStateLRoot.end()) {
			assert(mStateBoolSeen.find(it->second) == mStateBoolSeen.end());
			if (mStateBoolSeen.find(it->second) == mStateBoolSeen.end()) {
				
				State *stateTo = NULL;
				//float fWeightTo = 0.0;	
				State *stateRootL = it->second;	
				// put the auxiliar root state in the stack
				vStateStack.push_back(stateRootL);
				// add to the stack states coming from an epsilon transition
				bool bInserted = false;
				do {
					bInserted = false;
					State *stateHead = vStateStack.back();
					unsigned char iFound = 0;
					for(VTransition::iterator jt = stateHead->vTransition.begin() ; jt != stateHead->vTransition.end() ; ++jt) {
						if ((*jt)->iSymbol & EPSILON_TRANSITION) {
							vStateStack.push_back((*jt)->state);
							vWeight.push_back((*jt)->fWeight);
							bInserted = true;
							++iFound;
							//break;	// there cannot be more than one epsilon state leaving a state
						}
					}
					assert(iFound <= 1);
				} while(bInserted);
				
				// the first state in the stack (the one inserted the latter) has no epsilon transitions
				
				// process states in the stack
				assert(vStateStack.empty() == false);
				while(vStateStack.empty() == false) {
				
					State *stateRootLAux = vStateStack.back();
					vStateStack.pop_back();	
					
					State *stateLink = NULL;
					// create a state that will replace the "auxiliar root" state
					if (vStateStack.empty() == false) {
						stateLink = newState(iStateIdFinal++,NULL,NULL);
						++iStatesCreated;	
					}
					// use the root state (there is no need to create an auxiliar state) 
					else {
						stateLink = stateRoot;
					}	
					
					// create a new transition for each of the right contexts	
					bool bConnected = false;
					bool bEpsilon = false;
					bool bFinal = false;
					for(VTransition::iterator kt = stateRootLAux->vTransition.begin() ; kt != stateRootLAux->vTransition.end() ; ++kt) {
						// epsilon transition
						if ((*kt)->iSymbol & EPSILON_TRANSITION) {
							bEpsilon = true;
							continue;
						}
						// </s> transition
						else if ((*kt)->iSymbol & LEX_UNIT_TRANSITION) {
							
							assert((int)((*kt)->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit);
							Transition *transitionFinal = newTransition((*kt)->iSymbol,(*kt)->fWeight,(*kt)->state);
							stateLink->vTransition.push_back(transitionFinal);
							++iTransitionsCreated;
							++iLinks;
							bFinal = true;	
						} 
						// HMM-state transition
						else {
							assert((*kt)->iSymbol < LEX_UNIT_TRANSITION);
							// replicate the transition (the original cannot be reused since needs to remain unchanged)
							Transition *transitionLink = newTransition((*kt)->iSymbol,(*kt)->fWeight,(*kt)->state);
							stateLink->vTransition.push_back(transitionLink);
							++iTransitionsCreated;
							++iLinks;
							bConnected = true;	
						}
					}
					assert(bConnected || bEpsilon || bFinal);
					
					// add the epsilon transition if needed
					if (stateTo != NULL) {
						assert(vWeight.empty() == false);
						float fWeightEpsilon = vWeight.back();
						vWeight.pop_back();
						//Transition *transitionEpsilon = newTransition(EPSILON_TRANSITION,fWeightTo+fWeightEpsilon,stateTo);
						Transition *transitionEpsilon = newTransition(EPSILON_TRANSITION,fWeightEpsilon,stateTo);
						stateLink->vTransition.push_back(transitionEpsilon);
						++iTransitionsCreated;
						++iLinks;
					} else {
						//assert(fWeightTo == 0.0);
					}
					
					// try to reuse a state
					
					// apply weight pushing
					//float fWeightRest = 0.0;	
					
					if (vStateStack.empty() == false) {
					
						//fWeightRest = applyWeightPushing(stateLink);
					
						// try to reuse a state
						
						stateTo = stateLink;	
						//fWeightTo = fWeightRest;	
					}
					assert(stateLink->iId < iStateIdFinal);
				}
				assert(vWeight.empty());
				
				/*for(VTransition::iterator jt = stateRootL->vTransition.begin() ; jt != stateRootL->vTransition.end() ; ++jt) {
					if ((*jt)->iSymbol < LEX_UNIT_TRANSITION) {
						Transition *transitionStart = newTransition((*jt)->iSymbol,(*jt)->fWeight,(*jt)->state);
						stateRoot->vTransition.push_back(transitionStart);
						++iTransitionsCreated;
					} else {
						assert((*jt)->iSymbol == EPSILON_TRANSITION);
						// this transition goes to an auxiliar L root
						if (mStateBoolSeen.find((*jt)->state) == mStateBoolSeen.end()) {
							State *stateRootL2 = (*jt)->state;
							for(VTransition::iterator lt = stateRootL2->vTransition.begin() ; lt != stateRootL2->vTransition.end() ; ++lt) {
								assert((*lt)->iSymbol < LEX_UNIT_TRANSITION);
								Transition *transitionStart2 = newTransition((*lt)->iSymbol,(*lt)->fWeight+(*jt)->fWeight,(*lt)->state);
								stateRoot->vTransition.push_back(transitionStart2);
								++iTransitionsCreated;	
							}
							mStateBoolSeen.insert(MStateBool::value_type(stateRootL2,true));
						}
					}
				}*/
				
				mStateBoolSeen.insert(MStateBool::value_type(it->second,true));
			} 
		}
	}
	
	//delete [] statesG;
	
	printf("leaf transitions connected: %u, links created: %u\n",iLeavesConnected,iLinks);	
	printf("connection ends!!! (%u states created, %u transitions created)\n",iStatesCreated,iTransitionsCreated);
	
	unsigned int iTrue = 0;
	unsigned int iTrue2 = 0;
	for(int i=0 ; i< 100 ; ++i) {
		if (bPhones[i]) {
			++iTrue;
		}
		if (bPhones2[i]) {
			++iTrue2;
		}
	}
	
	// there can be duplicated values and they must be freed just once
	/*MStateBool mStateRemoved;
	for(MGStateLRoot::iterator it = mGStateLRoot.begin() ; it != mGStateLRoot.end() ; ++it) {
		if (mStateRemoved.find(it->second) == mStateRemoved.end()) {
			for(VTransition::iterator jt = it->second->vTransition.begin() ; jt != it->second->vTransition.end() ; ++jt) {	
				--iTransitionsCreated;
				delete *jt;
			}
			it->second->vTransition.clear();
			delete it->second;
			//--iStatesCreated;   not needed since it is a temporal state
			mStateRemoved.insert(MStateBool::value_type(it->second,true));
		}
	}
	mGStateLRoot.clear();*/
			
	//exit(-1);
	
	// sanity checks: 
	// (1) make sure there are not unconnected states
	// (2) make sure transition symbols are correct
	int iHMMStatesTotal = -1;
	m_hmmManager->getHMMStates(&iHMMStatesTotal);
	unsigned int iTransitionsSeen = 0;
	unsigned int iStatesSeen = 0;
	bool *bStateSeen = new bool[iStateIdFinal];
	for(unsigned int i=0 ; i < iStateIdFinal ; ++i) {
		bStateSeen[i] = false;
	}
	
	LState lState;
	lState.push_back(stateRoot);
	bStateSeen[stateRoot->iId] = true;
	++iStatesSeen;
	while(lState.empty() == false) {
		
		// get the next state to process
		State *state = lState.front();
		lState.pop_front();
		
		// (3.1) for each lexical-unit transition coming from the state, find the corresponding leaf in the lexical tree
		//assert(stateGFrom->vTransition.empty() == false);
		for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		
			checkTransitionSymbolCorrectness((*it)->iSymbol);
		
			++iTransitionsSeen;
			if (bStateSeen[(*it)->state->iId] == false) {
				lState.push_back((*it)->state);
				bStateSeen[(*it)->state->iId] = true;
				++iStatesSeen;
			}
		}	
	}
	printf("states created: %u, states seen: %u\n",iStatesCreated,iStatesSeen);
	for(unsigned int i=0 ; i < iStateIdFinal ; ++i) {
		if (bStateSeen[i] == false) {
			printf("state not seen: %d\n",i);
		}
		//assert(bStateSeen[i]);
	}
	delete [] bStateSeen;
	assert(iStatesCreated == iStatesSeen);
	
	iTransitionsCreated = iTransitionsSeen;
	printf("transitions created: %u, transitions seen: %u\n",iTransitionsCreated,iTransitionsSeen);
	assert(iTransitionsSeen == iTransitionsCreated);	
	
	unsigned int iStatesTotal = iStatesCreated;
	unsigned int iTransitionsTotal = iTransitionsCreated;
	float fNetworkSize = iStatesTotal*sizeof(StateX)+iTransitionsTotal*sizeof(TransitionX);
	printf("total states: %u, total transitions: %u\n",iStatesTotal,iTransitionsTotal);
	printf("Expected network size: %f MB\n",fNetworkSize/(1024.0*1024.0));
	
	//exit(-1);
	
	// (4) global minimization (includes global weight pushing)
	/*float fWeightInitial = 0;
	float fWeightFinal = 0;
	globalWeightPushing(stateRoot,stateFinalGDefinitive,iStatesCreated,&fWeightInitial,&fWeightFinal);*/
	
	// (5) equalize the acceptor input
	unsigned int iStatesTotalFinal = 0;
	unsigned int iTransitionsTotalFinal = 0;
	equalizeInput(stateRoot,iStatesTotal,iTransitionsTotal,&iStatesTotalFinal,&iTransitionsTotalFinal);
	
	// (6) transform the resulting graph to an acceptor that can be used for decoding
	WFSAcceptor *wfsAcceptor = optimize(stateRoot,iStatesTotalFinal,iTransitionsTotalFinal,stateFinalGDefinitive);

	return wfsAcceptor;
}



// perform global weight pushing, which is the first step of weighted minimization, the second step is the actual minimization
// note: it only works for the tropical semiring, for the log-semirint there are alternative methods
// note: there is always only one final state which is the destination state for end-of-sentence transitions
void WFSABuilder::globalWeightPushing(State *stateInitial, State *stateFinal, int iStates, 
	float *fWeightInitial, float *fWeightFinal) {

	// Dijkstra-based implementation (it is necessary to traverse the graph backwards so the list of predecessors
	// of each state needs to be available)
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// (1) for each state create an array containing its predecessors
	
	StatePredecessors *statePredecessors = new StatePredecessors[iStates];	
	bool *bState = new bool[iStates];	
	for(int i = 0; i < iStates; ++i) {
		statePredecessors[i].statePredecessor = NULL;
		statePredecessors[i].iPredecessors = 0;
		statePredecessors[i].transitionPredecessor = NULL;
		bState[i] = false;
	}
	bState[stateInitial->iId] = true;
	
	VState vState;
	vState.push_back(stateInitial);
	while(vState.empty() == false) {
	
		State *state = vState.back();
		vState.pop_back();
		
		for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		
			if (statePredecessors[(*it)->state->iId].statePredecessor == NULL) {
				statePredecessors[(*it)->state->iId].statePredecessor = new State*[1];
				statePredecessors[(*it)->state->iId].statePredecessor[0] = state;
				statePredecessors[(*it)->state->iId].transitionPredecessor = new Transition*[1];
				statePredecessors[(*it)->state->iId].transitionPredecessor[0] = *it;
				statePredecessors[(*it)->state->iId].iPredecessors = 1;
			} else {
				int iPredecessorsOriginal = statePredecessors[(*it)->state->iId].iPredecessors;	
				State **stateAux = new State*[iPredecessorsOriginal+1];
				Transition **transitionAux = new Transition*[iPredecessorsOriginal+1];	
				for(int i=0 ; i < iPredecessorsOriginal ; ++i) {
					stateAux[i] = statePredecessors[(*it)->state->iId].statePredecessor[i];
					transitionAux[i] = statePredecessors[(*it)->state->iId].transitionPredecessor[i];
				}
				stateAux[iPredecessorsOriginal] = state;
				transitionAux[iPredecessorsOriginal] = *it;
				delete [] statePredecessors[(*it)->state->iId].statePredecessor;
				statePredecessors[(*it)->state->iId].statePredecessor = stateAux;	
				delete [] statePredecessors[(*it)->state->iId].transitionPredecessor;
				statePredecessors[(*it)->state->iId].transitionPredecessor = transitionAux;
				statePredecessors[(*it)->state->iId].iPredecessors++;
			}	
			
			if (bState[(*it)->state->iId] == false) {
				vState.push_back((*it)->state);
				bState[(*it)->state->iId] = true;
			}
		}
	}
	printf("array of predecessors created\n");
	
	// (2) compute the minimum distance from each state to the final state using Dijkstra
	
	unsigned char *iStateStatus = new unsigned char [iStates];		// keeps the state status
	float *fStateWeightFinal = new float[iStates];						// keeps the weight from the state to the final state
	
	// initialization
	for(int i=0; i < iStates ; ++i) {
		iStateStatus[i] = STATE_STATUS_UNSEEN;
		fStateWeightFinal[i] = -FLT_MAX;
	}
	
	assert(vState.empty());
	vState.push_back(stateFinal);
	iStateStatus[stateFinal->iId] = STATE_STATUS_QUEUED;
	fStateWeightFinal[stateFinal->iId] = 0.0;
	int iIterations = 0;
	
	while(vState.empty() == false) {
	
		++iIterations;
		if ((iIterations % 1000000) == 0) {
			printf("iterations: %d\n",iIterations);
		}
	
		State *state = vState.back();
		vState.pop_back();	
		assert(fStateWeightFinal[state->iId] != -FLT_MAX);
		iStateStatus[state->iId] = STATE_STATUS_PROCESSED;
		
		for(int i=0 ; i < statePredecessors[state->iId].iPredecessors ; ++i) {
		
			State *statePredecessor = statePredecessors[state->iId].statePredecessor[i];
			float fWeight = statePredecessors[state->iId].transitionPredecessor[i]->fWeight;
			
			// is the new weight smaller?
			if ((fWeight + fStateWeightFinal[state->iId]) > fStateWeightFinal[statePredecessor->iId]) {
				fStateWeightFinal[statePredecessor->iId] = fWeight + fStateWeightFinal[state->iId];
				if (iStateStatus[statePredecessor->iId] != STATE_STATUS_QUEUED) {
					vState.push_back(statePredecessor);
					iStateStatus[statePredecessor->iId] = STATE_STATUS_QUEUED;
				}	
			} else if (iStateStatus[statePredecessor->iId] == STATE_STATUS_UNSEEN) {
				vState.push_back(statePredecessor);
				iStateStatus[statePredecessor->iId] = STATE_STATUS_QUEUED;
			}
		}
	}
	
	// sanity check
	for(int i=0 ; i < iStates ; ++i) {
		assert((fStateWeightFinal[i] != -FLT_MAX) || (i == (int)stateFinal->iId));
	}
	printf("Dijkstra completed\n");
	
	// (3) update the weights
	
	for(int i=0 ; i < iStates; ++i) {
		bState[i] = false;
	}	
	int iTransitions = 0;
	int iTransitionsUpdated = 0;
	
	assert(vState.empty());
	vState.push_back(stateInitial);
	while(vState.empty() == false) {
	
		State *state = vState.back();
		vState.pop_back();
		
		for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
		
			//printTransition(*it);
			float fWeightOriginal = (*it)->fWeight;
			(*it)->fWeight += fStateWeightFinal[(*it)->state->iId] - fStateWeightFinal[state->iId];
			if ((*it)->fWeight != fWeightOriginal) {
				//printf("weight update: %f -> %f\n",fWeightOriginal,(*it)->fWeight);
				++iTransitionsUpdated;
			}	
			++iTransitions;
			if (bState[(*it)->state->iId] == false) {
				vState.push_back((*it)->state);
				bState[(*it)->state->iId] = true;
			}
		}
	}	
	*fWeightInitial = fStateWeightFinal[stateInitial->iId];
	*fWeightFinal = 0;	
	
	float fPercent = 100.0*(((float)iTransitionsUpdated)/((float)iTransitions));	
	
	// clean-up
	delete [] iStateStatus;
	delete [] fStateWeightFinal;
	for(int i=0 ; i < iStates ; ++i) {
		if (statePredecessors[i].statePredecessor != NULL) {
			delete [] statePredecessors[i].statePredecessor;
			delete [] statePredecessors[i].transitionPredecessor;
		}
	}
	delete [] statePredecessors;
	delete [] bState;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();

	printf("-- global weight pushing ----------------------------\n");
	printf(" transitions: %d\n",iTransitions);
	printf(" transitions with updated weight: %d (%.2f%%)\n",iTransitionsUpdated,fPercent);
	printf(" initial state weight: %.2f\n",*fWeightInitial);
	printf(" final state weight:   %.2f\n",*fWeightFinal);
	printf(" time: %.2f seconds\n",(dTimeEnd-dTimeBegin)/1000.0);
	printf("-----------------------------------------------------\n");	
}

};	// end-of-namespace
