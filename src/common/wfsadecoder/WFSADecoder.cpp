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
#include "HMMManager.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "FeatureExtractor.h"
#include "TimeUtils.h"
#include "WFSADecoder.h"

namespace Bavieca {

// constructor
WFSADecoder::WFSADecoder(PhoneSet *phoneSet, HMMManager *hmmManager,
	LexiconManager *lexiconManager, WFSAcceptor *wfsAcceptor, int iPruningMaxActiveStates, float fPruningLikelihood, bool bLatticeGeneration, int iMaxWordSequencesState) {

	m_phoneSet = phoneSet;	
	m_hmmManager = hmmManager;
	m_lexiconManager = lexiconManager;
	m_wfsAcceptor = wfsAcceptor;
	
	// pruning
	m_iPruningMaxActiveStates = iPruningMaxActiveStates;
	m_fPruningLikelihood = fPruningLikelihood;
	
	// lattice generation
	m_bLatticeGeneration = bLatticeGeneration;
	if (bLatticeGeneration) {
		m_iMaxWordSequencesState = iMaxWordSequencesState;
		assert(m_iMaxWordSequencesState >= 2);
	}
}

// destructor
WFSADecoder::~WFSADecoder() {

	delete m_activeStateTable;
}

// initialization
void WFSADecoder::initialize() {

	// get the array of HMM-states
	m_hmmStatesDecoding = m_hmmManager->getHMMStatesDecoding((int*)&m_iHMMStatesDecoding);
	assert(m_hmmStatesDecoding != NULL);
	//m_hmmManager->initializeDecoding();
	
	// create and initialize the active state table
	m_activeStateTable = new ActiveStateTable(m_fPruningLikelihood,m_iPruningMaxActiveStates,100000,m_phoneSet,m_lexiconManager,m_hmmStatesDecoding,m_bLatticeGeneration,m_iMaxWordSequencesState);
	m_activeStateTable->initialize();
}


// initialize the Viterbi search
void WFSADecoder::initializeViterbi() {

}


// Viterbi search
void WFSADecoder::viterbi(Matrix<float> &mFeatures) {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	float fScore = 0.0;
	HMMStateDecoding *hmmStateDecoding;
	
	m_iFeatureVectors = mFeatures.getRows();
	
	m_hmmManager->resetHMMEmissionProbabilityComputation();
	
	m_activeStateTable->beginUtterance();	
	
	unsigned int iLexUnitEndSentence = m_lexiconManager->m_lexUnitEndSentence->iLexUnit;
		
	// create the initial set of active states (states coming from the initial state and non-epsilon arcs)	
	StateX *stateInitial = m_wfsAcceptor->getInitialState();
	m_activeStateTable->activateStateInitial(stateInitial);
	
	TransitionX *transition = NULL;
	TransitionX *transitionEnd = NULL;
	TransitionX *transitionAux = NULL;
	ActiveState *activeStatesCurrent = NULL;
	unsigned int iActiveStatesCurrent = 0;
	
	//m_activeStateTable->printSubgraph(stateInitial);
	//exit(-1);
	
	// process the feature vectors
	for(unsigned int t=0 ; t < mFeatures.getRows() ; ++t) {
		
		m_fScoreBest = -FLT_MAX;
		VectorStatic<float> vFeatureVector = mFeatures.getRow(t);
		
		activeStatesCurrent = m_activeStateTable->getActiveStatesCurrent(&iActiveStatesCurrent);
		
		m_activeStateTable->m_iStatesExpanded = 0;
		m_activeStateTable->m_iStatesActivated = 0;
		m_activeStateTable->m_iStatesPruned = 0;
		
		// (1) process all the active states (these are non-epsilon states)
		for(unsigned int i = 0 ; i < iActiveStatesCurrent ; ++i) {
		
			// skip pruned states
			if (activeStatesCurrent[i].state == NULL) {
				continue;
			}
			
			bool bExpanded = false;
			
			++m_activeStateTable->m_iStatesExpanded;
			
			ActiveState &activeState = activeStatesCurrent[i];	
			
			// (1.1) self-loop (this is a simulated transition)
			
			// compute emission probability
			fScore = activeState.hmmStateDecoding->computeEmissionProbability(vFeatureVector.getData(),t);	
			
			// preventive pruning goes here
			if (activeState.fScore+fScore < (m_fScoreBest-m_fPruningLikelihood)) {
				continue;
			}
			
			// lattice generation
			int iWGToken = -1;
			if (m_bLatticeGeneration) {
				assert(activeState.iWGToken != -1);
				iWGToken = m_activeStateTable->newWGToken(activeState.iWGToken);
				WGToken *wgToken = iWGToken+m_activeStateTable->m_wgTokens;
				// update the scores
				for(int i=0 ; ((i < m_activeStateTable->m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1)) ; ++i) {
					wgToken[i].fScore += fScore;
				}
			}
					
			// activate the state
			m_activeStateTable->activateState(activeState.state,activeState.fScore+fScore,&m_fScoreBest,
				activeState.hmmStateDecoding,activeState.iHistoryItem,iWGToken,0.0);	
				
			// (1.2) standard transitions
			
			// get the first and last transitions from the state
			transition = *activeState.state;
			transitionEnd = *(activeState.state+1);
			assert(transition <= transitionEnd);
			
			unsigned int iTransitions = 0;
			TransitionX *transitions = m_wfsAcceptor->getTransitions(&iTransitions);
			assert(transition < (transitions+iTransitions));
			assert(transitionEnd <= (transitions+iTransitions));
			
			while(transition != transitionEnd) {
			
				assert(transition != NULL);
			
				// epsilon-transition
				if (transition->iSymbol & EPSILON_TRANSITION) {
				
					// lattice generation
					int iWGToken = -1;
					if (m_bLatticeGeneration) {
						assert(activeState.iWGToken != -1);
						iWGToken = m_activeStateTable->newWGToken(activeState.iWGToken);
						WGToken *wgToken = iWGToken+m_activeStateTable->m_wgTokens;
						// update the scores
						for(int i=0 ; ((i < m_activeStateTable->m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1)) ; ++i) {
							wgToken[i].fScore += transition->fWeight;
						}
					}
				
					// activate the state
					m_activeStateTable->activateStateEpsilon(transition->state,activeState.fScore+transition->fWeight,
						activeState.iHistoryItem,iWGToken,0.0);
				}
				// fake-transition
				else if (transition->iSymbol & FAKE_TRANSITION) {
					//bool bStop = true;
				}	
				// lexical-unit transition
				else if (transition->iSymbol & LEX_UNIT_TRANSITION) {
				
					// create a new history item
					int iHistoryItem = m_activeStateTable->newHistoryItem();
					HistoryItem *historyItem = m_activeStateTable->m_historyItems+iHistoryItem;
					historyItem->iPrev = activeState.iHistoryItem;
					historyItem->iLexUnitPron = transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
					historyItem->iEndFrame = t-1;
					historyItem->fScore = activeState.fScore+transition->fWeight;
					historyItem->iActive = -1;
					
					// lattice generation
					int iWGToken = -1;
					int iWordSequence = -1;
					if (m_bLatticeGeneration) {
						// keep the N word sequences arriving to the state
						assert(activeState.iWGToken != -1);
						historyItem->iWGToken = activeState.iWGToken;
						// checks
						for(int i=0 ; i < m_iMaxWordSequencesState ; ++i) {
							if ((historyItem->iWGToken+m_activeStateTable->m_wgTokens)[i].iWordSequence == -1) {
								break;
							}
							assert(historyItem->iEndFrame > (m_activeStateTable->m_historyItems+ (historyItem->iWGToken+m_activeStateTable->m_wgTokens)[i].iHistoryItem)->iEndFrame);
						}
						// generate a new hash value for the new word sequence
						iWordSequence = m_activeStateTable->hashWordSequence(historyItem);	
					}
					
					if ((transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == iLexUnitEndSentence) {

					}
					
					for(transitionAux = *(transition->state) ; transitionAux != *(transition->state+1) ; ++transitionAux) {
						
						// epsilon-transition
						if (transitionAux->iSymbol & EPSILON_TRANSITION) {
						
							// preventive pruning
							if (activeState.fScore+transition->fWeight+transitionAux->fWeight < (m_fScoreBest-m_fPruningLikelihood)) {
								continue;
							}
							
							// lattice generation
							if (m_bLatticeGeneration) {
								bExpanded = true;
								iWGToken = m_activeStateTable->newWGToken(iWordSequence,
									activeState.fScore+transition->fWeight+transitionAux->fWeight,iHistoryItem);
							}
						
							// activate the state
							m_activeStateTable->activateStateEpsilon(transitionAux->state,
								activeState.fScore+transition->fWeight+transitionAux->fWeight,iHistoryItem,iWGToken,0.0);
						} 
						// fake transition
						else if (transitionAux->iSymbol & FAKE_TRANSITION) {
							//bool bStop = true;
						} 
						// lex-unit transition (end-of-sentence)
						else if (transitionAux->iSymbol & LEX_UNIT_TRANSITION) {
							if ((transitionAux->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == iLexUnitEndSentence) {
								//bool bStop = true;
							}	
						}
						// leaf-transition
						else {
						
							// compute emission probability
							hmmStateDecoding = &m_hmmStatesDecoding[transitionAux->iSymbol];
							fScore = hmmStateDecoding->computeEmissionProbability(vFeatureVector.getData(),t);	
							
							// preventive pruning
							if (activeState.fScore+transition->fWeight+transitionAux->fWeight+fScore < (m_fScoreBest-m_fPruningLikelihood)) {
								continue;
							}
						
							// lattice generation
							if (m_bLatticeGeneration) {
								bExpanded = true;
								iWGToken = m_activeStateTable->newWGToken(iWordSequence,
									activeState.fScore+transition->fWeight+transitionAux->fWeight+fScore,iHistoryItem);
							}
							
							// activate the state	
							m_activeStateTable->activateState(transitionAux->state,
								activeState.fScore+transition->fWeight+transitionAux->fWeight+fScore,
								&m_fScoreBest,hmmStateDecoding,iHistoryItem,iWGToken,0.0);
						}
					}
				}
				// leaf-transition
				else {
				
					// compute emission probability
					hmmStateDecoding = &m_hmmStatesDecoding[transition->iSymbol];	
					fScore = hmmStateDecoding->computeEmissionProbability(vFeatureVector.getData(),t);
					
					// preventive pruning goes here
					if (activeState.fScore+transition->fWeight+fScore < (m_fScoreBest-m_fPruningLikelihood)) {
						transition++;
						continue;
					}
					
					// lattice generation
					int iWGToken = -1;
					if (m_bLatticeGeneration) {
						assert(activeState.iWGToken != -1);
						iWGToken = m_activeStateTable->newWGToken(activeState.iWGToken);
						WGToken *wgToken = iWGToken+m_activeStateTable->m_wgTokens;
						// update the scores
						for(int i=0 ; ((i < m_activeStateTable->m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1)) ; ++i) {
							wgToken[i].fScore += transition->fWeight+fScore;
						}
					}
					
					// activate the state	
					m_activeStateTable->activateState(transition->state,activeState.fScore+transition->fWeight+fScore,
						&m_fScoreBest,hmmStateDecoding,activeState.iHistoryItem,iWGToken,0.0);
				}
				
				transition++;
			}
		}
		//printf("# active states before processing epsilon-transitions: %d\n",);
		
		// sanity check 
		if (m_bLatticeGeneration) {
			// active nodes of next time frame
			ActiveState *activeState = m_activeStateTable->m_activeStatesNext;
			while(activeState != m_activeStateTable->m_activeStateAvailable) {
			
				// make sure the WGToken structure keeps the historyItem
				if (activeState->iWGToken != -1) {
					assert(activeState->iHistoryItem != -1);
					WGToken *wgToken = activeState->iWGToken+m_activeStateTable->m_wgTokens;
					bool bFound = false;
					for(int i=0 ; ((i < m_activeStateTable->m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1)) ; ++i) {
						if (wgToken[i].iHistoryItem == activeState->iHistoryItem) {
							bFound = true;
						}
					}
					assert(bFound);
				}
			
				++activeState;
			}
		}
		
		// (3) process epsilon transitions in topological order
		m_activeStateTable->processEpsilonTransitions(vFeatureVector.getData(),&m_fScoreBest);
		
		if (t != mFeatures.getRows()-1) { 
		
			// (4) apply beam pruning
			m_activeStateTable->beamPruning(&m_fScoreBest);	
			
			// move to the next time frame
			m_activeStateTable->nextTimeFrame();
		} else {
			cout << "t=" << setw(5) << t << " bestScore: " << FLT(14,6) << m_fScoreBest << endl;
		}
	}
	
	//m_activeStateTable->printInfo();
	m_activeStateTable->endUtterance();
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	
	BVC_VERB << "decoding time: " << FLT(8,2) << dTimeSeconds << " seconds (RTF: " <<  
		FLT(5,2) << dTimeSeconds/(((float)mFeatures.getRows())/100.0) << ")";	
}

// return the best path
BestPath *WFSADecoder::getBestPath() {

	return m_activeStateTable->getBestPath(m_iFeatureVectors);;
}

// return the hypothesis lattice
HypothesisLattice *WFSADecoder::getHypothesisLattice() {

	return m_activeStateTable->getHypothesisLattice();
}

};	// end-of-namespace









