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


#include "ActiveStateTable.h"
#include "BestPath.h"
#include "LMManager.h"
#include "LexiconManager.h"
#include "TimeUtils.h"

namespace Bavieca {

// contructor
ActiveStateTable::ActiveStateTable(float fPruningLikelihood, int iPruningMaxStates, unsigned int iActiveStatesMax, PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMStateDecoding *hmmStatesDecoding, bool bLatticeGeneration, int iMaxWordSequencesState)
{
	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
	m_hmmStatesDecoding = hmmStatesDecoding;	

	m_fPruningLikelihood = fPruningLikelihood;
	m_iPruningMaxStates = iPruningMaxStates;
	m_iActiveStatesMax = iActiveStatesMax;
	m_iHashBuckets = 2*iActiveStatesMax;							
	m_iHashEntries = m_iHashBuckets+m_iActiveStatesMax;		// this guarantees that in the worst case scenario (all collisions), there are enough entries to store collisions
		
	// lattice generation
	m_bLatticeGeneration = bLatticeGeneration;
	m_iMaxWordSequencesState = iMaxWordSequencesState;
	m_iWSHashBuckets = UINT_MAX;
	m_iWSHashEntries = UINT_MAX;	
	m_wshashEntries = NULL;
	m_iWSHashEntryCollisionAvailable = -1;
	m_iLexUnitPronSilence = m_lexiconManager->getLexUnitSilence()->iLexUnitPron;
	m_iLexUnitPronUnknown = m_lexiconManager->m_lexUnitUnknown->iLexUnitPron;
		
	m_activeStatesCurrent = NULL;
	m_activeStatesNext = NULL;
	m_hashEntries = NULL;
	
	// history items
	m_iHistoryItems = 0;
	m_historyItems = NULL;
	m_iHistoryItemAvailable = -1;
	m_iTimeGarbageCollectionLast = -1;
	
	// wrod-graph tokens
	m_iWGTokens = 0;
	m_wgTokens = NULL;
	m_iWGTokenAvailable = -1;
	m_iWGTokensUsed = 0;
	m_iWordSequenceInitial = -1;
}

ActiveStateTable::~ActiveStateTable()
{
	if (m_activeStatesCurrent) {
		delete [] m_activeStatesCurrent;
	}
	if (m_activeStatesNext) {
		delete [] m_activeStatesNext;
	}
	if (m_activeStatesEpsilon) {
		delete [] m_activeStatesEpsilon;
	}	
	if (m_hashEntries) {	
		delete [] m_hashEntries;
	}
	if (m_historyItems) {
		delete [] m_historyItems;
	}
	if (m_wshashEntries) {
		cleanHashWordSequences();	
		delete [] m_wshashEntries;
	}
	if (m_wgTokens) {
		delete [] m_wgTokens;
	}
}

// hash table related functions
void ActiveStateTable::initialize() {
	
	// allocate memory for the hash table
	m_hashEntries = new HashEntry[m_iHashEntries];
	// set all the entries as "old"
	for(int i=0 ; i < (int)m_iHashEntries ; ++i) {
		m_hashEntries[i].iTime = -1;	
	}
	
	m_iTimeCurrent = 0;
	
	// allocate memory for the hash tables of active states
	m_activeStatesCurrent = new ActiveState[m_iActiveStatesMax];
	m_activeStatesNext = new ActiveState[m_iActiveStatesMax];
	m_activeStatesEpsilon = new ActiveStateEpsilon[m_iActiveStatesMax];
	m_activeStateAvailable = m_activeStatesNext;
	m_activeStateEpsilonHead = m_activeStatesEpsilon;
	m_activeStateEpsilonTail = m_activeStatesEpsilon;	
	m_iHashEntryCollisionAvailable = m_iHashBuckets;
	
	m_iActiveStatesCurrent = 0;
	
	// history item management
	//m_iHistoryItems = m_iActiveStatesMax; 
	m_iHistoryItems = 100; 
	//m_iHistoryItems = 10000000; 
	m_historyItems = new HistoryItem[m_iHistoryItems];
	for(unsigned int i=0 ; i < m_iHistoryItems-1 ; ++i) {
		m_historyItems[i].iActive = -1;
		m_historyItems[i].iPrev = i+1;
	}
	m_historyItems[m_iHistoryItems-1].iActive = -1;
	m_historyItems[m_iHistoryItems-1].iPrev = -1;
	m_iHistoryItemAvailable = 0;
	m_iTimeGarbageCollectionLast = -1;
	
	// lattice generation
	if (m_bLatticeGeneration) {
		// allocate memory for the hash table
		m_iWSHashBuckets = 100000;
		m_iWSHashEntries = 2*m_iWSHashBuckets;
		m_wshashEntries = new WSHashEntry[m_iWSHashEntries];
		// set all the entries as "old"
		for(int i=0 ; i < (int)m_iWSHashEntries ; ++i) {
			m_wshashEntries[i].iTime = -1;
		}	
		m_iWSHashEntryCollisionAvailable = m_iWSHashBuckets;	
		
		// lattice token management
		m_iWGTokens = m_iActiveStatesMax*m_iMaxWordSequencesState; 
		//m_iWGTokens = 30*m_iMaxWordSequencesState; 
		//m_iWGTokens = 10000000;
		m_wgTokens = new WGToken[m_iWGTokens];
		for(unsigned int i=0 ; i < m_iWGTokens-1 ; i += m_iMaxWordSequencesState) {
			m_wgTokens[i].iActive = -1;
			m_wgTokens[i].iPrev = i+m_iMaxWordSequencesState;
		}	
		m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iActive = -1;
		m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iPrev = -1;
		m_iWGTokenAvailable = 0;
	}
}


// compute the load factor of the hash table containing active states
float ActiveStateTable::computeLoadFactor() {

	// count 
	unsigned int iValid = 0;
	for(int i=0 ; i < (int)m_iHashBuckets ; ++i) {
		if (m_hashEntries[i].iTime == m_iTimeCurrent) {
			++iValid;
		}
	}
	
	return ((float)iValid)/((float)m_iHashBuckets);
}

// compute the collision factor (# collisions / # valid entries in the table)
float ActiveStateTable::computeCollisionFactor() {

	// count 
	unsigned int iValid = 0;
	unsigned int iCollisions = 0;
	for(unsigned int i=0 ; i < m_iHashEntries ; ++i) {
		if (m_hashEntries[i].iTime == m_iTimeCurrent) {
			if (i < m_iHashBuckets) {
				++iValid;
			} else {
				++iCollisions;
			}
		}
	}
	
	if (iValid == 0) {
		assert(iCollisions == 0);
		return 0.0;
	}
	
	return ((float)iCollisions)/((float)iValid);
}

// compute the load factor of the hash table containing unque word sequences
float ActiveStateTable::computeLoadFactorHashWordSequences(int *iBucketsUsed, int *iCollisions) {

	// count 
	*iBucketsUsed = 0;
	*iCollisions = 0;
	for(unsigned int i=0 ; i < m_iWSHashEntries ; ++i) {
		if (m_wshashEntries[i].iTime > -1) {
			if (i < m_iWSHashBuckets) {
				++(*iBucketsUsed);
			} else {
				++(*iCollisions);
			}
		}
	}

	// load factor
	return ((float)(*iBucketsUsed))/((float)m_iWSHashBuckets);
}

// shows object information
void ActiveStateTable::printInfo() {

	float fLoadFactor = computeLoadFactor();
	float fCollisionFactor = computeCollisionFactor();
	float fSize = m_iActiveStatesMax*3*sizeof(ActiveState);			// there are three tables of active states
	fSize += m_iHashEntries*sizeof(HashEntry);							// there is a hash table
	fSize += m_iHistoryItems*sizeof(HistoryItem);
	// convert to MB
	fSize /= 1024*1024;

	printf("-----------------------------------------\n");
	printf("load factor:        %8.4f\n",fLoadFactor);
	printf("collision factor:   %8.4f\n",fCollisionFactor);
	printf("size:               %8.2f MB\n",fSize);	
	printf("-----------------------------------------\n");
}

// activates the initial state
void ActiveStateTable::activateStateInitial(StateX *state) {

	// create the <s> history item
	int iHistoryItem = newHistoryItem();
	HistoryItem *historyItem = m_historyItems+iHistoryItem;
	historyItem->iLexUnitPron = m_lexiconManager->m_lexUnitBegSentence->iLexUnitPron;
	historyItem->iEndFrame = -1;
	historyItem->iPrev = -1;
	historyItem->iActive = -1;
	historyItem->iWGToken = -1;

	// treat it like an epsilon state
	m_activeStateEpsilonTail->state = state;	
	m_activeStateEpsilonTail->fScore = 0.0;	
	m_activeStateEpsilonTail->iHistoryItem = iHistoryItem;	
	if (m_bLatticeGeneration) {	
		// create the root lattice token
		m_activeStateEpsilonTail->iWGToken = newWGToken();	
		(m_activeStateEpsilonTail->iWGToken+m_wgTokens)[0].iWordSequence = -2;
		(m_activeStateEpsilonTail->iWGToken+m_wgTokens)[0].iLexUnitPron = m_iLexUnitPronUnknown;
		(m_activeStateEpsilonTail->iWGToken+m_wgTokens)[0].fScore = 0.0;
		(m_activeStateEpsilonTail->iWGToken+m_wgTokens)[0].iHistoryItem = iHistoryItem;
		(m_activeStateEpsilonTail->iWGToken+m_wgTokens)[1].iWordSequence = -1;
		m_iWordSequenceInitial = (m_activeStateEpsilonTail->iWGToken+m_wgTokens)[0].iWordSequence;
	} else {
		m_activeStateEpsilonTail->iWGToken = -1;
	}
	++m_activeStateEpsilonTail;
}

// process epsilon transitions in topological order
void ActiveStateTable::processEpsilonTransitions(float *fFeatureVector, float *fScoreBest) {

	TransitionX *transition = NULL;
	TransitionX *transitionEnd = NULL;
	TransitionX *transitionAux = NULL;
	float fScore = 0.0;
	HMMStateDecoding *hmmStateDecoding = NULL;
	
	unsigned int iLexUnitEndSentence = m_lexiconManager->m_lexUnitEndSentence->iLexUnit;	

	assert(m_activeStateEpsilonHead <= m_activeStateEpsilonTail);
	while(m_activeStateEpsilonHead != m_activeStateEpsilonTail) {
	
		//printf("# active epsilon states: %d\n",m_activeStateEpsilonTail-m_activeStateEpsilonHead);
		
		transition = *m_activeStateEpsilonHead->state;
		transitionEnd = *(m_activeStateEpsilonHead->state+1);
		
		//printf("# transitions: %d\n",transitionEnd-transition);
		
		while(transition != transitionEnd) {
		
			//printf("%u %f %x\n",transition->iSymbol,transition->fWeight,transition->state);
			//printTransition(transition);
		
			// epsilon transition
			if (transition->iSymbol & EPSILON_TRANSITION) {
			
				// lattice generation
				int iWGToken = -1;
				if (m_bLatticeGeneration) {
				 	assert(m_activeStateEpsilonHead->iWGToken != -1);
				 	iWGToken = newWGToken(m_activeStateEpsilonHead->iWGToken);
					WGToken *wgToken = iWGToken+m_wgTokens;
					// update the scores
					for(int i=0 ; ((i < m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1)) ; ++i) {
						wgToken[i].fScore += transition->fWeight;
					}
				}	
			
				// activate the state
				activateStateEpsilon(transition->state,m_activeStateEpsilonHead->fScore+transition->fWeight,
					m_activeStateEpsilonHead->iHistoryItem,iWGToken,0.0);
			} 
			// fake transitions
			else if (transition->iSymbol & FAKE_TRANSITION) {
			
			} 
			// lexical unit transition
			else if (transition->iSymbol & LEX_UNIT_TRANSITION) {
			
				if ((transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == iLexUnitEndSentence) {
					++transition;
					continue;
				}
			
				// create a new history item
				int iHistoryItem = newHistoryItem();
				HistoryItem *historyItem = m_historyItems+iHistoryItem;
				historyItem->iPrev = m_activeStateEpsilonHead->iHistoryItem;
				historyItem->iLexUnitPron = transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				historyItem->iEndFrame = m_iTimeCurrent-1;
				historyItem->fScore = m_activeStateEpsilonHead->fScore+transition->fWeight;
				historyItem->iActive = -1;
				
				// lattice generation
				int iWGToken = -1;
				int iWordSequence = -1;
				if (m_bLatticeGeneration) {
					// keep the best N word sequences arriving to the state
					assert(m_activeStateEpsilonHead->iWGToken != -1);
					historyItem->iWGToken = m_activeStateEpsilonHead->iWGToken;
					// checks
					for(int i=0 ; i < m_iMaxWordSequencesState ; ++i) {
						if ((historyItem->iWGToken+m_wgTokens)[i].iWordSequence == -1) {
							break;
						}
						m_lexiconManager->print(m_lexiconManager->getLexUnitPron(historyItem->iLexUnitPron));
						assert(historyItem->iEndFrame > m_historyItems[(historyItem->iWGToken+m_wgTokens)[i].iHistoryItem].iEndFrame);
					}	
					// generate a new hash value for the new word sequence
					iWordSequence = hashWordSequence(historyItem);	
				}
			
				for(transitionAux = *(transition->state) ; *(transition->state+1) ; ++transitionAux) {
					
					if ((transitionAux->iSymbol & LEX_UNIT_TRANSITION) == iLexUnitEndSentence) {
						continue;
					} else {
						//printTransition(transitionAux);
						assert((transitionAux->iSymbol & LEX_UNIT_TRANSITION) == 0);
					}
					
					// epsilon-transition
					if (transitionAux->iSymbol & EPSILON_TRANSITION) {
					
						// preventive pruning
						if (m_activeStateEpsilonHead->fScore+transition->fWeight+transitionAux->fWeight < (*fScoreBest-m_fPruningLikelihood)) {
							continue;
						}
						
						// lattice generation
						if (m_bLatticeGeneration) {
							iWGToken = newWGToken(iWordSequence,m_activeStateEpsilonHead->fScore+transition->fWeight+transitionAux->fWeight,iHistoryItem);
						}
					
						// activate the state
						activateStateEpsilon(transitionAux->state,
							m_activeStateEpsilonHead->fScore+transition->fWeight+transitionAux->fWeight,iHistoryItem,iWGToken,0.0);
					} 
					// fake transition
					else if (transitionAux->iSymbol & FAKE_TRANSITION) {
					
					} 
					// leaf-transition
					else {
					
						// compute emission probability
						hmmStateDecoding = &m_hmmStatesDecoding[transitionAux->iSymbol];
					#ifdef SIMD	
						fScore = hmmStateDecoding->computeEmissionProbabilityNearestNeighborSIMD(fFeatureVector,m_iTimeCurrent);	
					#else
						fScore = hmmStateDecoding->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,m_iTimeCurrent);	
					#endif
						
						// preventive pruning
						if (m_activeStateEpsilonHead->fScore+transition->fWeight+transitionAux->fWeight+fScore < (*fScoreBest-m_fPruningLikelihood)) {
							continue;
						}	
						
						// lattice generation
						if (m_bLatticeGeneration) {
							iWGToken = newWGToken(iWordSequence,m_activeStateEpsilonHead->fScore+transition->fWeight+transitionAux->fWeight+fScore,iHistoryItem);
						}	
						
						// activate the state	
						activateState(transitionAux->state,
							m_activeStateEpsilonHead->fScore+transition->fWeight+transitionAux->fWeight+fScore,
							fScoreBest,hmmStateDecoding,iHistoryItem,iWGToken,0.0);
					}
				}
			} 
			// leaf transition
			else {
			
				// compute emission probability
				hmmStateDecoding = &m_hmmStatesDecoding[transition->iSymbol];
			#ifdef SIMD
				float fScore = hmmStateDecoding->computeEmissionProbabilityNearestNeighborSIMD(fFeatureVector,0);
			#else
				float fScore = hmmStateDecoding->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,0);
			#endif
			
				// preventive pruning goes here
				
				// lattice generation
				int iWGToken = -1;
				if (m_bLatticeGeneration) {
				 	assert(m_activeStateEpsilonHead->iWGToken != -1);
				 	iWGToken = newWGToken(m_activeStateEpsilonHead->iWGToken);	
					WGToken *wgToken = iWGToken+m_wgTokens;
					// update the scores
					for(int i=0 ; ((i < m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1)) ; ++i) {
						wgToken[i].fScore += transition->fWeight+fScore;
					} 	
				}
			
				// activate the state	
				activateState(transition->state,m_activeStateEpsilonHead->fScore+transition->fWeight+fScore,fScoreBest,
					hmmStateDecoding,m_activeStateEpsilonHead->iHistoryItem,iWGToken,0.0);
			}
				
			++transition;
		}
			
		++m_activeStateEpsilonHead;
		if (m_activeStateEpsilonHead-m_activeStatesEpsilon == (int)m_iActiveStatesMax) {
			m_activeStateEpsilonHead = m_activeStatesEpsilon;
		}
	}
}

// apply beam pruning to the set of active states
void ActiveStateTable::beamPruning(float *fScoreBest) {

	unsigned int iActive = m_activeStateAvailable-m_activeStatesNext; 
	unsigned int iRemoved = 0;
	
	int iNumberBins = BEAM_PRUNING_NUMBER_BINS;
	float fScoreWorst = FLT_MAX;
	
	//checkTokens();
	
	//printf("t=%d pruning starts -------------------------------------------\n",m_iTimeCurrent);
	
	// (1) standard beam-pruning O(n)
	float fMinimumScore = *fScoreBest-m_fPruningLikelihood;
	for(ActiveState *activeState = m_activeStatesNext ; activeState != m_activeStateAvailable ; ++activeState) {
		if (activeState->fScore < fMinimumScore) {
			activeState->state = NULL;
			++iRemoved;
		} else if (activeState->fScore < fScoreWorst) {
			fScoreWorst = activeState->fScore;
		}
		//printf("%10.2f\n",activeState->fScore);
	}
	
   // (2) Histogram pruning O(n)
   
   // check if it applies (only if the number of active states exceeds the maximum allowed)
   if (iActive-iRemoved > (unsigned int)m_iPruningMaxStates) {
   
		// (2.1) compute the size of each bin and initialize them
		float fLength = *fScoreBest-fScoreWorst+1;
		float fBinSize = ((float)fLength)/((float)iNumberBins);
		assert(fBinSize > 0);
		int iBins[iNumberBins];
		for(int i = 0 ; i < iNumberBins ; ++i) {
			iBins[i] = 0;
		}
		// (2.2) fill the bins, the first bin keeps the best tokens
		int iBin;
		float fAux = ((float)iNumberBins)/fLength;
		for(ActiveState *activeState = m_activeStatesNext ; activeState != m_activeStateAvailable ; ++activeState) {
			if (activeState->state != NULL) {
				iBin = (int)(fabs(*fScoreBest-activeState->fScore)*fAux);
				assert((iBin >= 0) && (iBin < iNumberBins));
				iBins[iBin]++;	
			}
		}
		//for(int i = 0 ; i < iNumberBins ; ++i) {
		//	printf("bin %3d: %d\n",i,iBins[i]);
		//}		
		// (2.3) actual pruning O(n), tokens in the last bin always survive
		int iSurvivors = 0;
		for(int i = 0 ; i < iNumberBins-1 ; ++i) {
			iSurvivors += iBins[i];
			// remove tokens in the remaining bins
			if (iSurvivors >= m_iPruningMaxStates) {
				float fThreshold = *fScoreBest-(((float)(i+1))*(fLength/((float)iNumberBins)));
				for(ActiveState *activeState = m_activeStatesNext ; activeState != m_activeStateAvailable ; ++activeState) {
					if (activeState->state != NULL) {
						if (activeState->fScore < fThreshold) {
							activeState->state = NULL;
							++iRemoved;
						}
					}
				}
				break;
			}
		}		
	}
	
	if (m_iTimeCurrent%100 == 0) {
		printf("t=%5d # active states: %6u (removed: %6u, bestScore: %14.6f)\n",m_iTimeCurrent,iActive-iRemoved,iRemoved,*fScoreBest);
		if (m_bLatticeGeneration) {
			//printHashTableWordSequences();
		}
	}
	
   m_iStatesPruned = iRemoved;   
   
   //printf("expanded: %10u active: %10u pruned: %10u\n",m_iStatesExpanded,m_iStatesActivated,m_iStatesPruned);
}

// garbage collection of history items and word-lattice tokens (lattice generation mode)
// (1) it starts by marking the active items by traversing back items from the active states
// (2) it adds inactive items to the queue of available items

// note: garbage collection of lattice tokens is necessary since every time there is 
// a merge of two list of tokens (only the n-best tokens are kept at any state) the tree of
// lattice tokens belonging to the removed elements is lost
void ActiveStateTable::historyItemGarbageCollection() {
	
	unsigned int iItemsActive = 0;
	unsigned int iTokensActive = 0;
	
	// (1) check if garbage collection was already run within the current time frame
	// note: this is an undesirable situation because it requires an extra pass over the complete array of items
	// it should be avoided by allocating a larger number of entries from the beginning
	if (m_iTimeGarbageCollectionLast == m_iTimeCurrent) {
		// mark all the history items and wg-tokens as inactive
		for(unsigned int i=0 ; i < m_iHistoryItems ; ++i) {
			m_historyItems[i].iActive = -1;
		}
		for(unsigned int i=0 ; i < m_iWGTokens ; i += m_iMaxWordSequencesState) {
			m_wgTokens[i].iActive = -1;
		}	
	}
	
	bool *bTokensActive = new bool[m_iWGTokens];
	for(unsigned int i=0 ; i < m_iWGTokens ; ++i) {
		bTokensActive[i] = false;
	}
	
	// array to keep the active history items
	int iHistoryItemActiveSize = 0;
	int *iHistoryItemActive = new int[m_iHistoryItems];	
	
	// (1) mark items coming from active nodes as active
	// (1.1) active nodes from current time frame
	ActiveState *activeState = m_activeStatesCurrent;
	ActiveState *activeStateCurrentLast = m_activeStatesCurrent+m_iActiveStatesCurrent;	
	while(activeState != activeStateCurrentLast) {
	
		// skip pruned states
		if (activeState->state == NULL) {
			++activeState;
			continue;
		}	
		int iHistoryItem = activeState->iHistoryItem;
		while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
			(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;	
			iHistoryItemActive[iHistoryItemActiveSize++] = iHistoryItem;
			iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
			++iItemsActive;
		}
		
		if (m_bLatticeGeneration) {
			// mark the wg-token as active
			assert(activeState->iWGToken != -1);
			(activeState->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
			assert((bTokensActive[activeState->iWGToken] == false));
			bTokensActive[activeState->iWGToken] = true;
			++iTokensActive;
			
			// keep history items at the current wg-token
			for(int i=0 ; (i < m_iMaxWordSequencesState) && ((activeState->iWGToken+m_wgTokens)[i].iWordSequence != -1) ; ++i) {
				HistoryItem *historyItem = m_historyItems+(activeState->iWGToken+m_wgTokens)[i].iHistoryItem;
				if (historyItem->iActive != m_iTimeCurrent) {
					historyItem->iActive = m_iTimeCurrent;
					assert(historyItem->iWGToken >= -1);
					iHistoryItemActive[iHistoryItemActiveSize++] = historyItem-m_historyItems;
					++iItemsActive;
				}
			}	
		}
			
		++activeState;
	}	
	
	// (1.2) epsilon states	
	assert(m_activeStateEpsilonHead <= m_activeStateEpsilonTail);
	ActiveStateEpsilon *activeStateEpsilon = m_activeStateEpsilonHead;
	while(activeStateEpsilon != m_activeStateEpsilonTail) {
	
		assert(activeStateEpsilon->state != NULL);	
		
		int iHistoryItem = activeStateEpsilon->iHistoryItem;
		while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
			(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;	
			iHistoryItemActive[iHistoryItemActiveSize++] = iHistoryItem;
			iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
			++iItemsActive;
		}
		
		if (m_bLatticeGeneration) {
			// mark the wg-token as active
			assert(activeStateEpsilon->iWGToken != -1);
			(activeStateEpsilon->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
			assert((bTokensActive[activeStateEpsilon->iWGToken] == false));
			bTokensActive[activeStateEpsilon->iWGToken] = true;
			++iTokensActive;
			
			// keep history items at the current wg-token
			for(int i=0 ; (i < m_iMaxWordSequencesState) && ((activeStateEpsilon->iWGToken+m_wgTokens)[i].iWordSequence != -1) ; ++i) {
				HistoryItem *historyItem = m_historyItems+(activeStateEpsilon->iWGToken+m_wgTokens)[i].iHistoryItem;
				if (historyItem->iActive != m_iTimeCurrent) {
					historyItem->iActive = m_iTimeCurrent;
					assert(historyItem->iWGToken >= -1);
					iHistoryItemActive[iHistoryItemActiveSize++] = historyItem-m_historyItems;
					++iItemsActive;
				}
			}
		}
		
		++activeStateEpsilon;
		
		if (activeStateEpsilon-m_activeStatesEpsilon == (int)m_iActiveStatesMax) {
			activeStateEpsilon = m_activeStatesEpsilon;
		}
	}	
	
	// (1.3) active nodes from next time frame
	activeState = m_activeStatesNext;
	while(activeState != m_activeStateAvailable) {
	
		assert(activeState->state != NULL);
	
		int iHistoryItem = activeState->iHistoryItem;
		while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
			(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;	
			iHistoryItemActive[iHistoryItemActiveSize++] = iHistoryItem;
			iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
			++iItemsActive;
		}
		
		if (m_bLatticeGeneration) {
			// mark the wg-token as active
			assert(activeState->iWGToken != -1);
			(activeState->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
			assert((bTokensActive[activeState->iWGToken] == false));
			bTokensActive[activeState->iWGToken] = true;
			++iTokensActive;
			
			// keep history items at the current wg-token
			for(int i=0 ; (i < m_iMaxWordSequencesState) && ((activeState->iWGToken+m_wgTokens)[i].iWordSequence != -1) ; ++i) {
				HistoryItem *historyItem = m_historyItems+(activeState->iWGToken+m_wgTokens)[i].iHistoryItem;
				if (historyItem->iActive != m_iTimeCurrent) {
					historyItem->iActive = m_iTimeCurrent;
					assert(historyItem->iWGToken >= -1);
					iHistoryItemActive[iHistoryItemActiveSize++] = historyItem-m_historyItems;
					++iItemsActive;
				}
			}	
		}	
		++activeState;
	}
	
	if (m_bLatticeGeneration) {
		while(iHistoryItemActiveSize > 0) {
		
			HistoryItem *historyItem = m_historyItems+iHistoryItemActive[iHistoryItemActiveSize-1];
			--iHistoryItemActiveSize;
			
			if (historyItem->iWGToken != -1) {				
				(historyItem->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;				
				if (bTokensActive[historyItem->iWGToken] == false) {
					bTokensActive[historyItem->iWGToken] = true;
					++iTokensActive;	
				
					for(int i=0 ; ((i < m_iMaxWordSequencesState) && ((historyItem->iWGToken+m_wgTokens)[i].iWordSequence != -1)) ; ++i) {					
						HistoryItem *historyItem2 = m_historyItems+(historyItem->iWGToken+m_wgTokens)[i].iHistoryItem;	
						if (historyItem2->iActive != m_iTimeCurrent) {
							historyItem2->iActive = m_iTimeCurrent;
							iHistoryItemActive[iHistoryItemActiveSize++] = historyItem2-m_historyItems;	
							++iItemsActive;
						}
					}
				}
			}	
		}	
	}
	
	//printf("time: %d\n",m_iTimeCurrent);
	//printf("history item garbage collection...\n");	
	//printf("(items used: %d existing: %d)\n",iItemsActive,m_iHistoryItems);	
	
	// a bigger data structure to keep the new history items
	// (* if we wait until all the items are active there will be many calls to the garbage collector
	// when the array of reach a high occupation, that would introduce substantial overhead)
	assert(iItemsActive <= m_iHistoryItems);
	if (iItemsActive >= (0.20*m_iHistoryItems)) {

		//printf("history item garbage collection...\n");	
		//printf("allocating space for new items (item sused: %d existing: %d)\n",iItemsActive,m_iHistoryItems);	
		
		// allocate a new data structure with double capacity	
		HistoryItem *historyItems = new HistoryItem[m_iHistoryItems*2];
		
		// copy the active items from the old data structure
		for(unsigned int i=0 ; i < m_iHistoryItems ; ++i) {
			memcpy(historyItems+i,m_historyItems+i,sizeof(HistoryItem));
			historyItems[i].iActive = m_iTimeCurrent;
		}
		
		// create the linked list of available items
		for(unsigned int i=m_iHistoryItems ; i < (2*m_iHistoryItems)-1 ; ++i) {
			historyItems[i].iPrev = i+1;
			historyItems[i].iActive = -1;	
		}
		historyItems[(2*m_iHistoryItems)-1].iPrev = -1;
		historyItems[(2*m_iHistoryItems)-1].iActive = -1;
		
		delete [] m_historyItems;
		m_historyItems = historyItems;
		m_iHistoryItemAvailable = m_iHistoryItems;
		m_iHistoryItems *= 2;
	}
	// (2') there are inactive items: create a linked list with them
	else {
		int *iHistoryItemAux = &m_iHistoryItemAvailable;
		for(unsigned int i = 0 ; i < m_iHistoryItems ; ++i) {
			if (m_historyItems[i].iActive != m_iTimeCurrent) {
				m_historyItems[i].iActive = -1;
				m_historyItems[i].iEndFrame = -1;
				*iHistoryItemAux = i;
				iHistoryItemAux = &m_historyItems[i].iPrev;	
			}
		}
		*iHistoryItemAux = -1;
	}
	
	delete [] iHistoryItemActive;	
	
	//printf("wg-token garbage collection...\n");	
	//printf("(wg-tokens used: %d existing: %d)\n",iTokensActive,m_iWGTokens/m_iMaxWordSequencesState);	
	
	// lattice token garbage collection
	if (m_bLatticeGeneration) {
		assert(iTokensActive <= (m_iWGTokens/m_iMaxWordSequencesState));
		if (iTokensActive >= 0.20*(m_iWGTokens/m_iMaxWordSequencesState)) {
		
			//assert(0);
			
			// allocate a new data structure with double capacity	
			WGToken *wgTokens = NULL;
			try {
				//printf("allocating wg-tokens: %d\n",m_iWGTokens*2);
				wgTokens = new WGToken[m_iWGTokens*2];
			}
			catch (const std::bad_alloc&) {
				int iMB = m_iWGTokens*2*sizeof(WGToken);
				BVC_ERROR << "unable to allocate memory for the lattice tokens, " << iMB << "MBs needed";
			}		
			//printf("WGTokens: from %d to %d\n",m_iWGTokens*sizeof(WGToken),m_iWGTokens*2*sizeof(WGToken));
			
			// copy the active items from the old data structure
			for(unsigned int i=0 ; i < m_iWGTokens ; ++i) {
				memcpy(wgTokens+i,m_wgTokens+i,sizeof(WGToken));
				wgTokens[i].iPrev = -1;
				wgTokens[i].iActive = m_iTimeCurrent;
			}
			// create the linked list of available items
			unsigned int iLastIndex = ((2*m_iWGTokens)-m_iMaxWordSequencesState);
			for(unsigned int i=m_iWGTokens ; i < iLastIndex ; i += m_iMaxWordSequencesState) {
				wgTokens[i].iPrev = i+m_iMaxWordSequencesState;
				wgTokens[i].iActive = -1;	
			}
			wgTokens[iLastIndex].iPrev = -1;
			wgTokens[iLastIndex].iActive = -1;
			
			// swap structures and delete the old one
			delete [] m_wgTokens;
			m_wgTokens = wgTokens;
			m_iWGTokenAvailable = m_iWGTokens;
			m_iWGTokens *= 2;
		}
		// there are inactive tokens: create a linked list with them
		else {
			//printf("WGTokens: relinking\n");
			// the linked list goes from lower to higher memory addresses
			int *iAux = &m_iWGTokenAvailable;
			for(unsigned int i = 0 ; i < m_iWGTokens ; i += m_iMaxWordSequencesState) {
				if (m_wgTokens[i].iActive != m_iTimeCurrent) {	
					*iAux = i;
					iAux = &m_wgTokens[i].iPrev;	
				}	
			}
			*iAux = -1;
		}
	}
		
	delete [] bTokensActive;
	
	m_iTimeGarbageCollectionLast = m_iTimeCurrent;
}

// recovers the best path from the list of active states
BestPath *ActiveStateTable::getBestPath(int iFeatureVectors) {

	BestPath *bestPath = new BestPath(m_lexiconManager,-1.0);
	
	// (1) get the best scoring active state
	ActiveState *activeStateBest = NULL;
	int iLexUnitLastBest = -1;
	int iLexUnitLast = -1;
	ActiveState *activeState = m_activeStatesNext;
	float fScoreBest = -FLT_MAX;
	while(activeState != m_activeStateAvailable) {
		// get the last lexical unit
		iLexUnitLast = -1;	
		// get the best transition to a lexical unit (the transition with smaller weight including the transition to </s>)
		TransitionX *transitionBest = NULL;
		float fScoreBestTransition = -FLT_MAX;
		for(TransitionX *transition = *activeState->state ; transition != *(activeState->state+1) ; ++transition) {
			if (transition->iSymbol & LEX_UNIT_TRANSITION) {
				float fScore = transition->fWeight;
				// find the transition to the end of sentence </s>
				bool bFound = false;
				float fScoreFinalBest = -FLT_MAX;
				for(TransitionX *transition2 = *transition->state ; transition2 != *(transition->state+1) ; ++transition2) {
					if (transition2->iSymbol & LEX_UNIT_TRANSITION) {
						if ((int)(transition2->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
							if (transition2->fWeight > fScoreFinalBest) {
								fScoreFinalBest = transition2->fWeight;	
							}
							bFound = true;
							//break;
						}
					} else if (transition2->iSymbol & EPSILON_TRANSITION) {
						for(TransitionX *transition3 = *transition2->state ; transition3 != *(transition2->state+1) ; ++transition3) {
							if (transition3->iSymbol & LEX_UNIT_TRANSITION) {
								if ((int)(transition3->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
									if (transition2->fWeight + transition3->fWeight > fScoreFinalBest) {
										fScoreFinalBest = transition2->fWeight + transition3->fWeight;
									}
									bFound = true;
									break;	// there can't be multiple epsilon symbols coming from a given state
								}
							} else if (transition3->iSymbol & EPSILON_TRANSITION) {
								for(TransitionX *transition4 = *transition3->state ; transition4 != *(transition3->state+1) ; ++transition4) {
									if (transition4->iSymbol & LEX_UNIT_TRANSITION) {
										if ((int)(transition4->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
											if (transition2->fWeight + transition3->fWeight + transition4->fWeight > fScoreFinalBest) {
												fScoreFinalBest = transition2->fWeight + transition3->fWeight + transition4->fWeight;
											}
											bFound = true;
											break;	// there can't be multiple epsilon symbols coming from a given state
										}
									}
								}
							}
						}
						/*if (bFound) {
							break;	
						}*/
					}
				}
				assert(bFound);
				fScore += fScoreFinalBest;
				if ((transitionBest == NULL) || (fScore > fScoreBestTransition)) {
					transitionBest = transition;
					fScoreBestTransition = fScore;
				}
			} else if (transition->iSymbol & EPSILON_TRANSITION) {
				assert(0);
			}
		}
		// only consider states that go to a lexical unit (end-of-word states)
		if (transitionBest != NULL) {
			
			activeState->fScore += fScoreBestTransition;
			iLexUnitLast = transitionBest->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;

			if (activeState->fScore > fScoreBest) {
				activeStateBest = activeState;
				iLexUnitLastBest = iLexUnitLast;
				fScoreBest = activeState->fScore;
			}	
		}
		++activeState;
	}
	if (activeStateBest == NULL) {
		return NULL;
	}
	
	bestPath->setScore(fScoreBest);
	
	// (2) get the sequence of lexical units from the best scoring active state
	int iHistoryItem = activeStateBest->iHistoryItem;
	int iFrameStart = 0;
	while (iHistoryItem != -1) {
		HistoryItem *historyItem = m_historyItems+iHistoryItem;
		if (historyItem->iPrev != -1) {
			iFrameStart = (m_historyItems+historyItem->iPrev)->iEndFrame+1;
		} else {
			iFrameStart = 0;
		}
		// get the observed lexical unit
		LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(historyItem->iLexUnitPron);
		// add a new element
		bestPath->newElementFront(iFrameStart,historyItem->iEndFrame,0.0,0.0,0.0,0.0,lexUnit,0.0);
		//printf("%4d %4d %-20s \n",iFrameStart,historyItem->iEndFrame,m_lexiconManager->getStrLexUnit(historyItem->iLexUnit));	
		iHistoryItem = historyItem->iPrev;
	}
	
	// add the final lexical unit if any
	if (iLexUnitLastBest != -1) {
		LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(iLexUnitLastBest);
		if (activeStateBest->iHistoryItem != -1) {
			iFrameStart = (m_historyItems+activeStateBest->iHistoryItem)->iEndFrame+1;
		} else {
			iFrameStart = 0;
		}
		bestPath->newElementBack(iFrameStart,iFeatureVectors-1,0.0,0.0,0.0,0.0,lexUnit,0.0);
	}
	
	// add the end of sentence
	bestPath->newElementBack(iFeatureVectors,iFeatureVectors,0.0,0.0,0.0,0.0,m_lexiconManager->m_lexUnitEndSentence,0.0);
	
	//bestPath->print(false);
	
	//printf("best scoring token: %f\n",activeStateBest->fScore);
	
	return bestPath;
}

// print the given best paths (useful for debugging)
void ActiveStateTable::printBestPaths(unsigned int iPaths) {

	LActiveState lActiveState;

	// put all the active states in a list
	ActiveState *activeState = m_activeStatesNext;
	while(activeState != m_activeStateAvailable) {
		lActiveState.push_back(activeState);
		++activeState;
	}
	
	// sort the list by score
	lActiveState.sort(ActiveStateTable::activeStateComparisonScore);	
	
	// print the active states with their corresponding histories
	unsigned int i = 0;
	for(LActiveState::iterator it = lActiveState.begin() ; it != lActiveState.end() ; ++it, ++i) {
		
		if (i == iPaths) {
			break;
		}	
		
		printf("%4u %.4f ",i,(*it)->fScore);
		
		int iHistoryItem = (*it)->iHistoryItem;
		while(iHistoryItem != -1) {
			printf("%s ",m_lexiconManager->getStrLexUnitPron((m_historyItems+iHistoryItem)->iLexUnitPron));	
			iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
		}
		
		// does the state go to any lex-unit transition?
		int iLexUnitFinal = -1;
		for(TransitionX *transition = *(*it)->state ; transition != *((*it)->state+1) ; ++transition) {
			if ((transition->iSymbol & LEX_UNIT_TRANSITION) != 0) {
				iLexUnitFinal = transition->iSymbol-LEX_UNIT_TRANSITION;
				printf("%s ",m_lexiconManager->getStrLexUnitPron(iLexUnitFinal));
				break;
			}
		}
		printf("\n");		
	}	
}

// print stats connected to the hash table containing unique word-sequences
void ActiveStateTable::printHashTableWordSequences() {

	int iBucketsUsed = 0;
	int iCollisions = 0;
	float fLoadFactor = computeLoadFactorHashWordSequences(&iBucketsUsed,&iCollisions);
	
	printf("- hash table containing word-sequences ---\n");
	printf(" # buckets:      %8d\n",m_iWSHashBuckets);
	printf(" # buckets used: %8d\n",iBucketsUsed);
	printf(" # collisions:   %8d\n",iCollisions);
	printf(" # elements:     %8d\n",iBucketsUsed+iCollisions);
	printf(" load factor: %8.4f\n",fLoadFactor);
	printf("------------------------------------------\n");
}

// build a hypothesis lattice for the utterance
HypothesisLattice *ActiveStateTable::getHypothesisLattice() {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	VHistoryItem vHistoryItem;

	// (1) create the list of history items from the active states 
	// note: active states that are final must go to a lexical-unit transition	
	int iLexUnitLast = -1;
	ActiveState *activeState = m_activeStatesNext;
	while(activeState != m_activeStateAvailable) {
		// get the last lexical unit
		iLexUnitLast = -1;	
		// get the best transition to a lexical unit (the transition with smaller weight including the transition to </s>)
		for(TransitionX *transition = *activeState->state ; transition != *(activeState->state+1) ; ++transition) {
			if (transition->iSymbol & LEX_UNIT_TRANSITION) {
				float fScore = transition->fWeight;
				// find the transition to the end of sentence </s>
				bool bFound = false;
				float fScoreFinalBest = -FLT_MAX;
				for(TransitionX *transition2 = *transition->state ; transition2 != *(transition->state+1) ; ++transition2) {
					// direct transition to: </s>
					if (transition2->iSymbol & LEX_UNIT_TRANSITION) {
						if ((int)(transition2->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
							if (transition2->fWeight > fScoreFinalBest) {
								fScoreFinalBest = transition2->fWeight;	
							}
							bFound = true;
						}
					}
					// epsilon transition before transition to: </s> 
					else if (transition2->iSymbol & EPSILON_TRANSITION) {
						for(TransitionX *transition3 = *transition2->state ; transition3 != *(transition2->state+1) ; ++transition3) {
							if (transition3->iSymbol & LEX_UNIT_TRANSITION) {
								if ((int)(transition3->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
									if (transition2->fWeight + transition3->fWeight > fScoreFinalBest) {
										fScoreFinalBest = transition2->fWeight + transition3->fWeight;
									}
									bFound = true;
									break;	// there can't be multiple epsilon symbols coming from a given state
								}
							} else if (transition3->iSymbol & EPSILON_TRANSITION) {
								for(TransitionX *transition4 = *transition3->state ; transition4 != *(transition3->state+1) ; ++transition4) {
									if (transition4->iSymbol & LEX_UNIT_TRANSITION) {
										if ((int)(transition4->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT) == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
											if (transition2->fWeight + transition3->fWeight + transition4->fWeight > fScoreFinalBest) {
												fScoreFinalBest = transition2->fWeight + transition3->fWeight + transition4->fWeight;
											}
											bFound = true;
											break;	// there can't be multiple epsilon symbols coming from a given state
										}
									}
								}
							}
						}
					}
				}
				assert(bFound);
				
				HistoryItem *historyItem = new HistoryItem();
				historyItem->iLexUnitPron = transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				historyItem->iEndFrame = m_iTimeCurrent;
				historyItem->fScore = activeState->fScore + fScore;
				historyItem->iPrev = activeState->iHistoryItem;
				historyItem->iActive = m_iTimeCurrent;
				historyItem->iWGToken = activeState->iWGToken;
				vHistoryItem.push_back(historyItem);
			} 
			assert((transition->iSymbol & EPSILON_TRANSITION) == 0);
		}
		++activeState;
	}

	// if the vector is empty the lattice cannot be built
	if (vHistoryItem.empty()) {
		return NULL;
	}
	
	// create the initial node 
	LNode *lnodeInitial = HypothesisLattice::newNode(-1);
	
	// create the final node 
	LNode *lnodeFinal = HypothesisLattice::newNode(m_iTimeCurrent);
	
	MHistoryItemLNode mHistoryItemLNode;	
	int iNodes = 0;
	int iEdges = 0;
	
	// put all the history items in the list of items to process
	for(VHistoryItem::iterator it = vHistoryItem.begin() ; it != vHistoryItem.end() ; ++it) {
		
		HistoryItem *historyItem = *it;
		
		assert(historyItem != NULL);
		assert(historyItem->iWGToken != -1);
	
		mHistoryItemLNode.insert(MHistoryItemLNode::value_type(historyItem,lnodeFinal));
	}
	
	// process all the history items
	while(vHistoryItem.empty() == false) {
		
		HistoryItem *historyItemAux = vHistoryItem.back();
		vHistoryItem.pop_back();
		
		// get the graph node for this history item
		MHistoryItemLNode::iterator it = mHistoryItemLNode.find(historyItemAux);
		assert(it != mHistoryItemLNode.end());
		LNode *lnodeAux = it->second;	
	
		//printHistoryItem(historyItemAux);
		
		// items that go to this item
		//printf("----------------------------------------------------\n");
		int iTokens = 0;
		for(int i=0 ; ((i < m_iMaxWordSequencesState) && ((historyItemAux->iWGToken+m_wgTokens)[i].iWordSequence != -1)) ; ++i, ++iTokens) {
			if ((historyItemAux->iWGToken+m_wgTokens)[i].iWordSequence == m_iWordSequenceInitial) {
				if (iTokens == 0) {
					//printf("connects to: <s>\n");
				}
				break;
			}
			//HistoryItem *historyItemPrev = historyItemAux->wgToken[i].historyItem;	
			HistoryItem *historyItemPrev = m_historyItems+(historyItemAux->iWGToken+m_wgTokens)[i].iHistoryItem;	
			LNode *lnodePrev = NULL;
			LEdge *ledgePrev = NULL;
			LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(historyItemAux->iLexUnitPron);
			MHistoryItemLNode::iterator jt = mHistoryItemLNode.find(historyItemPrev);
			// the history item is not in the graph: create a graph node for it
			if (jt == mHistoryItemLNode.end()) {
				//printHistoryItem(historyItemPrev);
				lnodePrev = HypothesisLattice::newNode(historyItemPrev->iEndFrame);
				ledgePrev = HypothesisLattice::newEdge(historyItemPrev->iEndFrame+1,historyItemAux->iEndFrame,lexUnit,0.0,0.0,0.0);
				mHistoryItemLNode.insert(MHistoryItemLNode::value_type(historyItemPrev,lnodePrev));
				vHistoryItem.push_back(historyItemPrev);
				++iNodes;
			}
			// the history item is in the graph: create a link
			else {
				ledgePrev = HypothesisLattice::newEdge(historyItemPrev->iEndFrame+1,historyItemAux->iEndFrame,lexUnit,0.0,0.0,0.0);	
				lnodePrev = jt->second;
			}
			++iEdges;
			// make the connection
			HypothesisLattice::connectEdge(lnodePrev,ledgePrev,lnodeAux);
		}
		// if no tokens then the node should connect to the initial node <s>
		if (iTokens == 0) {	
			// make the connection
			LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(historyItemAux->iLexUnitPron);
			LEdge *ledge = HypothesisLattice::newEdge(0,historyItemAux->iEndFrame,lexUnit,0.0,0.0,0.0);
			// make the connection
			HypothesisLattice::connectEdge(lnodeInitial,ledge,lnodeAux);
		}
		
		// delete only the items that were just created
		if (historyItemAux->iEndFrame == m_iTimeCurrent) {
			delete historyItemAux;
		}
		
		//printf("----------------------------------------------------\n");
 	}
 	
 	assert(lnodeInitial->edgeNext != NULL);
 	
 	HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager);
 	hypothesisLattice->buildContainer(lnodeInitial,lnodeFinal);
 	
 	double dTimeEnd = TimeUtils::getTimeMilliseconds();
 	double dTime = (dTimeEnd-dTimeBegin)/1000.0;
 	
 	printf("Lattice building time: %.4fs\n",dTime);
	
	//exit(-1);

	return hypothesisLattice;
}


// merge two sets of word sequences by keeping the N best unique word sequences in wgToken1 (not commutative)
// 1) sorting: both sets are sorted so it is very efficient: O(n)
// 2) unique: this is linear too
bool ActiveStateTable::mergeWordSequences(int iWGToken1, int iWGToken2) {

	WGToken *wgTokenTable = NULL;
	WGToken *wgToken1 = iWGToken1+m_wgTokens;
	WGToken *wgToken2 = iWGToken2+m_wgTokens;
	
	int iLength1 = 0;
	for( ; ((iLength1 < m_iMaxWordSequencesState) && (wgToken1[iLength1].iWordSequence != -1)) ; ++iLength1);
	int iLength2 = 0;
	for( ; ((iLength2 < m_iMaxWordSequencesState) && (wgToken2[iLength2].iWordSequence != -1)) ; ++iLength2);
	
	assert(iWGToken1 != iWGToken2);	
	assert((iLength1 > 0) && (iLength2 > 0));
		
	int k = iLength1-1;	
	int j = iLength2-1;
	int iTotal = iLength1+iLength2;
	// using a table of elements instead of pointers simplifies the "unique" process and makes it more efficient
	wgTokenTable = new WGToken[iTotal];
	
	// sort all the elements (not just the top N since there can be duplicated word-sequences)
	bool bReturn = false;
	for(int i=iTotal-1 ; i >= 0 ; --i) {
		if (wgToken1[k].fScore < wgToken2[j].fScore) {
			memcpy(&wgTokenTable[i],&wgToken1[k],sizeof(WGToken));
			--k;
			if (k < 0) {
				for(int h=i-1 ; h >= 0 ; --h) {
					assert(j>=0);
					memcpy(&wgTokenTable[h],&wgToken2[j--],sizeof(WGToken));
				}	
				bReturn = false;
				break;
			}
		} else {
			memcpy(&wgTokenTable[i],&wgToken2[j],sizeof(WGToken));
			--j;
			if (j < 0) {
				for(int h=i-1 ; h >= 0 ; --h) {
					assert(k>=0);
					memcpy(&wgTokenTable[h],&wgToken1[k--],sizeof(WGToken));
				}
				bReturn = true;
				break;
			}
		}	
	}
	
	// unique 
	int iUniqueElements = 0;
	memcpy(&wgToken1[0],&wgTokenTable[0],sizeof(WGToken));
	for(int i=1 ; ((i < iTotal) && (iUniqueElements+1 < m_iMaxWordSequencesState)) ; ++i) {
		bool bDuplicated = false;
		for(int j=0 ; j < iUniqueElements+1 ; ++j) {
			if (wgTokenTable[i].iWordSequence == wgToken1[j].iWordSequence) {
				bDuplicated = true;
				break;
			}
		}
		if (bDuplicated == false) {	
			memcpy(&wgToken1[++iUniqueElements],&wgTokenTable[i],sizeof(WGToken));
		}
	}	
	
	if (iUniqueElements+1 < m_iMaxWordSequencesState) {
		wgToken1[iUniqueElements+1].iWordSequence = -1;
	}
	
	delete [] wgTokenTable;

	return bReturn;
}

};	// end-of-namespace
