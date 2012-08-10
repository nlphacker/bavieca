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

#ifndef ACTIVESTATEHASH_H
#define ACTIVESTATEHASH_H

#include <stdio.h>
#include <limits.h>
#include <float.h>

#include "BestPath.h"
#include "HMMStateDecoding.h"
#include "LexiconManager.h"
#include "LMManager.h"
#include "WFSAcceptor.h"
#include "HypothesisLattice.h"

#define BEAM_PRUNING_NUMBER_BINS 25

struct WGToken_;

// structure to keep the lexical unit history
typedef struct _HistoryItem {
   int iLexUnitPron;			// lexical unit (including alternative pronunciations)
   int iEndFrame;				// ending frame (the start frame can  be obtained from the previous lexical unit)
   float fScore;				// global score (accumulated across the whole utterance) includes lm and am scores
   int iPrev;					// previous element 
   int iActive;				// last time the item was active (part of an active token's history) (garbage collection)
   int iWGToken;				// best N paths that arrive at this word history item (each of them has a backpointer)
} HistoryItem;

typedef struct vector<HistoryItem*> VHistoryItem;
typedef struct map<HistoryItem*,bool> MHistoryItem;
typedef struct map<HistoryItem*,LNode*> MHistoryItemLNode;

// word-graph token (word-graph generation)
typedef struct _WGToken {
	int iWordSequence;		// index in the hash table containing unique word sequences
	int iLexUnitPron;			// lexical-unit at the time this token was created (in case it was known)
	float fScore;				// path score
	int iHistoryItem;			// pointer to the previous item in the history
	int iActive;				// last time the token was active (part of an active token's history) (garbage collection)
	int iPrev;					// previous token in the list of available tokens (memory management and garbage collection)	
} WGToken;

typedef list<WGToken*> LWGToken;
typedef vector<WGToken*> VWGToken;
typedef map<WGToken*,bool> MWGToken;

// active state in the search
typedef struct {
	StateX *state;									// state in the acceptor (decoding network)
	float fScore;									// accumulated path score
	HMMStateDecoding *hmmStateDecoding;		// pointer to the HMM-state
	int iHistoryItem;								// lexical unit history
	int iWGToken;									// word-graph token (word-graph generation)
} ActiveState;

typedef list<ActiveState*> LActiveState;

// epsilon active state in the search
typedef struct {
	StateX *state;				// state in the acceptor (decoding network)
	float fScore;				// accumulated path score
	int iHistoryItem;			// lexical unit history
	int iWGToken;				// word-graph token (word-graph generation)
} ActiveStateEpsilon;

typedef struct _HashEntry {
	int iTime;							// time frame of last insertion (-1 initially)
	StateX *state;						// acceptor state
	ActiveState *activeState;		// pointer to the active state
	int iNext;							// next table entry (to handle collisions)
} HashEntry;

typedef struct _WSHashEntry {
	int iTime;					// time frame of last insertion (-1 initially)
	int iLexUnits;				// number of lexical units in the word sequence
	int *iLexUnit;				// lexical units in the word-sequence
	int iNext;					// next table entry (to handle collisions)
} WSHashEntry;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ActiveStateTable {

	public:
	
		unsigned int m_iActiveStatesMax;					// max # of active states = number of buckets in the hash table
		unsigned int m_iHashBuckets;						// # buckets in the hash table
		unsigned int m_iHashEntries;						//	# entries in the hash table (#entries = #buckets + #activeStatesMax)
		HashEntry *m_hashEntries;							// hash table
		int m_iHashEntryCollisionAvailable; 			// available entry in the hash table to store collisions 
		int m_iTimeCurrent;									// current time frame
		ActiveState *m_activeStatesCurrent;				// table of active states
		ActiveState *m_activeStatesNext;					// table of active states
		ActiveStateEpsilon *m_activeStatesEpsilon;		// table of epsilon active states
		ActiveState *m_activeStateAvailable;				// next available entry in the table of active states
		ActiveStateEpsilon *m_activeStateEpsilonHead;	// head of the table of active epsilon states
		ActiveStateEpsilon *m_activeStateEpsilonTail;	// tail of the table of active epsilon states
		unsigned int m_iActiveStatesCurrent;
		int m_iLexUnitPronSilence;
		int m_iLexUnitPronUnknown;
			
		// Note: the table is organized so all the buckets are at the beginning (amounting for 1/2 of the table) and
		// additional entries (1/2 of the table) are placed right after in order to resolve collisions
		
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		HMMStateDecoding *m_hmmStatesDecoding;
		
		// history item management
		unsigned int m_iHistoryItems;					// number of history items allocated
		HistoryItem *m_historyItems;					// history items allocated
		int m_iHistoryItemAvailable;					//	next history item available to be used
		int m_iTimeGarbageCollectionLast;			// last time frame the garbage collection was run
		
		// word-graph token management (lattice generation)
		unsigned int m_iWGTokens;						// number of word-graph tokens allocated
		WGToken *m_wgTokens;								// word-graph tokens allocated
		int m_iWGTokenAvailable;
		unsigned int m_iWGTokensUsed;					// number of word-graph tokens used
		int m_iTimeWGTokenGarbageCollectionLast;	// last time frame the garbage collectionn was run
		int m_iWordSequenceInitial;
		
		// pruning
		float m_fPruningLikelihood;		// beam width
		int m_iPruningMaxStates;			// maximum number of active states per frame
		
		// word-graph generation
		bool m_bLatticeGeneration;				// whether to generate a word-graph	
		int m_iMaxWordSequencesState;			// maximum number of word sequences arriving at any state		
		
		// hast table for hashing word sequences
		unsigned int m_iWSHashBuckets;		// # buckets in the hash table
		unsigned int m_iWSHashEntries;		//	# entries in the hash table (#entries = #buckets + #unique word sequences)
		WSHashEntry *m_wshashEntries;							// hash table
		int m_iWSHashEntryCollisionAvailable;		// next available entry in the hash table to store collisions		
		
		MWGToken m_mWGTokenDeleted;
		map<int,bool> m_mWGTokenUsed;
		
	public:
		
		// debug
		unsigned int m_iStatesExpanded;
		unsigned int m_iStatesActivated;
		unsigned int m_iStatesPruned;

	public:
	
		MHistoryItem m_mHistoryItem;

		// constructor
		ActiveStateTable(float fPruningLikelihood, int iPruningMaxStates, unsigned int iActiveStatesMax, PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMStateDecoding *hmmStatesDecoding, bool bLatticeGeneration, int iMaxWordSequencesState);

		// destructor
		~ActiveStateTable();

		// initialization
		bool initialize();
		
		// return the current active states
		inline ActiveState *getActiveStatesCurrent(unsigned int *iActiveStatesCurrent) {
		
			*iActiveStatesCurrent = m_iActiveStatesCurrent;
		
			return m_activeStatesCurrent;
		}
		
		// return the next active states
		inline ActiveState *getActiveStatesNext(unsigned int *iActiveStatesNext) {
			
			*iActiveStatesNext = m_activeStateAvailable-m_activeStatesNext;
		
			return m_activeStatesNext;
		}
		
		// move to the next time frame 
		inline void nextTimeFrame() {
		
			// increase the current time
			m_iTimeCurrent++;
			
			// swap pointers
			m_iActiveStatesCurrent = m_activeStateAvailable-m_activeStatesNext;
			assert((m_iActiveStatesCurrent > 0) && (m_iActiveStatesCurrent < m_iActiveStatesMax));
			ActiveState *aux = m_activeStatesCurrent;
			m_activeStatesCurrent = m_activeStatesNext;
			m_activeStatesNext = aux;
			
			// reset indices
			m_iHashEntryCollisionAvailable = m_iHashBuckets;
			m_activeStateAvailable = m_activeStatesNext;
			m_activeStateEpsilonHead = m_activeStatesEpsilon;
			m_activeStateEpsilonTail = m_activeStatesEpsilon;	
		}
		
		// reset the time frame
		inline void resetTimeFrame() {
		
			m_iTimeCurrent = 0;
			
			// set all the entries as "old"
			for(int i=0 ; i<m_iHashEntries ; ++i) {
				m_hashEntries[i].iTime = -1;	
			}
				
			m_iHashEntryCollisionAvailable = m_iHashBuckets;
		}
		
		// operations needed for processing a new utterance
		inline void beginUtterance() {
		
			resetTimeFrame();
		
			// reset indices
			m_iHashEntryCollisionAvailable = m_iHashBuckets;
			m_activeStateAvailable = m_activeStatesNext;
			m_activeStateEpsilonHead = m_activeStatesEpsilon;
			m_activeStateEpsilonTail = m_activeStatesEpsilon;
		
			// reset history items
			for(unsigned int i=0 ; i < m_iHistoryItems-1 ; ++i) {
				m_historyItems[i].iActive = -1;
				m_historyItems[i].iPrev = i+1;
			}
			m_historyItems[m_iHistoryItems-1].iActive = -1;
			m_historyItems[m_iHistoryItems-1].iPrev = -1;
			m_iHistoryItemAvailable = 0;
			m_iTimeGarbageCollectionLast = -1;
			
			// lattice generation
			if (m_bLatticeGeneration == true) {
				
				// reset word-graph tokens
				for(unsigned int i=0 ; i < m_iWGTokens-1 ; i += m_iMaxWordSequencesState) {
					m_wgTokens[i].iActive = -1;
					m_wgTokens[i].iPrev = i+m_iMaxWordSequencesState;
				}
				m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iActive = -1;
				m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iPrev = -1;
				m_iWGTokenAvailable = 0;	
				
				// set all the entries as "old" (hash table of word-sequences)
				cleanHashWordSequences();
				//for(int i=0 ; i < m_iWSHashEntries ; ++i) {
				//	m_wshashEntries[i].iTime = -1;	
				//}	
				m_iWSHashEntryCollisionAvailable = m_iWSHashBuckets;
			}
		
			return;
		}
		
		// operations needed to end the processing of an utterance
		inline void endUtterance() {
		
			m_iActiveStatesCurrent = 0;	
		
			return;
		}
		
		// activates the initial state
		void activateStateInitial(StateX *state);
		
		// activates a state if not active, otherwise updates score
		inline void activateState(StateX *state, float fScore, float *fScoreBest, HMMStateDecoding *hmmStateDecoding, 
			int iHistoryItem, int iWGToken, float fScoreAdded) {
		
			assert(fScoreAdded == 0.0);
		
			unsigned int iEntry = (((unsigned long)state)/sizeof(StateX))%m_iHashBuckets;
			HashEntry &entry = m_hashEntries[iEntry];
			// old entry
			if (entry.iTime < m_iTimeCurrent) {
				entry.iTime = m_iTimeCurrent;
				entry.state = state;
				entry.activeState = m_activeStateAvailable;
				if (m_bLatticeGeneration == true) {
					entry.activeState->iWGToken = iWGToken;	
				}
				entry.iNext = -1;
				// overwrite values
				m_activeStateAvailable->state = state;
				m_activeStateAvailable->fScore = fScore;
				if (*fScoreBest < fScore) {
					*fScoreBest = fScore;
				}
				m_activeStateAvailable->hmmStateDecoding = hmmStateDecoding;
				m_activeStateAvailable->iHistoryItem = iHistoryItem;
				++m_activeStateAvailable;
				assert(m_activeStateAvailable-m_activeStatesNext < m_iActiveStatesMax);
				++m_iStatesActivated;
				return;
			} 
			// hit
			else if (entry.activeState->state == state) {
				assert(entry.iTime == m_iTimeCurrent);
				assert(entry.activeState->hmmStateDecoding == hmmStateDecoding);
				if (m_bLatticeGeneration == false) {	
					// update the score if necessary (the history is updated too)
					if (entry.activeState->fScore < fScore) {	
						entry.activeState->fScore = fScore;
						entry.activeState->iHistoryItem = iHistoryItem;
						if (*fScoreBest < fScore) {
							*fScoreBest = fScore;	
						}
					}
				}	
				// merge word sequences (lattice generation)
				else {
					mergeWordSequences(entry.activeState->iWGToken,iWGToken);
					// update the score based on the result from the merging process
					entry.activeState->fScore = (entry.activeState->iWGToken+m_wgTokens)[0].fScore;
					entry.activeState->iHistoryItem = (entry.activeState->iWGToken+m_wgTokens)[0].iHistoryItem;
					if (entry.activeState->fScore < fScore) {	
						if (*fScoreBest < fScore) {
							*fScoreBest = fScore;	
						}
					}
					deleteWGToken(iWGToken);
				}
				return;
			} 
			// collision
			else {
				//printf("collision\n");
				int *iHashEntryAux = &entry.iNext;
				while(*iHashEntryAux != -1) {
					HashEntry *hashEntryAux = m_hashEntries+*iHashEntryAux;
					assert(hashEntryAux->iTime == m_iTimeCurrent);
					// hit
					if (hashEntryAux->state == state) {
						assert(hashEntryAux->activeState->hmmStateDecoding == hmmStateDecoding);
						// update the score if necessary (the history is updated too)
						if (hashEntryAux->activeState->fScore < fScore) {
							hashEntryAux->activeState->fScore = fScore;
							hashEntryAux->activeState->iHistoryItem = iHistoryItem;
							if (*fScoreBest < fScore) {
								*fScoreBest = fScore;
							}
						}
						// merge word sequences (word-graph generation)
						if (m_bLatticeGeneration == true) {
							mergeWordSequences(hashEntryAux->activeState->iWGToken,iWGToken);
							hashEntryAux->activeState->fScore = (hashEntryAux->activeState->iWGToken+m_wgTokens)[0].fScore;
							hashEntryAux->activeState->iHistoryItem = (hashEntryAux->activeState->iWGToken+m_wgTokens)[0].iHistoryItem;
							deleteWGToken(iWGToken);
						}
						return;
					}
					iHashEntryAux = &hashEntryAux->iNext;
				}
				assert(*iHashEntryAux == -1);
				// the state was not found, create it as a collision
				*iHashEntryAux = m_iHashEntryCollisionAvailable;
				assert(m_iHashEntryCollisionAvailable != m_iHashEntries);
				(m_hashEntries+*iHashEntryAux)->iTime = m_iTimeCurrent;
				(m_hashEntries+*iHashEntryAux)->state = state;
				(m_hashEntries+*iHashEntryAux)->activeState = m_activeStateAvailable;
				(m_hashEntries+*iHashEntryAux)->iNext = -1;
				++m_iHashEntryCollisionAvailable;
				m_activeStateAvailable->state = state;
				m_activeStateAvailable->fScore = fScore;
				m_activeStateAvailable->iHistoryItem = iHistoryItem;
				if (*fScoreBest < fScore) {
					*fScoreBest = fScore;
				}
				m_activeStateAvailable->hmmStateDecoding = hmmStateDecoding;
				++m_activeStateAvailable;
				assert(m_activeStateAvailable-m_activeStatesNext < m_iActiveStatesMax);
				++m_iStatesActivated;
				if (m_bLatticeGeneration) {
					(m_hashEntries+*iHashEntryAux)->activeState->iWGToken = iWGToken;
				}
				return;
			}
		}

		// activates a state if not active, otherwise updates the score
		// the overall best score is not updated from here
		inline void activateStateEpsilon(StateX *state, float fScore, int iHistoryItem, int iWGToken, float fScoreAdded) {
			
			assert(fScoreAdded == 0.0);	
			
			//printf("activating epsilon (head: %x)\n",m_activeStateEpsilonHead);
		
			unsigned int iEntry = (((unsigned long)state)/sizeof(StateX))%m_iHashBuckets;
			HashEntry &entry = m_hashEntries[iEntry];
			// old entry
			if (entry.iTime < m_iTimeCurrent) {
				entry.iTime = m_iTimeCurrent;
				entry.state = state;
				entry.activeState = (ActiveState*)m_activeStateEpsilonTail;
				entry.iNext = -1;
				//printf("epsilon activated: %x\n",m_activeStateEpsilonTail);
				// overwrite values
				m_activeStateEpsilonTail->state = state;
				m_activeStateEpsilonTail->fScore = fScore;
				m_activeStateEpsilonTail->iHistoryItem = iHistoryItem;
				if (m_bLatticeGeneration == true) {
					m_activeStateEpsilonTail->iWGToken = iWGToken;	
				}
				// treat it like a circular buffer
				++m_activeStateEpsilonTail;
				if (m_activeStateEpsilonTail == m_activeStatesEpsilon+m_iActiveStatesMax) {
					m_activeStateEpsilonTail = m_activeStatesEpsilon;
				}
				assert(m_activeStateEpsilonTail != m_activeStateEpsilonHead);

				return;
			} 
			// hit
			else if (entry.activeState->state == state) {
				assert(entry.iTime == m_iTimeCurrent);
				ActiveStateEpsilon *activeStateEpsilon = ((ActiveStateEpsilon*)entry.activeState);
				if (m_bLatticeGeneration == false) {
					// update the score if necessary (the history is updated too)
					if (activeStateEpsilon->fScore < fScore) {	
						activeStateEpsilon->fScore = fScore;
						activeStateEpsilon->iHistoryItem = iHistoryItem;
					}
				} 
				// lattice generation
				else {
					// merge word sequences (word-graph generation)
					mergeWordSequences(activeStateEpsilon->iWGToken,iWGToken);
					activeStateEpsilon->fScore = (activeStateEpsilon->iWGToken+m_wgTokens)[0].fScore;
					activeStateEpsilon->iHistoryItem = (activeStateEpsilon->iWGToken+m_wgTokens)[0].iHistoryItem;
					deleteWGToken(iWGToken);
				}
				return;
			} 
			// collision
			else {
				int *iHashEntryAux = &entry.iNext;
				while(*iHashEntryAux != -1) {
					HashEntry *hashEntryAux = m_hashEntries+*iHashEntryAux;
					assert(hashEntryAux->iTime == m_iTimeCurrent);
					// hit
					if (hashEntryAux->state == state) {
						ActiveStateEpsilon *activeStateEpsilon = ((ActiveStateEpsilon*)hashEntryAux->activeState);
						// update the score if necessary (the history is updated too)
						if (activeStateEpsilon->fScore < fScore) {
							activeStateEpsilon->fScore = fScore;
							activeStateEpsilon->iHistoryItem = iHistoryItem;
						}
						// merge word sequences (word-graph generation)
						if (m_bLatticeGeneration == true) {
							mergeWordSequences(activeStateEpsilon->iWGToken,iWGToken);
							activeStateEpsilon->fScore = (activeStateEpsilon->iWGToken+m_wgTokens)[0].fScore;
							activeStateEpsilon->iHistoryItem = (activeStateEpsilon->iWGToken+m_wgTokens)[0].iHistoryItem;
							deleteWGToken(iWGToken);
						}
						return;
					}
					iHashEntryAux = &hashEntryAux->iNext;
				}
				assert(*iHashEntryAux == -1);
				// the state was not found, create it as a collision
				*iHashEntryAux = m_iHashEntryCollisionAvailable;
				assert(m_iHashEntryCollisionAvailable != m_iHashEntries);
				(m_hashEntries+*iHashEntryAux)->iTime = m_iTimeCurrent;
				(m_hashEntries+*iHashEntryAux)->state = state;
				(m_hashEntries+*iHashEntryAux)->activeState = ((ActiveState*)m_activeStateEpsilonTail);
				(m_hashEntries+*iHashEntryAux)->iNext = -1;
				++m_iHashEntryCollisionAvailable;
				// active state values
				m_activeStateEpsilonTail->state = state;
				m_activeStateEpsilonTail->fScore = fScore;
				m_activeStateEpsilonTail->iHistoryItem = iHistoryItem;
				// treat it like a circular buffer
				++m_activeStateEpsilonTail;
				if (m_activeStateEpsilonTail == m_activeStatesEpsilon+m_iActiveStatesMax) {
					m_activeStateEpsilonTail = m_activeStatesEpsilon;
				}
				assert(m_activeStateEpsilonTail != m_activeStateEpsilonHead);
				
				if (m_bLatticeGeneration) {
					((ActiveStateEpsilon*)(m_hashEntries+*iHashEntryAux)->activeState)->iWGToken = iWGToken;	
				}
				return;
			}
		}
		
		// apply beam pruning to the set of active states
		void beamPruning(float *fScoreBest);
		
		// compute the load factor of the hash table containing active states
		float computeLoadFactor();
		
		// compute the collision factor (# collisions / # valid entries in the table)
		float computeCollisionFactor();
		
		// compute the load factor of the hash table containing unque word sequences
		float computeLoadFactorHashWordSequences(int *iBucketsUsed, int *iCollisions);
		
		// print stats connected to the hash table containing unique word-sequences
		void printHashTableWordSequences();	
		
		// shows object information
		void printInfo();	
		
		// process epsilon transitions in topological order
		void processEpsilonTransitions(float *fFeatureVector, float *fScoreBest);	
		
		// garbage collection of history items
		// (1) it starts by marking the active items by traversing back items from the active states
		// (2) it adds inactive items to the queue of available items
		void historyItemGarbageCollection();
		
		// mark a history item as active
		int markActive(int iHistoryItem);	
					
		// return an unused history item
		inline int newHistoryItem() {
		
			if (m_iHistoryItemAvailable == -1) {
				historyItemGarbageCollection();
			}
			assert(m_iHistoryItemAvailable != -1);
		
			int iReturn = m_iHistoryItemAvailable;
			m_iHistoryItemAvailable = m_historyItems[m_iHistoryItemAvailable].iPrev; 
		
			return iReturn;
		}	
		
		// return a new word-graph token structure
		inline int newWGToken() {
		
			if (m_iWGTokenAvailable == -1) {
				historyItemGarbageCollection();
			} 
			assert(m_iWGTokenAvailable != -1);
			
			int iWGTokenAux = m_iWGTokenAvailable;
			m_iWGTokenAvailable = (m_iWGTokenAvailable+m_wgTokens)->iPrev;
			
			(m_wgTokens+iWGTokenAux)->iActive = m_iTimeCurrent;
			
			return iWGTokenAux;
		}	
		
		// return a new word-graph token structure
		inline int newWGToken(int iWordSequence, float fScore, int iHistoryItem) {
		
			int iWGTokenAux = newWGToken();
		
			(iWGTokenAux+m_wgTokens)[0].iWordSequence = iWordSequence;
			(iWGTokenAux+m_wgTokens)[0].iLexUnitPron = m_iLexUnitPronUnknown;
			(iWGTokenAux+m_wgTokens)[0].fScore = fScore;
			(iWGTokenAux+m_wgTokens)[0].iHistoryItem = iHistoryItem;
			(iWGTokenAux+m_wgTokens)[1].iWordSequence = -1;
		
			return iWGTokenAux;
		}
		
		// return a new word-graph token structure that is a copy of an existing one
		inline int newWGToken(int iWGToken) {
		
			int iWGTokenAux = newWGToken();
			
			WGToken *wgToken = (iWGToken+m_wgTokens);
			WGToken *wgTokenAux = (iWGTokenAux+m_wgTokens);
			
			for(int i=0 ; i < m_iMaxWordSequencesState ; ++i) {
				if (wgToken[i].iWordSequence == -1) {
					wgTokenAux[i].iWordSequence = -1;	
					break;
				}
				wgTokenAux[i].iWordSequence = wgToken[i].iWordSequence;
				wgTokenAux[i].iLexUnitPron = wgToken[i].iLexUnitPron;
				wgTokenAux[i].fScore = wgToken[i].fScore;
				wgTokenAux[i].iHistoryItem = wgToken[i].iHistoryItem;
			}
			
			return iWGTokenAux;
		}
		
		// delete a word-graph token structure
		inline void deleteWGToken(int iWGToken) {
			
			assert(iWGToken >= 0);
			(iWGToken+m_wgTokens)->iPrev = m_iWGTokenAvailable;
			m_iWGTokenAvailable = iWGToken;
		}
		
		// recovers the best path from the list of active states
		BestPath *getBestPath(int iFeatureVectors);	
		
		// print the given best paths (useful for debugging)
		void printBestPaths(unsigned int iPaths);	
		
		// comparision function used to sort active states by their score
		inline static bool activeStateComparisonScore(const ActiveState *activeState1, const ActiveState *activeState2) {
			
			return (activeState1->fScore > activeState2->fScore);
		}
		
		// print a transition
		void printTransition(TransitionX *transitionX) {
		
			// epsilon:
			if (transitionX->iSymbol == EPSILON_TRANSITION) {
				printf("eps %.2f %x\n",transitionX->fWeight,transitionX->state);
			}
			// fake transition
			else if (transitionX->iSymbol & FAKE_TRANSITION) {
				printf("eps %.2f %x\n",transitionX->fWeight,transitionX->state);	
			} 
			// lexical unit
			else if (transitionX->iSymbol & LEX_UNIT_TRANSITION) {
				int iLexUnit = transitionX->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				printf("%-20s %.2f %x\n",m_lexiconManager->getStrLexUnit(iLexUnit),transitionX->fWeight,transitionX->state);
			}
			// HMM-state
			else {
				assert(transitionX->iSymbol < LEX_UNIT_TRANSITION);
				int iHMMState = transitionX->iSymbol;
				printf("hmm-state: %10u %.2f %x\n",transitionX->iSymbol,transitionX->fWeight,transitionX->state);
			}	
		}
		
		// print a subgraph
		void printSubgraph(StateX *stateX) {
		
			MStateX mStateX;
			mStateX.insert(MStateX::value_type(stateX,true));
			VStateX vStateX;
			vStateX.push_back(stateX);
			
			while(vStateX.empty() == false) {
		
				StateX *stateXAux = vStateX.back();
				vStateX.pop_back();	
				assert(stateXAux != NULL);
				printf("-> state: %x\n",stateXAux);
				
				TransitionX *transition = *stateXAux;
				TransitionX *transitionEnd = *(stateXAux+1);
				while(transition != transitionEnd) {
					if (!(transition->iSymbol & FAKE_TRANSITION)) {
						printTransition(transition);	
						if (!(transition->iSymbol & LEX_UNIT_TRANSITION)) {
							if (mStateX.find(transition->state) == mStateX.end()) {
								mStateX.insert(MStateX::value_type(transition->state,true));
								vStateX.push_back(transition->state);
							}
						}						
					}
					transition++;
				}
			}
			
			return;
		}
		
		// hash a word sequence and returns a pointer to the bucket
		inline int hashWordSequence(HistoryItem *historyItem) {
		
			assert(historyItem != NULL);
		
			// apply the hash function and get the bucket (xor is commutative)
			int i = 0;
			int iAux;
			int iEntry = 0;
			HistoryItem *historyItemAux = historyItem;
			do {
				LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(historyItemAux->iLexUnitPron);
				iAux = lexUnit->iLexUnit;
				// ignore silence and filler lexical units
				if (lexUnit->iType == LEX_UNIT_TYPE_STANDARD) {
					iAux <<= 8*i;
					iEntry ^= iAux;
					++i;
					i %= sizeof(int);
				}
				if (historyItemAux->iPrev == -1) {
					break;
				}
				historyItemAux = m_historyItems+historyItemAux->iPrev;
			} while(1);
			iEntry %= m_iWSHashBuckets;	
			WSHashEntry &entry = m_wshashEntries[iEntry];
			
			// old entry
			if (entry.iTime == -1) {
				entry.iTime = m_iTimeCurrent;
				entry.iNext = -1;
				storeWordSequence(&entry,historyItem);
				return iEntry;
			} 
			// hit
			else if (compareWordSequencesStandard(&entry,historyItem)) {
				return iEntry;
			} 
			// collision
			else {
				// collision: try to find a hit and otherwise append a new entry to the bucket
				int *iEntryAux = &entry.iNext;
				while(*iEntryAux != -1) {
					// hit
					if (compareWordSequencesStandard(m_wshashEntries+(*iEntryAux),historyItem)) {	
						assert(*iEntryAux >= 0);
						assert(*iEntryAux < m_iWSHashEntries);
						return *iEntryAux;
					}
					iEntryAux = &(m_wshashEntries+(*iEntryAux))->iNext;
				}
				
				// the word sequence was not found, create it as a collision
				*iEntryAux = m_iWSHashEntryCollisionAvailable;
				assert(*iEntryAux != m_iWSHashEntries);
				(m_wshashEntries+(*iEntryAux))->iTime = m_iTimeCurrent;
				storeWordSequence(m_wshashEntries+(*iEntryAux),historyItem); 
				(m_wshashEntries+(*iEntryAux))->iNext = -1;
				++m_iWSHashEntryCollisionAvailable;
				return *iEntryAux;	
			}
		}
		
		// clean the hash-table
		void cleanHashWordSequences() {		
			
			// set all the entries as "old" and free memory
			for(int i=0 ; i < m_iWSHashEntries ; ++i) {
				if (m_wshashEntries[i].iTime != -1) {
					delete [] m_wshashEntries[i].iLexUnit;
					m_wshashEntries[i].iTime = -1;
				}	
			}	
		}
				
		// copy a word sequence from the history-element ot 
		inline void storeWordSequence(WSHashEntry *entry, HistoryItem *historyItem) {
		
			// (a) count lexical-units
			entry->iLexUnits = 0;
			int iHistoryItemAux = historyItem-m_historyItems;
			do {
				if (isStandard((iHistoryItemAux+m_historyItems)->iLexUnitPron)) {
					entry->iLexUnits++;
				}
				iHistoryItemAux = (iHistoryItemAux+m_historyItems)->iPrev;
			} while(iHistoryItemAux != -1);
			
			// (b) copy lexical units
			entry->iLexUnit = new int[entry->iLexUnits];
			int iIndex = 0;
			iHistoryItemAux = historyItem-m_historyItems;
			do {
				if (isStandard((iHistoryItemAux+m_historyItems)->iLexUnitPron)) {
					entry->iLexUnit[iIndex++] = m_lexiconManager->getLexUnitPron((iHistoryItemAux+m_historyItems)->iLexUnitPron)->iLexUnit;
				}
				iHistoryItemAux = (iHistoryItemAux+m_historyItems)->iPrev;
			} while(iHistoryItemAux != -1);
			assert(iIndex == entry->iLexUnits);	
		
			return;
		}
		
		// return whether the lexical unit is a silence or filler
		inline bool isStandard(int iLexUnitPron) {
		
			return m_lexiconManager->isStandard(m_lexiconManager->getLexUnitPron(iLexUnitPron));
		}
		
		// return whether two word sequences are identical 
		// - ignores silence/filler symbols
		// - ignores alternative pronunciations
		inline bool compareWordSequencesStandard(WSHashEntry *entry, HistoryItem *historyItem) {
			
			int iHistoryItem = historyItem-m_historyItems;
			int iIndex = 0;
			
			do {
			
				// get the next "real" lexical unit
				while(iHistoryItem != -1) {
					if (isStandard((iHistoryItem+m_historyItems)->iLexUnitPron)) {
						break;
					}
					iHistoryItem = (iHistoryItem+m_historyItems)->iPrev;
				}
				
				if (iHistoryItem == -1) {
					return (iIndex == entry->iLexUnits);	
				} else if (iIndex == entry->iLexUnits) {
					return false;
				} else if (entry->iLexUnit[iIndex++] != 
					m_lexiconManager->getLexUnitNoPron((iHistoryItem+m_historyItems)->iLexUnitPron)) {
					return false;
				}
				
				iHistoryItem = (iHistoryItem+m_historyItems)->iPrev;
				
			} while(1);
		}	

		
		// return whether the lexical unit is a silence or filler
		inline bool isSilenceFiller(int iLexUnitPron) {
		
			return m_lexiconManager->isFiller(m_lexiconManager->getLexUnitPron(iLexUnitPron));
		}
		
		// merge two word sequences by keeping the N best in wgToken2
		bool mergeWordSequences(int iWGToken1, int iWGToken2);
				
		inline void swap(WGToken *wgToken1, WGToken *wgToken2) {
			
			WGToken wgTokenAux;
			memcpy(&wgTokenAux,wgToken1,sizeof(WGToken));
			memcpy(wgToken1,wgToken2,sizeof(WGToken));
			memcpy(wgToken2,&wgTokenAux,sizeof(WGToken));
			
			return;
		}
		
		// comparison function
		inline static bool wgTokenEqual(WGToken *wgToken1, WGToken *wgToken2) {
			
			return (wgToken1->iWordSequence == wgToken2->iWordSequence);
		}
		
		// comparison function
		inline static bool wgTokenBetterPath(WGToken *wgToken1, WGToken *wgToken2) {
			
			return (wgToken1->fScore > wgToken2->fScore);
		}	
		
		// build a hyothesis lattice for the utterance
		HypothesisLattice *getHypothesisLattice();
		
		// build a word graph from 
		HypothesisLattice *generateLattice(HistoryItem *historyItem);
		
		// print a history item
		inline void printHistoryItem(HistoryItem *historyItem) {
		
			printf("%5d %20s %10.2f %x\n",historyItem->iEndFrame,m_lexiconManager->getStrLexUnitPron(historyItem->iLexUnitPron),
				historyItem->fScore,historyItem);
		
			return;
		}
		
		// print a token set
		inline void printWGToken(int iWGToken) {
		
			printWGToken(iWGToken+m_wgTokens);
		}
		
		// print a token set
		inline void printWGToken(WGToken *wgToken) {
		
			printf("------------------------------------------------------\n");
			for(int i=0 ; i < m_iMaxWordSequencesState ; ++i) {
				if (wgToken[i].iWordSequence == -1) {
					break;
				}
				printf("%10d %30.20f %d\n",wgToken[i].iWordSequence,wgToken[i].fScore,wgToken[i].iHistoryItem);
			}
			printf("------------------------------------------------------\n");
		
			return;
		}
		
		void checkTokens();
	
};

#endif
