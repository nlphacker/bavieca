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


#ifndef DYNAMICDECODERX_H
#define DYNAMICDECODERX_H

#include "DynamicNetworkX.h"
#include "HypothesisLattice.h"
#include "LexiconManager.h"

namespace Bavieca {

class BestPath;
class HMMManager;
class PhoneSet;
class LMLookAhead;
class LMFSM;
class LMManager;

#define NUMBER_BINS_HISTOGRAM							50
#define NUMBER_BINS_HISTOGRAM_WITHIN_NODE			100

// language model transition (auxiliar structure)
typedef struct {
	int iLMState;				// next lm-state
	float fScoreLM;			// associated score
} LMTransition;

// entry in the hash table that keeps unique word sequences (lattice generation)
typedef struct _WSHashEntry {
	int iTime;					// time frame of last insertion (-1 initially)
	int iLexUnits;				// number of lexical units in the word sequence
	int *iLexUnit;				// lexical units in the word-sequence
	int iNext;					// next table entry (to handle collisions)
} WSHashEntry;

// structure to keep the lexical unit history
typedef struct _HistoryItem {
   int iLexUnitPron;			// lexical unit (including alternative pronunciations)
   int iEndFrame;				// ending frame (the start frame can  be obtained from the previous lexical unit)
   float fScore;				// global score (accumulated across the whole utterance) includes lm and am scores
   int iPrev;					// previous history item
   int iActive;				// last time the item was active (part of an active token's history) (garbage collection)
   int iWGToken;				// best N paths that arrive at this word history item (each of them has a backpointer)
} HistoryItem;

typedef vector<HistoryItem*> VHistoryItem;
typedef list<HistoryItem*> LHistoryItem;
typedef map<HistoryItem*,bool> MHistoryItem;
typedef map<HistoryItem*,LNode*> MHistoryItemLNode;

// word-graph token (word-graph generation)
typedef struct _WGToken {
	int iWordSequence;				// index in the hash table containing unique word sequences
	int iLexUnitPron;					// lexical-unit at the time this token was created (in case it was known)
	//add support for this field and remove added code used for debugging
	float fScore;						// path score
	int iHistoryItem;					// index of the previous item in the history
	int iActive;						// last time the token was active (part of an active token's history) (garbage collection)
	int iPrev;							// previous token in the list of available tokens (memory management and garbage collection)	
} WGToken;

typedef list<WGToken*> LWGToken;
typedef vector<WGToken*> VWGToken;
typedef map<WGToken*,bool> MWGToken;

// token
typedef struct {
	float fScore;						// path score
	HMMStateDecoding *state;		// hmm-state
	int iLMState;						// language model state
	int iLexUnitPron;					// lexical unit (including alternative pronunciations)
	int iNode;							// index of the node where the token is
	int iHistoryItem;					// index of the history item
	int iWGToken;
	// look-ahead
	int iLANode;						// node in the look ahead tree
	float *fLAScores;					// language model look-ahead scores	
} Token;

// active token
typedef struct {
	int iLMState;						// language model state
	int iToken;							// active token
} ActiveToken;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class DynamicDecoderX {

	private:
	
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManager;
		LexiconManager *m_lexiconManager;
		LMFSM *m_lmFSM;
		LMManager *m_lmManager;
		unsigned char m_iNGram;
		DynamicNetworkX *m_dynamicNetwork;
		
		// network properties
		int m_iArcs;
		DArc *m_arcs;
		int m_iNodes;
		DNode *m_nodes;
		
		// pruning parameters
		int m_iMaxActiveNodes;				// maximum number of active arcs
		int m_iMaxActiveNodesWE;			// maximum number of active arcs at word-ends
		int m_iMaxActiveTokensNode;		// maximum number of active tokens within an arc
		float m_fBeamWidthNodes;			// beam width for all arcs
		float m_fBeamWidthNodesWE;			// beam width for all arcs at word-ends
		float m_fBeamWidthTokensNode;		// beam width for all tokens within an arc
		
		// scaling factor
		float m_fLMScalingFactor;
		
		// current time frame
		int m_iTimeCurrent;
		
		bool m_bInitialized;
		
		// tokens
		int m_iTokensMax;						// number of tokens allocated
		Token *m_tokensCurrent;
		Token *m_tokensNext;
		int m_iTokensNext;
		
		// active nodes
		DNode **m_nodesActiveCurrent;
		DNode **m_nodesActiveNext;
		int m_iNodesActiveCurrent;
		int m_iNodesActiveNext;
		int m_iNodesActiveCurrentMax;
		int m_iNodesActiveNextMax;
		
		// active tokens
		int m_iTokensNodeMax;
		ActiveToken *m_activeTokenCurrent;
		ActiveToken *m_activeTokenNext;
		int m_iActiveTokenTables;
		int m_iActiveTokenMax;
		
		float m_fScoreBest;
		float m_fScoreBestWE;
		
		// history item management
		unsigned int m_iHistoryItems;					// number of history items allocated
		HistoryItem *m_historyItems;					// history items allocated
		int m_iHistoryItemBegSentence;				// initial history item
		int m_iHistoryItemAvailable;					//	next history item available to be used
		int m_iTimeGarbageCollectionLast;			// last time frame the garbage collection was run
		// auxiliar arrays (used for token expansion)
		int *m_iHistoryItemsAuxBuffer;
		int *m_iHistoryItemsAux;
		int m_iHistoryItemsAuxSize;
		
		// word-graph generation
		bool m_bLatticeGeneration;				// whether to generate a word-graph	
		int m_iMaxWordSequencesState;			// maximum number of word sequences arriving at any state	
		
		// word-graph token management (lattice generation)
		unsigned int m_iWGTokens;						// number of word-graph tokens allocated
		WGToken *m_wgTokens;								// word-graph tokens allocated
		int m_iWGTokenAvailable;
		int *m_iWordSequenceAux;
				
		// hash table for hashing word sequences
		unsigned int m_iWSHashBuckets;		// # buckets in the hash table
		unsigned int m_iWSHashEntries;		//	# entries in the hash table (#entries = #buckets + #unique word sequences)
		WSHashEntry *m_wshashEntries;							// hash table
		int m_iWSHashEntryCollisionAvailable;		// next available entry in the hash table to store collisions
		
		// utterance information
		int m_iFeatureVectorsUtterance;	
		
		// unknown lexical unit
		int m_iLexUnitPronUnknown;
		
		// language model look-ahead
		LMLookAhead *m_lmLookAhead;
		
		// create a new token
		inline int newToken() {	
		
			assert(m_iTokensNext < m_iTokensMax);
			int iToken = m_iTokensNext;	
			++m_iTokensNext;
		
			return iToken;
		}
		
		inline int newActiveTokenTable() {
		
			assert(m_iActiveTokenTables+m_iTokensNodeMax < m_iActiveTokenMax);
			m_iActiveTokenTables += m_iTokensNodeMax; 
		
			return m_iActiveTokenTables-m_iTokensNodeMax;
		}
		
		// return an unused history item
		inline int newHistoryItem() {
		
			if (m_iHistoryItemAvailable == -1) {
				if (m_bLatticeGeneration == false) {
					historyItemGarbageCollection();
				} else {
					historyItemGarbageCollectionLattice();
				}
			}
			assert(m_iHistoryItemAvailable != -1);
			
			int iReturn = m_iHistoryItemAvailable;
			m_historyItems[m_iHistoryItemAvailable].iWGToken = -1;
			m_iHistoryItemAvailable = m_historyItems[m_iHistoryItemAvailable].iPrev; 
		
			return iReturn;
		}
		
		// return whether a history item is inactive (debugging)
		inline bool inactive(int iHistoryItem) {
		
			int iAux = m_iHistoryItemAvailable;
			while(iAux != -1) {
				if (iAux == iHistoryItem) {
					return true;
				}
				iAux = m_historyItems[iAux].iPrev;
			}
			
			return false;
		}
		
		// swap token tables (after processing a feature vector)
		inline void swapTokenTables() {
				
			// swap token tables
			Token *tokenAux = m_tokensCurrent;
			m_tokensCurrent = m_tokensNext;
			m_tokensNext = tokenAux;
			m_iTokensNext = 0;
			
			// swap active-token tables
			ActiveToken *activeTokenAux = m_activeTokenNext;
			m_activeTokenNext = m_activeTokenCurrent;
			m_activeTokenCurrent = activeTokenAux;
			m_iActiveTokenTables = 0;
		}		
		
		// root-node expansion
		void expandRoot(VectorBase<float> &vFeatureVector);
		
		// regular expansion
		void expand(VectorBase<float> &vFeatureVector, int t);
		
		// expand a series of tokens to a hmm-state
		void expandToHMM(DNode *node, DArc *arcNext, VectorBase<float> &vFeatureVector, int t);
		
		// expand a series of tokens to a hmm-state after obsering a new word
		void expandToHMMNewWord(DNode *node, DArc *arcNext, LexUnit *lexUnit, LMTransition *lmTransition, 
			VectorBase<float> &vFeatureVector, int t);	
		
		// pruning (token based)
		void pruningOriginal();
		
		// pruning (token based)
		void pruning();		
		
		// marks unused history items as available
		void historyItemGarbageCollection();
				
		// marks unused history items as available (lattice generation)
		void historyItemGarbageCollectionLattice();
		
		// get the destination lexical units
		void getDestinationMonophoneLexUnits(DNode *node, VLexUnit &vLexUnitDest);
		
		// prune active tokens that are not in the top-N within a node 
		void pruneExtraTokens(DNode *node);
				
		// word-graph generation ---------------------------------------------------------------------
		
		// return a new word-graph token structure
		inline int newWGToken() {
		
			if (m_iWGTokenAvailable == -1) {
				historyItemGarbageCollectionLattice();
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
		
		// return a new word-graph token structure that is a copy of an existing one but 
		// with updated scores
		inline int newWGToken(int iWGToken, float fScore) {
		
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
				wgTokenAux[i].fScore = wgToken[i].fScore+fScore;
				wgTokenAux[i].iHistoryItem = wgToken[i].iHistoryItem;
			}
			
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
		
		// attach a lexical unit to all the wg tokens
		inline void attachLexUnit(int iWGToken, LexUnit *lexUnit) {
		
			WGToken *wgToken = m_wgTokens+iWGToken;
			for(int i=0 ; (i < m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1) ; ++i) {
				wgToken[i].iLexUnitPron = lexUnit->iLexUnitPron;
			}
		}
		
		
		// delete a word-graph token structure
		inline void deleteWGToken(int iWGToken) {
		
			assert(iWGToken >= 0);
			(iWGToken+m_wgTokens)->iPrev = m_iWGTokenAvailable;
			m_iWGTokenAvailable = iWGToken;
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
						assert(*iEntryAux < (int)m_iWSHashEntries);
						return *iEntryAux;
					}
					iEntryAux = &(m_wshashEntries+(*iEntryAux))->iNext;
				}
				
				// the word sequence was not found, create it as a collision
				*iEntryAux = m_iWSHashEntryCollisionAvailable;
				assert(*iEntryAux != (int)m_iWSHashEntries);
				(m_wshashEntries+(*iEntryAux))->iTime = m_iTimeCurrent;
				storeWordSequence(m_wshashEntries+(*iEntryAux),historyItem); 
				(m_wshashEntries+(*iEntryAux))->iNext = -1;
				++m_iWSHashEntryCollisionAvailable;
				return *iEntryAux;
			}
		}
		
		// copy a word sequence from the history-element ot 
		inline void storeWordSequence(WSHashEntry *entry, HistoryItem *historyItem) {
		
			// (a) count lexical-units
			entry->iLexUnits = 0;
			int iHistoryItemAux = (int)(historyItem-m_historyItems);
			do {
				if (isStandard((iHistoryItemAux+m_historyItems)->iLexUnitPron)) {
					entry->iLexUnits++;
				}
				iHistoryItemAux = (iHistoryItemAux+m_historyItems)->iPrev;
			} while(iHistoryItemAux != -1);
			
			// (b) copy lexical units
			entry->iLexUnit = new int[entry->iLexUnits];
			int iIndex = 0;
			iHistoryItemAux = (int)(historyItem-m_historyItems);
			do {
				if (isStandard((iHistoryItemAux+m_historyItems)->iLexUnitPron)) {
					entry->iLexUnit[iIndex++] = m_lexiconManager->getLexUnitPron((iHistoryItemAux+m_historyItems)->iLexUnitPron)->iLexUnit;
				}
				iHistoryItemAux = (iHistoryItemAux+m_historyItems)->iPrev;
			} while(iHistoryItemAux != -1);
			assert(iIndex == entry->iLexUnits);	
		}
		
		// return whether the lexical unit is a silence or filler
		inline bool isStandard(int iLexUnitPron) {
		
			return m_lexiconManager->isStandard(m_lexiconManager->getLexUnitPron(iLexUnitPron));
		}
		
		// return whether two word sequences are identical 
		// - ignores silence/filler symbols
		// - ignores alternative pronunciations
		inline bool compareWordSequencesStandard(WSHashEntry *entry, HistoryItem *historyItem) {
			
			int iHistoryItem = (int)(historyItem-m_historyItems);
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
		
		// merge two word sequences by keeping the N best in wgToken1
		bool mergeWordSequences(int iWGToken1, int iWGToken2);
				
		inline void swap(WGToken *wgToken1, WGToken *wgToken2) {
			
			WGToken wgTokenAux;
			memcpy(&wgTokenAux,wgToken1,sizeof(WGToken));
			memcpy(wgToken1,wgToken2,sizeof(WGToken));
			memcpy(wgToken2,&wgTokenAux,sizeof(WGToken));
		}
		
		// comparison function
		inline static bool wgTokenEqual(WGToken *wgToken1, WGToken *wgToken2) {
			
			return (wgToken1->iWordSequence == wgToken2->iWordSequence);
		}
		
		// comparison function
		inline static bool wgTokenBetterPath(WGToken *wgToken1, WGToken *wgToken2) {
			
			return (wgToken1->fScore > wgToken2->fScore);
		}	
		
		// print a history item
		inline void printHistoryItem(HistoryItem *historyItem) {
		
			string strLexUnit;
			LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(historyItem->iLexUnitPron);
			m_lexiconManager->getStrLexUnitPronunciation(lexUnit,strLexUnit);
			
			int iWordSequence = -1;
			if (m_bLatticeGeneration) {
				iWordSequence = hashWordSequence(historyItem);
			}
		
			cout << setw(5) << historyItem->iEndFrame << " " << setw(20) << strLexUnit << " " 
				<< setw(10) << historyItem->fScore << " " << historyItem << " ws: " 
				<< iWordSequence << endl;
		}
		
		// print a token set
		inline void printWGToken(int iWGToken) {
		
			printWGToken(iWGToken+m_wgTokens);
		}
		
		// print a token set
		inline void printWGToken(WGToken *wgToken) {
		
			cout << "------------------------------------------------------" << endl;
			for(int i=0 ; i < m_iMaxWordSequencesState ; ++i) {
				if (wgToken[i].iWordSequence == -1) {
					break;
				}
				cout << setw(10) << wgToken[i].iWordSequence << " " << 
					FLT(12,4) << wgToken[i].fScore << " " << wgToken[i].iHistoryItem << endl;
			}
			cout << "------------------------------------------------------" << endl;
		}		
		
		// compute the load factor of the hash table containing unque word sequences (debugging)
		float computeLoadFactorHashWordSequences(int *iBucketsUsed, int *iCollisions);		
		
		// shows hash-occupation information (debugging)
		void printHashsStats();
		
		// print the hash-contents 
		void printHashContents();	
		
		// keeps the best history item for each unique word-sequence (auxiliar method)
		void keepBestHistoryItem(map<int,pair<float,int> > &mWSHistoryItem, int iHistoryItem);	
		

	public:
		
		// constructor
		DynamicDecoderX(PhoneSet *phoneSet, HMMManager *hmmManager, 
			LexiconManager *lexiconManager, LMManager *lmManager, float fLMScalingFactor, 
			DynamicNetworkX *dynamicNetwork, int iMaxActiveNodes, 
			int iMaxActiveNodesWE, int iMaxActiveTokensNode, float fBeamWidthNodes, 
			float fBeamWidthNodesWE, float fBeamWidthTokensNode, bool bWordGraphGeneration, 
			int iMaxWordSequencesState);

		// destructor
		~DynamicDecoderX();
		
		// initialization
		void initialize();
		
		// uninitialize
		void uninitialize();
		
		// begin utterance
		void beginUtterance();
		
		// end utterance
		void endUtterance();
		
		// process input feature vectors
		void process(MatrixBase<float> &mFeatures);	
		
		// return the BestPath
		BestPath *getBestPath();
		
		// return a hypothesis lattice for the utterance
		HypothesisLattice *getHypothesisLattice();
		
		// return the active lm-states at the current time (lm-state in active tokens)
		void getActiveLMStates(map<int,bool> &mLMState);	
		
};

};	// end-of-namespace

#endif
