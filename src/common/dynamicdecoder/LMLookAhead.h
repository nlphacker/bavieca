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


#ifndef LMLOOKAHEAD_H
#define LMLOOKAHEAD_H

class DynamicDecoderX;
class DynamicNetworkX;
class LexiconManager;
class LMManager;

using namespace std;

#include <string>

namespace Bavieca {

// note: the hash table is organized so the first N entries of the table are the buckets and the
//       last M entries are used to handle collisions using a linked list from the bucket-entries
//       If S is the maximum number of elements the hash-table needs to hold, then M must be S-1
//       to deal with the worst case scenario (all elements have the same hash key)

// entry in the hash table
typedef struct _HashEntry {
	int iLMState;						// word-history (language model state)
	float *fLAScores;					// look-ahead scores
	int iNext;							// next table entry (to handle collisions)
} LAHashEntry;

/**
	@author root <dani.bolanos@gmail.com>
*/
class LMLookAhead {

	private:
	
		LexiconManager *m_lexiconManager;
		LMManager *m_lmManager;
		DynamicDecoderX *m_dynamicDecoderX;
		int m_iVocabularySize;
		bool m_bInitialized;
		
		// look-ahead tree
		int m_iLANodes;		// # nodes
		int *m_iLATree;		// actual tree (in topological order), each position keeps index of predecessor, or -1
		
		// cache (hash-table) of look-ahead scores for different word histories		
		unsigned int m_iCacheElementsMax;		// maximum number of look-ahead trees that can be kept simultaneously
		unsigned int m_iHashBucketEntries;		// number of bucket-entries in the hash-table
		unsigned int m_iHashCollisionEntries;	// number of collision-entries in the hash-table
		unsigned int m_iHashEntries;				// total number of entries in the hash-table (bucket + collision entries)
		int m_iHashEntryCollisionAvailable;		// index of next entry available for collisions
		LAHashEntry *m_hashEntries;				// hash table
		
		// cache stats
		unsigned int m_iCacheHits;
		unsigned int m_iCacheMisses;
		
		// create a tree of look-ahead scores for the given lm-state (word-history)
		float *computeLAScores(int iLMState);
		
		// make sure the cache is clean (debug)
		void checkCache();
		
		// clean the hash-table used as a cache
		void cacheGarbageCollection();
		
		// destroy
		void destroy();	

		// print cache stats (performance analysis)
		void printCacheStats();
		
		// clear a cache entry
		inline void clear(int iEntry) {
		
			m_hashEntries[iEntry].iLMState = -1;
			m_hashEntries[iEntry].fLAScores = NULL;
			m_hashEntries[iEntry].iNext = -1;
		}

	public:

		// constructor
		LMLookAhead(LexiconManager *lexiconManager, LMManager *lmManager, 
			DynamicNetworkX *dynamicNetwork, DynamicDecoderX *dynamicDecoderX, int iCacheElementsMax);

		// destructor
		~LMLookAhead();
		
		// initialize
		void initialize();		
		
		// return a look-ahead tree for the given lm-state (word-history)
		float *getLAScores(int iLMState);	

};

};	// end-of-namespace

#endif
