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


#include "DynamicDecoderX.h"
#include "DynamicNetworkX.h"
#include "LMLookAhead.h"
#include "LexiconManager.h"
#include "LMManager.h"

namespace Bavieca {

// constructor
LMLookAhead::LMLookAhead(LexiconManager *lexiconManager, LMManager *lmManager, 
	DynamicNetworkX *dynamicNetwork, DynamicDecoderX *dynamicDecoderX, int iCacheElementsMax)
{
	m_lexiconManager = lexiconManager;
	m_lmManager = lmManager;
	m_dynamicDecoderX = dynamicDecoderX;
	//m_iCacheElementsMax = iCacheElementsMax;	
	m_iCacheElementsMax = 40000;	
	m_iLANodes = -1;
	m_iLATree = dynamicNetwork->getLMLookAheadTree(&m_iLANodes);
	m_bInitialized = false;	
}

// destructor
LMLookAhead::~LMLookAhead()
{
	destroy();
}

// initialize
void LMLookAhead::initialize() {

	m_iVocabularySize = m_lexiconManager->getVocabularySize();
	
	// cache (hash-table) of look-ahead scores for different word histories		
	m_iHashBucketEntries = m_iCacheElementsMax*2;
	m_iHashCollisionEntries = m_iCacheElementsMax;	// worst scenario implies storing all elements-1 as collisions
	m_iHashEntries = m_iHashBucketEntries + m_iHashCollisionEntries;
	m_hashEntries = new LAHashEntry[m_iHashEntries];
	for(unsigned int i=0 ; i < m_iHashEntries ; ++i) {
		m_hashEntries[i].iLMState = -1;
		m_hashEntries[i].fLAScores = NULL;
		m_hashEntries[i].iNext = -1;
	}
	m_iHashEntryCollisionAvailable = m_iHashBucketEntries;
	int iCacheSizeBytes = m_iHashEntries*sizeof(LAHashEntry);
	BVC_VERB << "look-ahead cache size: " << iCacheSizeBytes << " bytes";
	
	// cache stats
	m_iCacheHits = 0;
	m_iCacheMisses = 0;	
	
	m_bInitialized = true;
}

// destroy
void LMLookAhead::destroy() {

	assert(m_bInitialized);
	
	//destroy the cache
	for(unsigned int i=0 ; i < m_iHashEntries ; ++i) {
		if (m_hashEntries[i].fLAScores) {
			delete [] m_hashEntries[i].fLAScores;
		}
	}
	delete [] m_hashEntries;
	
	m_bInitialized = false;
}

// create a tree of look-ahead scores for the given lm-state (word-history)
float *LMLookAhead::computeLAScores(int iLMState) {

	//return new float[1];

	float *fLAScores = new float[m_iVocabularySize+m_iLANodes];
	
	// get the LM-scores for the word history
	m_lmManager->getLMScores(iLMState,fLAScores,m_iVocabularySize);
	
	// compute the look-ahead scores (loop over the array, which contains a topological order of the tree)
	// - each position in the array keeps the index of the predecessor in the look-ahead tree
	// (1) keep the maximum lm-score at any node in the look-ahead tree
	/*for(int i=0 ; i < m_iLANodes ; ++i) {
		if (m_iLATree[i] != -1) {
			fLAScores[m_iLATree[i]] = max(fLAScores[m_iLATree[i]],fLAScores[i]);
		}
	}
	// (2) keep the delta lm-score between every node and its predecessor
	for(int i=0 ; i < m_iLANodes ; ++i) {
		if (m_iLATree[i] != -1) {
			//fLAScores[i] -= fLAScores[m_iLATree[i]];
		}
	}*/	

	return fLAScores;
}

// return a tree of look-ahead scores for the given lm-state (word-history)
float *LMLookAhead::getLAScores(int iLMState) {

	// compute the hash-key
	unsigned int iEntry = iLMState % m_iHashBucketEntries;
	LAHashEntry &entry = m_hashEntries[iEntry];
	
	if (((m_iCacheHits+m_iCacheMisses) % 10000) == 0) {
		printCacheStats();
	}	
	
	// (1) empty entry
	if (entry.iLMState == -1) {	
		entry.iLMState = iLMState;
		entry.fLAScores = computeLAScores(iLMState);
		++m_iCacheMisses;
		return entry.fLAScores;
	} 
	// (2) hit, reuse entry
	else if (entry.iLMState == iLMState) {	
		assert(entry.fLAScores);
		++m_iCacheHits;
		return entry.fLAScores;	
	} 
	// (3) collision, look for the lm-state in the list of collision-entries linked to the bucket-entry
	else {
		// clean cache?
		if (m_iHashEntryCollisionAvailable == (int)m_iHashEntries) {
			cacheGarbageCollection();
			assert(m_iHashEntryCollisionAvailable != (int)m_iHashEntries);
		}
		int *iHashEntryAux = &entry.iNext;
		while(*iHashEntryAux != -1) {
			LAHashEntry &entryAux = m_hashEntries[*iHashEntryAux];
			// hit, reuse entry
			if (entryAux.iLMState == iLMState) {
				assert(entryAux.fLAScores);
				++m_iCacheHits;
				return entryAux.fLAScores;
			}
			iHashEntryAux = &entryAux.iNext;
		}
		// add a new a collision to the linked list of collisions
		assert(*iHashEntryAux == -1);
		*iHashEntryAux = m_iHashEntryCollisionAvailable;	
		assert(m_iHashEntryCollisionAvailable != (int)m_iHashEntries);
		LAHashEntry *entryCollision = m_hashEntries+*iHashEntryAux;
		entryCollision->iLMState = iLMState;
		entryCollision->fLAScores = computeLAScores(iLMState);
		entryCollision->iNext = -1;
		++m_iHashEntryCollisionAvailable;
		++m_iCacheMisses;	
		return entryCollision->fLAScores;	
	}
}

// clean the hash-table used as a cache
// - garbage collection is only called when the meximum cache entries are used
// - lm-states corresponding to active tokens are retrieved and their cache entries cleared
// - this procedure should be faster than keeping track of used entries at all times
void LMLookAhead::cacheGarbageCollection() {

	// (1) get list of active lm-states (lm-states from active tokens)
	map<int,bool> mLMState;
	m_dynamicDecoderX->getActiveLMStates(mLMState);	
	
	// (2) clear unused entries
	int iEntriesCleared = 0;
	for(int i=0 ; i < m_iHashEntryCollisionAvailable ; ++i) {
		if ((m_hashEntries[i].iLMState != -1) && (mLMState.find(m_hashEntries[i].iLMState) == mLMState.end())) {
			delete [] m_hashEntries[i].fLAScores;
			clear(i);
			++iEntriesCleared;
		}
	}
	
	if (iEntriesCleared == 0) {
		// cache needs to be resized, so far just thrown an error
		BVC_ERROR << "look-ahead cache is full, space for more entries is needed";
	}
	
	// (3) reorganize surviving collisions
	int iEntryLast = m_iHashEntryCollisionAvailable;
	m_iHashEntryCollisionAvailable = m_iHashBucketEntries;
	for(int i=m_iHashBucketEntries ; i < iEntryLast ; ++i) {	
		// skip cleared collision entries
		if (m_hashEntries[i].iLMState == -1) {
			continue;
		}
		unsigned int iEntry = m_hashEntries[i].iLMState % m_iHashBucketEntries;
		// (a) use the bucket
		if (m_hashEntries[iEntry].iLMState == -1) {
			m_hashEntries[iEntry].iLMState = m_hashEntries[i].iLMState;
			m_hashEntries[iEntry].fLAScores = m_hashEntries[i].fLAScores;
			m_hashEntries[iEntry].iNext = -1;
			clear(i);
		}
		// (b) add to the linked list of collisions
		else {
			int *iHashEntryAux = &m_hashEntries[iEntry].iNext;
			while(*iHashEntryAux != -1) {
				iHashEntryAux = &m_hashEntries[*iHashEntryAux].iNext;
			}
			*iHashEntryAux = m_iHashEntryCollisionAvailable;
			assert(m_iHashEntryCollisionAvailable != (int)m_iHashEntries);
			LAHashEntry *entryCollision = m_hashEntries+*iHashEntryAux;
			entryCollision->iLMState = m_hashEntries[i].iLMState;
			entryCollision->fLAScores = m_hashEntries[i].fLAScores;
			entryCollision->iNext = -1;
			if (*iHashEntryAux != i) {
				clear(i);
			}
			++m_iHashEntryCollisionAvailable;
			if (*iHashEntryAux != i) {
				assert(m_hashEntries[m_iHashEntryCollisionAvailable].iLMState == -1);
			}
		}
	}
	
	// report error if there is no space for more collissions (although there are some buckets empty)
	if (m_iHashEntryCollisionAvailable == (int)m_iHashEntries) {
		// cache needs to be resized, so far just thrown an error
		BVC_ERROR << "look-ahead cache is full, " << iEntriesCleared << " buckets are empty, but there is no space to handle collisions";
	}	
	
	BVC_VERB << "entries cleared: " << iEntriesCleared;
	BVC_VERB << "collision available: " << m_iHashEntryCollisionAvailable;
}

// print cache stats (performance analysis)
void LMLookAhead::printCacheStats() {

	// count used entries
	unsigned int iBucketsUsed = 0;
	unsigned int iCollisions = 0;
	for(unsigned int i=0 ; i < m_iHashBucketEntries ; ++i) {
		if (m_hashEntries[i].iLMState != -1) {
			++iBucketsUsed;
		}
	}
	for(unsigned int i=m_iHashBucketEntries ; i < m_iHashEntries ; ++i) {
		if (m_hashEntries[i].iLMState != -1) {
			++iCollisions;
		}
	}
	
	BVC_VERB << "-- look-ahead cache stats --" ;
	BVC_VERB << "cache hits: " << m_iCacheHits;
	BVC_VERB << "cache miss: " << m_iCacheMisses;	
	BVC_VERB << "total entries:    " << m_iHashEntries;
	BVC_VERB << "- buckets used:   " << iBucketsUsed;
	BVC_VERB << "- collisions:     " << iCollisions;
	BVC_VERB << "load factor:      " << ((float)iBucketsUsed)/((float)m_iHashBucketEntries);
	BVC_VERB << "collision factor: " << ((float)iCollisions)/((float)iBucketsUsed);
	BVC_VERB << "----------------------------";
}

};	// end-of-namespace

