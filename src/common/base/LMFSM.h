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


#ifndef LMFSM_H
#define LMFSM_H

using namespace std;

#if defined __linux__ || defined __APPLE__ || __MINGW32__
#include <tr1/unordered_map>
#elif _MSC_VER
#include <hash_map>
#else 
	#error "unsupported platform"
#endif

#include <list>
#include <string>
#include <vector>

#include "LexiconManager.h"
#include "LMARPA.h"
#include "LogMessage.h"

namespace Bavieca {

struct _LMArcTemp;
typedef vector<_LMArcTemp*> VLMArcTemp;
typedef list<_LMArcTemp*> LLMArcTemp;

#define BACKOFF_ARC			INT_MAX		// backoff arcs are used to connect to back-off states

// temporal language model state (each state is connected to other states via epsilon or a word)
typedef struct {
	int iState;	
	LLMArcTemp lArc;
	//NGram *ngram;
} LMStateTemp;

typedef list<LMStateTemp*> LLMStateTemp;

// temporal language model arc (keeps the score)
typedef struct _LMArcTemp {
	int iLexUnit;								// lexical unit
	float fScore;								// language model score
	LMStateTemp *stateDest;
} LMArcTemp;

// language model state (each state is connected to other states via epsilon or a word)
typedef struct {	
	int iArcBase;								// index of the first outgoing arc
} LMState;

// language model arc (keeps the score)
typedef struct {
	int iLexUnit;								// lexical unit
	float fScore;								// language model score
	int iStateDest;
} LMArc;


#if defined __linux__ || defined __APPLE__ || __MINGW32__
typedef std::tr1::unordered_map<string,LMStateTemp*> MNGramState;
#elif _MSC_VER
typedef std::unordered_map<string,LMStateTemp*> MNGramState;
#endif

/**
	@author root <dani.bolanos@gmail.com>
*/
class LMFSM {

	private:
	
		LexiconManager *m_lexiconManager;
		bool m_bLoaded;
	
		// lm properties
		int m_iNGramOrder;
		int *m_iNGrams;
		
		// language model states and arcs
		int m_iLMStateInitial;
		int m_iLMStateFinal;
		LMState *m_states;
		LMArc *m_arcs;
		int m_iStates;
		int m_iArcs;
		int m_iArcsStandard;
		int m_iArcsBackoff;	
		
		// compare two arcs by lexical unit index
		static bool compareArcs(const LMArcTemp *arc1, const LMArcTemp *arc2) {
			
			return (arc1->iLexUnit < arc2->iLexUnit);
		}
		
		// return a new temporal language model arc
		LMArcTemp *newArc(int iLexUnit, float fScore, LMStateTemp *stateDest) {
		
			LMArcTemp *arc = new LMArcTemp;
			arc->iLexUnit = iLexUnit;
			arc->fScore = fScore;
			arc->stateDest = stateDest;
			
			if (iLexUnit != BACKOFF_ARC) {
				++m_iArcsStandard;
			} else {
				++m_iArcsBackoff;
			}
	
			return arc;
		}
		
		// compute a hash key
		string hashKey(NGram *ngram, bool bIgnoreLower = false, int iLexUnitInit = -1) {
		
			int iLen = sizeof(int)*(m_iNGramOrder-1)+2;
			char *strKey = new char[iLen];
			strKey[0] = '#';	// for the backoff
			int i=1;
			if (iLexUnitInit != -1) {
				memcpy(strKey+i,&iLexUnitInit,sizeof(int));
				i+=sizeof(int);
			}
			assert(ngram);
			NGram *ngramAux = ngram;
			while(ngramAux->ngramBase && (!bIgnoreLower || ngramAux->ngramBase->ngramBase)) {	
				memcpy(strKey+i,&ngramAux->iLexUnit,sizeof(int));
				i+=sizeof(int);
				ngramAux = ngramAux->ngramBase;
			}
			strKey[i] = 0;
			assert(i<iLen);
			string str(strKey,i);
			delete [] strKey;
			
			return str;
		}
	
	public:	

		//constructor
		LMFSM(LexiconManager *lexiconManager);

		// destructor
		~LMFSM();
		
		// return whether the lm is loaded
		bool lmLoaded() {
		
			return m_bLoaded;
		}	
		
		// build the FSM
		void build(LMARPA *lmARPA);
		
		// perform sanity checks to make sure that all the states/transitions created are connected
		void checkConnected(LMStateTemp *stateInitial, LMStateTemp *stateFinal, 
			LMStateTemp *stateBackoffZerogram, int iStates, int iArcs);	
		
		// compact the FSM to use less memory and speed-up lookups (better locality)
		void compact(LMStateTemp *states, int iStates, int iArcs, LMStateTemp *stateInitial, LMStateTemp *stateFinal);
		
		// store to disk
		void store(const char *strFile);
		
		// load from disk
		void load(const char *strFile);

		// get the initial state
		int getInitialState();		
		
		// update the language model state with the given lexical unit and returns the new lm state
		int updateLMState(int iLMStatePrev, int iLexUnit, float *fScore);
		
		// return the score resulting from moving to the given lm-state to the final state
		float toFinalState(int iLMState);
		
		// return language model scores for all words in the vocabulary for a given LM-state (word history)
		// typically used for language model look-ahead
		void getLMScores(int iLMState, float *fLMScores, int iVocabularySize);
		
		// compute the likelihood of the given sequence of word
		float computeLikelihood(const char *str);	
		
		// return the n-gram order
		int getNGramOrder() {
			
			return m_iNGramOrder;
		}
		
		// print
		void print();
		
};

};	// end-of-namespace

#endif
