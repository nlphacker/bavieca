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


#ifndef LMMANAGER_H
#define LMMANAGER_H

#include "LexiconManager.h"
#include "Global.h"

namespace Bavieca {

class FileInput;
class FileOutput;
class IOBase;
class LogMessage;

// language model format
#define LM_FILE_FORMAT_ARPA		"ARPA"
#define LM_FILE_FORMAT_FSM			"FSM"					// Finite State Machine

// language model types
#define LM_TYPE_NGRAM			0
#define LM_TYPE_CFG				1

// ngram types
#define LM_NGRAM_ZEROGRAM		0			// (uniform)
#define LM_NGRAM_UNIGRAM		1
#define LM_NGRAM_BIGRAM			2
#define LM_NGRAM_TRIGRAM		3
#define LM_NGRAM_FOURGRAM		4

// n-gram types (text format)
#define LM_NGRAM_TXT_ZEROGRAM		"zerogram"
#define LM_NGRAM_TXT_UNIGRAM		"unigram"
#define LM_NGRAM_TXT_BIGRAM		"bigram"
#define LM_NGRAM_TXT_TRIGRAM		"trigram"
#define LM_NGRAM_TXT_FOURGRAM		"fourgram"

using namespace std;

#include <string>

// look-ahead cache element
/*typedef struct _LACacheElement {
	float *fScores;							// look-ahead scores
	int iTimeLastAccess;						// time of the last access	
	int iTokenReferences;					// number of tokens that reference the element
	_LACacheElement *nextActive;			// next element in the list of active elements
	_LACacheElement *prevActive;			// next element in the list of active elements
	//_LACacheElement *next;				// next element in the pool
	_LACacheElement **lookupReference;	// address of the pointer that references this element from the look-up array
} LACacheElement;*/

typedef struct {
	int iLexUnit;								// lexical unit id
	float fProbability;						// probability
} Trigram;

typedef struct {
	int iLexUnit;								// lexical unit id
	float fProbability;						// probability
	float fProbabilityBackoff;				// backoff probability
	int iTrigrams;								// number of trigrams that come from this bigram
	Trigram *trigrams;						// linked list of trigrams that come from this bigram
} Bigram;

typedef struct {
	int iLexUnit;								// lexical unit id
	float fProbability;						// probability
	float fProbabilityBackoff;				// backoff probability
	int iBigrams;								// number of bigrams coming from this unigram
	Bigram *bigrams;							// bigrams coming from this unigram
} Unigram;

struct _LMArcTemp;
typedef vector<_LMArcTemp*> VLMArcTemp;
typedef list<_LMArcTemp*> LLMArcTemp;

#define BACKOFF_ARC			INT_MAX		// backoff arcs are used to connect to back-off states

// temporal language model state (each state is connected to other states via epsilon or a word)
typedef struct {
	int iState;	
	LLMArcTemp lArc;
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

/**
	@author Daniel Bolanos <bolanos@cslr.colorado.edu>
*/
class LMManager {

	private:
	
		LexiconManager *m_lexiconManager;
		string m_strFile;
		string m_strFormat;
		string m_strType;	
		int m_iType;						// language model type (ngram, CFG, etc)
		int m_iNGram;
		
		// language model data		
		int m_iUnigrams;
		int m_iBigrams;
		int m_iTrigrams;
		Unigram *m_unigrams;
		Bigram *m_bigrams;
		Trigram *m_trigrams;
		
		// droped n-grams (those that are not in the lexicon)
		int m_iUnigramsDropped;
		int m_iBigramsDropped;
		int m_iTrigramsDropped;
		
		// whether the lm is loaded
		bool m_bLMLoaded;		
		
		// language model states
		int m_iLMStateInitial;
		int m_iLMStateFinal;
		LMState *m_states;
		LMArc *m_arcs;
		int m_iStates;
		int m_iArcs;
		
		// compare two arcs by lexical unit index
		static bool compareArcs(LMArcTemp *arc1, LMArcTemp *arc2) {
			
			return (arc1->iLexUnit < arc2->iLexUnit);
		}
		
   public:
   
      // constructor
      LMManager(LexiconManager *lexiconManager, const char *strFile, const char *strFormat,
       const char *strType, const char *strNGram);

      // destructor 
      ~LMManager();
      
      int getNGram() {
      
      	return m_iNGram;
      }
      
      const char *getStrNGram() {
      
      	return getStrNGram(m_iNGram);
      }      
      
      static const char *getStrNGram(unsigned char iNGram) {
      
      	switch(iNGram) {
      		case LM_NGRAM_ZEROGRAM: {
      			return LM_NGRAM_TXT_ZEROGRAM;
      		}
      		case LM_NGRAM_UNIGRAM: {
      			return LM_NGRAM_TXT_UNIGRAM;
      		}
      		case LM_NGRAM_BIGRAM: {
      			return LM_NGRAM_TXT_BIGRAM;
      		}
      		case LM_NGRAM_TRIGRAM: {
      			return LM_NGRAM_TXT_TRIGRAM;
      		}
      		case LM_NGRAM_FOURGRAM: {
      			return LM_NGRAM_TXT_FOURGRAM;
      		}
      		default: {
      			return NULL;
      		}
      	}
      }
      
		// load the language model from disk
		void load(); 
      
		// load in ARPA format
		void loadARPA();
		
		// store as a Finite State Machine
		bool storeFSM(const char *strFile);	
		
		// load as a Finite State Machine
		void loadFSM(const char *strFile);	
		
		// desotry the language model
		void destroy();
		 
		// print stats
		void print();
      
		// compute a trigram language model score (log-likelihood)
		float computeTrigramScore(int iId1, int iId2, int iId3);	

		// compute a bigram language model score (log-likelihood)
		float computeBigramScore(int iId1, int iId2);	

		// compute a unigram language model score (log-likelihood)
		float computeUnigramScore(int iId);	

		// get a reference to the unigrams
		Unigram *getUnigrams(int &iUnigrams) {
		
			iUnigrams = m_iUnigrams;
		
			return m_unigrams;
		}
		
		// get the unigram
		Unigram *getUnigram(int iLexUnit) {
		
			return &(m_unigrams[iLexUnit]);
		}
		
		// get a reference to the bigrams
		Bigram *getBigrams(int &iBigrams) {
		
			iBigrams = m_iBigrams;
		
			return m_bigrams;
		}

		// get a bigram (if exists)
		Bigram *getBigram(int iLexUnit1, int iLexUnit2) {
		
			// binary search (the list of bigrams is sorted)
			Bigram *bigrams = m_unigrams[iLexUnit1].bigrams;
			if (bigrams == NULL) {
				return NULL;
			}
			
			int iFirst = 0;
			int iLast = m_unigrams[iLexUnit1].iBigrams-1;
			int iMiddle;
			while(iFirst <= iLast) {
				iMiddle = (iFirst+iLast)/2;
				if (bigrams[iMiddle].iLexUnit == iLexUnit2) {
					return &(bigrams[iMiddle]);
				} else if (bigrams[iMiddle].iLexUnit < iLexUnit2) {
					iFirst = iMiddle+1;
				} else {
					iLast = iMiddle-1;
				}		
			}		
			
			return NULL;			
		}

		// get a reference to the trigrams
		Trigram *getTrigrams(int &iTrigrams) {
		
			iTrigrams = m_iTrigrams;
		
			return m_trigrams;
		}
		
		// get a trigram (if it exists)
		Trigram *getTrigram(Bigram *bigram, int iLexUnit) {
		
			for(int i=0 ; i<bigram->iTrigrams ; ++i) {
				if (bigram->trigrams[i].iLexUnit == iLexUnit) {
					return &bigram->trigrams[i];
				}
			}
			
			return NULL;
		}
		
		// return whether the lm has already been loaded
		bool lmLoaded() {
		
			return m_bLMLoaded;
		}
		
		// sort the Bigrams by lexUnit-id (to prevent problems when spcial lexical units are not first in the language model)
		void sortBigrams();
		
		// convert n-gram descriptor to integer format
		static int getNGram(const char *strNGram);
		
		// convert lexical units read from the lm-file to upper case
		void lexUnitToUpperCase(string &strLexUnit) {	
			
			string str = "";
			if ((strLexUnit.compare("<s>") != 0) && (strLexUnit.compare("</s>") != 0)) { 
				for(unsigned int i=0 ; i < strLexUnit.length() ; ++i) {
					str += toupper(strLexUnit.c_str()[i]);
				}
				strLexUnit = str;
			}
		}
		
		// convert lexical units read from the lm-file to upper case
		void lexUnitToUpperCase(char *str) {	
		}
		
		// build the LM-graph
		void buildLMGraph();	
		
		// return a new temporal language model state
		LMStateTemp *newLMStateTemp(int iState) {
		
			LMStateTemp *state = new LMStateTemp;
			state->iState = iState;
		
			return state;
		}
		
		// return a new temporal language model arc
		LMArcTemp *newLMArcTemp(int iLexUnit, float fScore, LMStateTemp *stateDest) {
		
			LMArcTemp *arc = new LMArcTemp;
			arc->iLexUnit = iLexUnit;
			arc->fScore = fScore;
			arc->stateDest = stateDest;
		
			return arc;
		}
		
		// get the initial state
		int getInitialState();		
		
		// update the language model state with the given lexical unit and returns the new lm state
		int updateLMState(int iLMStatePrev, int iLexUnit, float *fScore);
		
		// return the score resulting from moving to the given lm-state to the final state
		float toFinalState(int iLMState);
		
		// return language model scores for all words in the vocabulary for a given LM-state (word history)
		void getLMScores(int iLMState, float *fLMScores, int iVocabularySize);		

};

};	// end-of-namespace

#endif
