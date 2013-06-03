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

#include "Global.h"
#include "LexiconManager.h"

namespace Bavieca {

class FileInput;
class FileOutput;
class IOBase;
class LogMessage;
class LMFSM;
class LMARPA;

// language model format
#define LM_FILE_FORMAT_ARPA		"ARPA"				// ARPA format
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

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class LMManager {

	private:
	
		LexiconManager *m_lexiconManager;	// pronunciation lexicon
		string m_strFile;							// lm file 
		string m_strFormat;						// lm format
		string m_strType;							// lm type
		int m_iType;								// lm type (ngram, CFG, etc)
		int m_iNGram;								// n-gram order
		
		// whether the lm is loaded
		bool m_bLMLoaded;
		
		// FSM
		LMFSM *m_lmFSM;
		// ARPA
		LMARPA *m_lmARPA;
		
   public:
   
      // constructor
      LMManager(LexiconManager *lexiconManager, const char *strFile, const char *strFormat, const char *strType);

      // destructor 
      ~LMManager();
      
		// load the language model
		void load();
      
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

		// convert n-gram descriptor to integer format
		static int getNGram(const char *strNGram) {
		
			if (strcmp(strNGram,LM_NGRAM_TXT_ZEROGRAM) == 0) {
				return LM_NGRAM_ZEROGRAM;
			} 
			else if (strcmp(strNGram,LM_NGRAM_TXT_UNIGRAM) == 0) {
				return LM_NGRAM_UNIGRAM;
			} 
			else if (strcmp(strNGram,LM_NGRAM_TXT_BIGRAM) == 0) {
				return LM_NGRAM_BIGRAM;
			} 
			else if (strcmp(strNGram,LM_NGRAM_TXT_TRIGRAM) == 0) {
				return LM_NGRAM_TRIGRAM;
			} 
			else {
				assert(strcmp(strNGram,LM_NGRAM_TXT_FOURGRAM) == 0);
				return LM_NGRAM_FOURGRAM;
			} 
		}
		
		// return the lm FSM
      LMFSM *getFSM() {
 
      	return m_lmFSM;
      }
		
		// return the lm FSM
      LMARPA *getARPA() {
 
      	return m_lmARPA;
      }
		
};

};	// end-of-namespace

#endif
