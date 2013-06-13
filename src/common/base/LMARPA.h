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


#ifndef LMARPA_H
#define LMARPA_H

#include <iostream>

using namespace std;

#include <algorithm>
#include <string>

#include "LexiconManager.h"

namespace Bavieca {

// n-gram
typedef struct _NGram {
	int iLexUnit;								// lexical unit id
	float fProbability;						// probability
	float fProbabilityBackoff;				// backoff probability
	int iNGrams;								// number of higher order n-grams that come from this n-gram
	_NGram *ngrams;							// higher order n-grams that come from this n-gram
	_NGram *ngramBase;						// lower order n-gram from which this n-gram comes from 
} NGram;


/**
	@author root <dani.bolanos@gmail.com>
*/
class LMARPA {

	private:
	
		LexiconManager *m_lexiconManager;
		
		// properties
		string m_strFile;
		int m_iNGramOrder;
		vector<int> m_vNGrams;
		NGram **m_ngrams;
		int *m_iNGramsDropped;
		
		// whether the lm is loaded
		bool m_bLoaded;
		
		// convert lexical units read from the lm-file to upper case
		void lexUnitToUpperCase(string &str) {	
			
			if (str.compare(LEX_UNIT_BEG_SENTENCE) && str.compare(LEX_UNIT_END_SENTENCE) && str.compare(LEX_UNIT_UNKNOWN)) { 
				std::transform(str.begin(),str.end(),str.begin(),::toupper);
			}
		}	

		// sort the n-grams by lexUnit-id
		// not-very efficient: O(n^2), not good for very large vocabularies
		void sort(NGram *ngrams, int n);	
	
	public:

		// constructor
		LMARPA(LexiconManager *lexiconManager, const char *strFile);

		// destructor
		~LMARPA();
		
		// load
		void load();	
		
		// get ngrams
		NGram *getNGrams(int iOrder, int *iNGrams) {
		
			assert((iOrder >= 0) && (iOrder < m_iNGramOrder));
			*iNGrams = m_vNGrams[iOrder];
			
			return m_ngrams[iOrder];
		}
		
		// get ngrams
		int getNGrams(int iOrder) {
		
			assert((iOrder >= 0) && (iOrder <= m_iNGramOrder));
			return m_vNGrams[iOrder];
		}
		
		// print a ngram
		void printNGram(int *iId, int n);
		
		// print a ngram
		void print(NGram *ngram, bool bText = true) {
		
			assert(ngram);
			ostringstream oss;
			oss << "( ";
			NGram *ngramAux = ngram;
			while(ngramAux->ngramBase != NULL) {
				if (bText) {
					oss << m_lexiconManager->getStrLexUnit(ngramAux->iLexUnit);
				} else {
					oss << ngramAux->iLexUnit;
				}
				oss << " ";
				ngramAux = ngramAux->ngramBase;
			}
			oss << ")";
			cout << oss.str() << endl;
		}
		
		// return a ngram
		NGram *getNGram(int *iId, int n);	
		
		// return whether the lm is loaded
		bool lmLoaded() {
		
			return m_bLoaded;
		}
		
		// print language model information
		void print();
		
		// return the n-gram order
		int getNGramOrder() {
		
			return m_iNGramOrder;
		}
};

};	// end-of-namespace

#endif
