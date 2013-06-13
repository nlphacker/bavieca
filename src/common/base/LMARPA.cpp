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


#include "LMARPA.h"
#include "LMManager.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "LogMessage.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
LMARPA::LMARPA(LexiconManager *lexiconManager, const char *strFile) {

	m_lexiconManager = lexiconManager;
	m_strFile = strFile;
	m_iNGramOrder = -1;
	m_bLoaded = false;
}

// destructor
LMARPA::~LMARPA()
{
	if (m_bLoaded) {
		for(int i=0 ; i <= m_iNGramOrder ; ++i) {
			delete [] m_ngrams[i]; 
		}
		delete [] m_ngrams;
		delete [] m_iNGramsDropped;
		m_bLoaded = false;
	}
}

// load from disk
void LMARPA::load() {

	string strLine;
	float fProbability;
	float fProbabilityBackoff;	
	int iLine = 0;
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	FileInput file(m_strFile.c_str(),false);
	file.open();
	
	// (1) skip all the lines until the "data" section is found
	bool bDataSection = false;
	while(std::getline(file.getStream(),strLine).good()) {	
		++iLine;
		if (strLine.compare("\\data\\") == 0) {
			bDataSection = true;
			break;
		}
	}
	if (bDataSection == false) {
		BVC_ERROR<< "no data section defined in language model file";	
	}
	
	// (2) read the number of n-grams
	m_iNGramOrder = 0;
	m_vNGrams.push_back(1);	// zerograms
	while(1) {
		++iLine;
		if (!std::getline(file.getStream(),strLine).good()) {
			BVC_ERROR << "ngram counter expected in file: \"" << m_strFile << "\" line: " << iLine;
		}
		const size_t iFound = strLine.find("ngram");
		if (iFound == string::npos) {
			break;
		}	
		std::stringstream s(strLine);
		string str;
		IOBase::readString(s,str);
		IOBase::readString(s,str);
		const size_t iEq = str.find("=");
		int iN = atoi(str.substr(0,iEq).c_str());
		int iElements = atoi(str.substr(iEq+1,str.length()-iEq).c_str());
		++m_iNGramOrder;
		if (iN != m_iNGramOrder) {
			BVC_ERROR << "unexpected n-gram order in file: \"" << m_strFile << "\" line: " << iLine;
		}
		if (iElements <= 0) {
			BVC_ERROR << "non-positive number of n-grams: \"" << m_strFile << "\" line: " << iLine;
		}
		m_vNGrams.push_back(iElements);
	}	
	// there has to be at least unigrams
	if (m_iNGramOrder == 0) {
		BVC_ERROR << "no n-grams found in language model file";
	}
	
	int *iId = new int[m_iNGramOrder+1];
	string *strLexUnit = new string[m_iNGramOrder+1];
	m_ngrams = new NGram*[m_iNGramOrder+1];
	m_iNGramsDropped = new int[m_iNGramOrder+1];
	for(int i=0 ; i <= m_iNGramOrder ; ++i) {
		m_iNGramsDropped[i] = 0;
	}
	int *iIdCurrent = new int[m_iNGramOrder+1];
	
	// deal with the zerogram
	m_ngrams[0] = new NGram[1];
	m_ngrams[0][0].iLexUnit = -1;
	m_ngrams[0][0].fProbability = 0.0;
	m_ngrams[0][0].fProbabilityBackoff = 0.0;
	m_ngrams[0][0].ngrams = NULL;
	m_ngrams[0][0].iNGrams = 0;	
	m_ngrams[0][0].iNGrams = m_vNGrams[1];
	m_ngrams[0][0].ngramBase = NULL;	
	
	// load the n-grams
	for(int i=1 ; i <= m_iNGramOrder ; ++i) {
		
		// find the "n-gram data section"
		bool bNGramSection = false;
		while(std::getline(file.getStream(),strLine).good()) {
			++iLine;
			ostringstream oss;
			oss << "\\" << i << "-grams:"; 
			if (strLine.compare(oss.str()) == 0) {
				bNGramSection = true;
				break;
			}
		}
		if (bNGramSection == false) {
			BVC_ERROR << "no n-gram (" << i << ") data section defined in language model file";
		}
		
		// allocate memory for the n-grams
		m_ngrams[i] = new NGram[m_vNGrams[i]];
		if (i == 1) {
			m_ngrams[0][0].ngrams = m_ngrams[1];
		}
		NGram *ngrams = m_ngrams[i];
		
		// read the n-grams
		//int iIdCurrent = NON_LEXUNIT_ID;
		for(int j = 0 ; j < m_iNGramOrder+1 ; ++j) {
			iIdCurrent[j] = NON_LEXUNIT_ID;
		}
		NGram *ngramBase = NULL;
		int iNGram = 0;
		int iSubNGrams = 0;
		for(int j=0 ; j < m_vNGrams[i] ; ++j) {
			
			++iLine;
			if (!std::getline(file.getStream(),strLine).good()) {
				BVC_ERROR << "n-gram expected at file \"" << m_strFile << "\", line: " << iLine;
			}
			std::stringstream s(strLine);
			IOBase::read(s,&fProbability,false);
			for(int k=0 ; k < i ; ++k) {
				IOBase::readString(s,strLexUnit[k]);
				lexUnitToUpperCase(strLexUnit[k]);			// convert to upper case (just in case)
			}
			fProbabilityBackoff = 0.0;
			IOBase::readWhiteSpaces(s);
			if (!s.eof()) {
				IOBase::read(s,&fProbabilityBackoff,false);
			}
			
			// are the words in the lexicon?
			bool bFound = true;
			for(int k=0 ; k < i ; ++k) {
				iId[k] = m_lexiconManager->getLexUnitId(strLexUnit[k].c_str());
				if (iId[k] == -1) {
					BVC_WARNING << "dropped n-gram \"" << strLexUnit[k] << "\" at line: " << iLine << " (it is not in the lexicon)";
					++m_iNGramsDropped[i+1];
					bFound = false;
					break;
				}
			}
			if (bFound) {	
				if (i > 1) { // except for unigrams	
					bool bDiff = false;
					for(int k=0 ; k <= i-2 ; ++k) {
						if (iIdCurrent[k] != iId[k]) {
							bDiff = true;
							break;
						}
					}
					if (bDiff) {
						if (ngramBase) {
							ngramBase->iNGrams = iSubNGrams;
							// sort n-grams by lexUnit-id (needed for binary search)
							sort(ngramBase->ngrams,ngramBase->iNGrams);
						}
						ngramBase = getNGram(iId,i-1);
						if (!ngramBase) {
							printNGram(iId,i-1);
						}
						assert(ngramBase);	
						ngramBase->ngrams = ngrams+iNGram;
						iSubNGrams = 0;
						for(int k=0 ; k <= i-2 ; ++k) {
							iIdCurrent[k] = iId[k];
						}
					}
					assert(ngramBase);
				}	
				ngrams[iNGram].iLexUnit = iId[i-1];
				ngrams[iNGram].fProbability = fProbability;
				ngrams[iNGram].fProbabilityBackoff = fProbabilityBackoff;	
				ngrams[iNGram].iNGrams = 0;
				ngrams[iNGram].ngrams = NULL;
				if (i > 1) {
					ngrams[iNGram].ngramBase = ngramBase;
				}
				// unigrams 
				else {
					ngrams[iNGram].ngramBase = m_ngrams[0];	// backoff
				}
				++iNGram;
				++iSubNGrams;
			}
		}
		if (i > 1) {
			ngramBase->iNGrams = iSubNGrams;
		} else {
			m_ngrams[0][0].iNGrams = m_vNGrams[1]-m_iNGramsDropped[1];
			// sort unigrams by lexUnit-id (needed for binary search)
			sort(m_ngrams[0][0].ngrams,m_ngrams[0][0].iNGrams);	
		}
		// check
		if (i > 1) {
			int iTotalNGrams = 0;
			for(int j = 0 ; j < m_vNGrams[i-2] ; ++j) {
				iTotalNGrams += m_ngrams[i-2][j].iNGrams;
				if (m_ngrams[i-2][j].iNGrams) {
					assert(m_ngrams[i-2][j].ngrams);
				}
			}
			assert(iTotalNGrams+m_iNGramsDropped[i-1] == m_vNGrams[i-1]);	
		}
	}
	file.close();	
	
	// make sure there are initial and final unigram states (<s> and </s>)
	if (getNGram(&m_lexiconManager->m_lexUnitBegSentence->iLexUnit,1) == NULL) {
		BVC_ERROR << "beginning of sentence unigram (<s>) was not found";
	}
	if (getNGram(&m_lexiconManager->m_lexUnitEndSentence->iLexUnit,1) == NULL) {
		BVC_ERROR << "end of sentence unigram (</s>) was not found";
	}	
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	
	BVC_VERB << "-- language model ARPA -------------------";
	BVC_VERB << " ngram order: " << LMManager::getStrNGram(m_iNGramOrder);
	BVC_VERB << " file: " << m_strFile.c_str();
	int iNGrams = 0;
	for(int i=0 ; i <= m_iNGramOrder ; ++i) {
		BVC_VERB << " # " << LMManager::getStrNGram(i) << "s: " << m_vNGrams[i];
		iNGrams += m_vNGrams[i];
		if (m_iNGramsDropped[i] > 0) {
			BVC_VERB << " (" << m_iNGramsDropped[i] << " not in lexicon)"; 
		}
	}
	BVC_VERB << " # ngrams: " << iNGrams;	
	BVC_VERB << " loading time: " << (dTimeEnd-dTimeBegin)/1000.0 << " seconds";	
	BVC_VERB << "------------------------------------------";
	
	delete [] iId;
	delete [] strLexUnit;
	delete [] iIdCurrent;

	m_bLoaded = true;
}

// sort the n-grams by lexUnit-id
// it is not efficient: O(n^2), not good for very large vocabularies
// however is lexicon and language model entries are sorted alphabetically
// the time spent in sorting is minimal
void LMARPA::sort(NGram *ngrams, int n) {

	bool bSwapped;
	do {
		bSwapped = false;
		for(int i=0 ; i < n-1 ; ++i) {
			if (ngrams[i].iLexUnit > ngrams[i+1].iLexUnit) {
				NGram ngramAux;
				memcpy(&ngramAux,&ngrams[i],sizeof(NGram));
				memcpy(&ngrams[i],&ngrams[i+1],sizeof(NGram));
				memcpy(&ngrams[i+1],&ngramAux,sizeof(NGram));
				bSwapped = true;
			}
		}
	} while(bSwapped);
}

// print a ngram
void LMARPA::printNGram(int *iId, int n) {

	cout << "( ";
	for(int i=0 ; i < n ; ++i) {	
		cout << m_lexiconManager->getStrLexUnit(iId[i]) << " ";
	}
	cout << ")" << endl;
}

// return a ngram (binary search)
NGram *LMARPA::getNGram(int *iId, int n) {

	assert((iId) && (n));
	
	NGram *ngram = m_ngrams[0];
	for(int i=0 ; i < n ; ++i) {
		
		bool bFound = false;
		int iLexUnit = iId[i];
		int iFirst = 0;
		int iLast = ngram->iNGrams-1;
		int iMiddle;
		while(iFirst <= iLast) {
			iMiddle = (iFirst+iLast)/2;
			if (ngram->ngrams[iMiddle].iLexUnit == iLexUnit) {
				ngram = &(ngram->ngrams[iMiddle]);
				bFound = true;
				break;
			} else if (ngram->ngrams[iMiddle].iLexUnit < iLexUnit) {
				iFirst = iMiddle+1;
			} else {
				iLast = iMiddle-1;
			}
		}	
		
		if (!bFound) {
			return NULL;	
		}
	}
	
	return ngram;
}

// print language model information
void LMARPA::print() {
	
	cout << "---- LM stats -----------------------------" << endl;
	cout << " file: " << m_strFile.c_str() << endl;
	cout << " ngram order: " << LMManager::getStrNGram(m_iNGramOrder) << endl;
	int iNGrams = 0;
	for(int i=0 ; i <= m_iNGramOrder ; ++i) {
		cout << " # " << LMManager::getStrNGram(i+1) << "s: " << m_vNGrams[i];
		iNGrams += m_vNGrams[i];
		if (m_iNGramsDropped[i] > 0) {
			cout << " (" << m_iNGramsDropped[i] << " not in lexicon)"; 
		}
	}
	cout << " # ngrams: " << iNGrams << endl;
	cout << "-------------------------------------------------------------" << endl;
}

};	// end-of-namespace