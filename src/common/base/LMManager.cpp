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


#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "LMManager.h"
#include "LogMessage.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
LMManager::LMManager(LexiconManager *lexiconManager, const char *strFile, 
	const char *strFormat, const char *strType, const char *strNGram) {

	m_lexiconManager = lexiconManager;
	m_strFile = strFile;
	m_strFormat = strFormat;
	m_strType = strType;
	
	// n-gram order
	if (strNGram != NULL) {
		m_iNGram = getNGram(strNGram);
	} else {
		m_iNGram = -1;
	}
	
	m_unigrams = NULL;
	m_bigrams = NULL;
	m_trigrams = NULL;
	
	m_iUnigrams = 0;
	m_iBigrams = 0;
	m_iTrigrams = 0;
	
	m_iUnigramsDropped = 0;
	m_iBigramsDropped = 0;
	m_iTrigramsDropped = 0;	
	
	m_states = NULL;
	m_arcs = NULL;
	
	m_bLMLoaded = false;
}

// destructor
LMManager::~LMManager() {
	destroy();
}

// load the language model from disk
void LMManager::load() {
	
	double dStartTime = TimeUtils::getTimeMilliseconds();

	// ARPA
	if (m_strFormat.compare(LM_FILE_FORMAT_ARPA) == 0) {
		loadARPA();
	} 
	// binary
	else {
		assert(m_strFormat.compare(LM_FILE_FORMAT_FSM) == 0);	
	} 
	
	double dEndTime = TimeUtils::getTimeMilliseconds();
	double dSeconds = (dEndTime-dStartTime)/1000;	
	cout << "LM loading time: " << dSeconds << "seconds\n";	
	
	m_bLMLoaded = true;
}

// load the language model from disk (ARPA format)
void LMManager::loadARPA() {

	string strLine;
	string strLexUnit;
	string strLexUnitPrev;
	string strLexUnitPrevPrev;
	string strLexUnitPrevPrevPrev;
	string strProbability;
	string strProbabilityBackoff;

	float fProbability;
	float fProbabilityBackoff;	
	int iId;
	int iIdPrev;
	int iIdPrevPrev;
	
	int iLine = 0;
	
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
	m_iUnigrams = 0;
	m_iBigrams = 0;
	m_iTrigrams = 0;
	
	int iNGram = 0;
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
		switch(iN) {
			case 1: {
				m_iUnigrams = iElements;
				break;
			} 
			case 2: {
				m_iBigrams = iElements;
				break;
			} 
			case 3: {
				m_iTrigrams = iElements;
				break;
			} 
			default: {
				BVC_ERROR << "ngram order not supported";
			}
		}
		++iNGram;
	}	
	
	// check that there are ngrams of the required order
	if (iNGram < m_iNGram) {
		BVC_ERROR << "not " << getStrNGram() << " entries found";
	} 
	
	if (m_iNGram == -1) {
		m_iNGram = iNGram;
	}
	
	// there has to be at least unigrams
	if (m_iUnigrams == 0) {
		BVC_ERROR << "no unigrams found in language model file";
	}	
	
	// (3) load the unigrams	
	
	// find the "unigram data section"
	bool bUnigramSection = false;
	while(std::getline(file.getStream(),strLine).good()) {
		++iLine;
		if (strLine.compare("\\1-grams:") == 0) {
			bUnigramSection = true;
			break;
		}
	}
	if (bUnigramSection == false) {
		BVC_ERROR << "no unigram data section defined in language model file";
	}
	
	// allocate memory for the unigrams
	m_unigrams = new Unigram[m_iUnigrams];
	for(int i = 0 ; i < m_iUnigrams ; ++i) {
		m_unigrams[i].iLexUnit = -1;
		m_unigrams[i].fProbability = -1;
		m_unigrams[i].fProbabilityBackoff = -1;
		m_unigrams[i].iBigrams = 0;
		m_unigrams[i].bigrams = NULL;
	}
	// read the unigrams one by one
	vector<string> vStrLexUnit;
	int iUnigramsInLexicon = m_iUnigrams;
	m_iUnigramsDropped = 0;
	for(int i = 0 ; i < iUnigramsInLexicon ; ) {
   	++iLine;
   	if (!std::getline(file.getStream(),strLine).good()) {
   		BVC_ERROR << "unigram expected at file \"" << m_strFile << "\", line: " << iLine;
		}
		std::stringstream s(strLine);
		IOBase::read(s,&fProbability,false);
		IOBase::readString(s,strLexUnit);
		fProbabilityBackoff = 0.0;
		IOBase::readWhiteSpaces(s);
		if (!s.eof()) {
			IOBase::read(s,&fProbabilityBackoff,false);
		}	
		// convert the lexical unit to upper case (just in case)
		lexUnitToUpperCase(strLexUnit);
		
		// is the word in the lexicon?
		int iId = m_lexiconManager->getLexUnitId(strLexUnit.c_str());
		// no: drop the unigram
		if (iId == -1) {
			BVC_WARNING << "dropped unigram \"" << strLexUnit << "\" at line: " << iLine << " (it is not in the lexicon)";
			++m_iUnigramsDropped;
			--iUnigramsInLexicon;
		}
		// yes: keep the unigram
		else {
			vStrLexUnit.push_back(strLexUnit);
			
			// in this case the lexical unit id can be used as the unigram index 
			// because all the ids are ordered and different
			m_unigrams[i].iLexUnit = i;
			m_unigrams[i].fProbability = fProbability;
			m_unigrams[i].fProbabilityBackoff = fProbabilityBackoff;	
			m_unigrams[i].iBigrams = 0;
			m_unigrams[i].bigrams = NULL;
			++i;
		}
	} 
	
	// make sure there is a unigram for <s> and </s>, otherwise create them
	Unigram *unigramBeg = &m_unigrams[m_lexiconManager->m_lexUnitBegSentence->iLexUnit];
	if (unigramBeg->iLexUnit == -1) {
		unigramBeg->iLexUnit = m_lexiconManager->m_lexUnitBegSentence->iLexUnit;
		unigramBeg->fProbability = 0.0;
		unigramBeg->fProbabilityBackoff = 0.0;
		unigramBeg->iBigrams = 0;
		unigramBeg->bigrams = NULL;
		++m_iUnigrams;
		BVC_WARNING << "unigram <s> not found, it was added to the language model";
	}
	Unigram *unigramEnd = &m_unigrams[m_lexiconManager->m_lexUnitEndSentence->iLexUnit];
	if (unigramEnd->iLexUnit == -1) {
		unigramEnd->iLexUnit = m_lexiconManager->m_lexUnitEndSentence->iLexUnit;
		unigramEnd->fProbability = 0.0;
		unigramEnd->fProbabilityBackoff = 0.0;
		unigramEnd->iBigrams = 0;
		unigramEnd->bigrams = NULL;
		++m_iUnigrams;
		BVC_WARNING << "unigram </s> not found, it was added to the language model";
	}	
	
	// check that lexical units in the lexicon fit the unigrams (they will if the lexicon was previously formatted)
	if (m_lexiconManager->checkConsistency(vStrLexUnit) == false) {
		BVC_WARNING << "lexicon and language model are inconsistent: lexicon entries will be reordered";
		// arrange the lexicon to fit the unigrams: (unigram i == lexical unit i)
		m_lexiconManager->arrangeLexUnits(vStrLexUnit);
	}
	vStrLexUnit.clear();
	
	// (4) load the bigrams (if any and if the system is configured to apply them) 	
	if ((m_iBigrams > 0) && (m_iNGram >= LM_NGRAM_BIGRAM)) {	
	
		// find the "bigram data section"
		bool bBigramSection = false;
		while(std::getline(file.getStream(),strLine).good()) {
			++iLine;
			if (strLine.compare("\\2-grams:") == 0) {
				bBigramSection = true;
				break;
			}
		}		
		if (bBigramSection == false) {
			BVC_ERROR << "no bigram data section defined in language model file";
		}
		
		// allocate memory for the bigrams
		m_bigrams = new Bigram[m_iBigrams];
		
		int iIdUnigram = NON_LEXUNIT_ID;		// initialized to a non-lexical unit id
		int iUnigramBigrams = 0;
		int iBigramsInLexicon = m_iBigrams;
		m_iBigramsDropped = 0;
		for(int i = 0 ; i < iBigramsInLexicon ; ++i) {
			++iLine;	
			if (!std::getline(file.getStream(),strLine).good()) {
				BVC_ERROR << "bigram expected at file \"" << m_strFile << "\", line: " << iLine;
			}
			std::stringstream s(strLine);
			IOBase::read(s,&fProbability,false);
			IOBase::readString(s,strLexUnitPrev);
			IOBase::readString(s,strLexUnit);
			fProbabilityBackoff = 0.0;
			IOBase::readWhiteSpaces(s);
			if (!s.eof()) {
				IOBase::read(s,&fProbabilityBackoff,false);
			}	
			
			// get the lexical unit id	
			iId = m_lexiconManager->getLexUnitId(strLexUnit);
			iIdPrev = m_lexiconManager->getLexUnitId(strLexUnitPrev);	
			
			//printf("%s %s %d %d\n",strLexUnitPrev,strLexUnit,iId,iIdPrev);
			
			// check that the lexical units are in the lexicon
			if ((iId == -1) || (iIdPrev == -1)) {
				BVC_WARNING << "dropped bigram (\"" << strLexUnitPrev << "\",\"" << strLexUnit << "\") at line: " << iLine << " (it is not in the lexicon)";
				++m_iBigramsDropped;
				--iBigramsInLexicon;
				--i;
				continue;
			}
			
			m_bigrams[i].iLexUnit = iId;
			m_bigrams[i].fProbability = fProbability;
			m_bigrams[i].fProbabilityBackoff = fProbabilityBackoff;
			m_bigrams[i].iTrigrams = 0;
			m_bigrams[i].trigrams = NULL;
			
			++iUnigramBigrams;
			
			if (iIdUnigram == NON_LEXUNIT_ID) {
				m_unigrams[iIdPrev].bigrams = &m_bigrams[i];
				iIdUnigram = iIdPrev;
			} else if (iIdPrev != iIdUnigram) {
				m_unigrams[iIdUnigram].iBigrams = iUnigramBigrams-1;
				iUnigramBigrams = 1;
				m_unigrams[iIdPrev].bigrams = &m_bigrams[i];
				iIdUnigram = iIdPrev;	
			}
		}
		if (iUnigramBigrams > 0) {
			m_unigrams[iIdPrev].iBigrams = iUnigramBigrams;
		}
		// check
		int iTotalBigrams = 0;
		for(int i = 0 ; i < m_iUnigrams ; ++i) {
			iTotalBigrams += m_unigrams[i].iBigrams;
		}
		assert(iTotalBigrams+m_iBigramsDropped == m_iBigrams);
	}
	
	// bigrams need to be sorted by lexUnitId, because binary search is used to locate them
	sortBigrams();
	
	// (4) load the trigrams (if any and the system is configured to apply them)	
	m_iTrigramsDropped = 0;
	if ((m_iTrigrams > 0) && (m_iNGram >= LM_NGRAM_TRIGRAM)) {	
	
		// find the "trigram data section"
		bool bTrigramSection = false;
		while(std::getline(file.getStream(),strLine).good()) {
			++iLine;
			if (strLine.compare("\\3-grams:") == 0) {
				bTrigramSection = true;
				break;
			}
		}		
		if (bTrigramSection == false) {
			BVC_ERROR << "no trigram data section defined in language model file: \"" << m_strFile << "\" line: " + iLine;
		}
		
		// allocate memory for the trigrams
		m_trigrams = new Trigram[m_iTrigrams];
		
		int iIdUnigram = NON_LEXUNIT_ID;		// initialized to a non-lexical unit id
		int iIdBigram = NON_LEXUNIT_ID;		// initialized to a non-lexical unit id
		int iBigramTrigrams = 0;	
		Bigram *bigram = NULL;
		
		int iTrigramsInLexicon = m_iTrigrams;
		m_iTrigramsDropped = 0;
		for(int i = 0 ; i < iTrigramsInLexicon ; ++i) {
			++iLine;
			if (!std::getline(file.getStream(),strLine).good()) {
				BVC_ERROR << "bigram expected at file \"" << m_strFile << "\", line: " << iLine;
			}
			std::stringstream s(strLine);
			IOBase::read(s,&fProbability,false);
			IOBase::readString(s,strLexUnitPrevPrev);
			IOBase::readString(s,strLexUnitPrev);
			IOBase::readString(s,strLexUnit);
			fProbabilityBackoff = 0.0;
			IOBase::readWhiteSpaces(s);
			if (!s.eof()) {
				IOBase::read(s,&fProbabilityBackoff,false);
			}	
			
			// get the lexical unit id	
			iId = m_lexiconManager->getLexUnitId(strLexUnit);
			iIdPrev = m_lexiconManager->getLexUnitId(strLexUnitPrev);
			iIdPrevPrev = m_lexiconManager->getLexUnitId(strLexUnitPrevPrev);			
			
			// check that the lexical units are in the lexicon
			if ((iId == -1) || (iIdPrev == -1) ||(iIdPrevPrev == -1)) {
				BVC_WARNING << "dropped trigram (\"" << strLexUnitPrevPrev << "\",\"" << strLexUnitPrev
				 << "\",\"" << strLexUnit << "\") at line: " << iLine << "(it is not in the lexicon)";
				++m_iTrigramsDropped;
				--i;
				--iTrigramsInLexicon;
				continue;
			}

			m_trigrams[i].iLexUnit = iId;
			m_trigrams[i].fProbability = fProbability;
			
			++iBigramTrigrams;
			
			if ((iIdUnigram == NON_LEXUNIT_ID) && (iIdBigram == NON_LEXUNIT_ID)) {
				bigram = getBigram(iIdPrevPrev,iIdPrev);
				if (bigram == NULL) {
					BVC_ERROR << "bigram (" << strLexUnitPrevPrev << "," << strLexUnitPrev 
					<< ") not defined but it is part of the trigram (" << strLexUnitPrevPrev 
					<< "," << strLexUnitPrev << "," << strLexUnit << ")";
				}
				bigram->trigrams = &m_trigrams[i];
				iIdUnigram = iIdPrevPrev;
				iIdBigram = iIdPrev;
			} else if ((iIdPrevPrev != iIdUnigram) || (iIdPrev != iIdBigram)) {
				bigram->iTrigrams = iBigramTrigrams-1;
				bigram = getBigram(iIdPrevPrev,iIdPrev);
				if (bigram == NULL) {
					BVC_ERROR << "bigram (" << strLexUnitPrevPrev << "," << strLexUnitPrev 
					<< ") not defined but it is part of the trigram (" << strLexUnitPrevPrev 
					<< "," << strLexUnitPrev << "," << strLexUnit << ")";
				}
				iBigramTrigrams = 1;
				bigram->trigrams = &m_trigrams[i];
				iIdUnigram = iIdPrevPrev;
				iIdBigram = iIdPrev;
			}	
		} 		
		if (iBigramTrigrams > 0) {
			bigram->iTrigrams = iBigramTrigrams;
		}
		// check
		int iTotalTrigrams = 0;
		for(int i = 0 ; i < m_iBigrams-m_iBigramsDropped ; ++i) {
			iTotalTrigrams += m_bigrams[i].iTrigrams;
		}
		assert(iTotalTrigrams+m_iTrigramsDropped == m_iTrigrams);
	}
	
	file.close();
}

// destroy the resources allocated for the language model
void LMManager::destroy() {

	// deallocate memory
	if (m_trigrams != NULL) {
		delete [] m_trigrams;
		this->m_trigrams = NULL;
	}
	if (m_bigrams != NULL) {
		delete [] m_bigrams;
		this->m_bigrams = NULL;
	}
	if (m_unigrams != NULL) {
		delete [] m_unigrams;
		this->m_unigrams = NULL;
	}
	
	if (m_states != NULL) {
		delete [] m_states;
		delete [] m_arcs;
	}	
	
	m_iUnigrams = 0;
	m_iBigrams = 0;
	m_iTrigrams = 0;	
	
	m_bLMLoaded = false;
}

// print language model information
void LMManager::print() {
	
	cout << "---- LM stats -----------------------------\n";
	cout << " file:        " << m_strFile.c_str() << "\n";
	cout << " format:      " << m_strFormat.c_str() << "\n";
	cout << " ngram order: " << getStrNGram() << "\n";
	if (m_iNGram >= LM_NGRAM_UNIGRAM) {
		cout << "  # unigrams: " << m_iUnigrams;
		if (m_iUnigramsDropped > 0) {
			cout << " (" << m_iUnigramsDropped << " not in lexicon)"; 
		}
		cout << "\n";
	} 
	if (m_iNGram >= LM_NGRAM_BIGRAM) {
		cout << "  # bigrams:  " << m_iBigrams;
		if (m_iBigramsDropped > 0) {
			cout << " (" << m_iBigramsDropped << " not in lexicon)"; 
		}
		cout << "\n";
	} 
	if (m_iNGram >= LM_NGRAM_TRIGRAM) {
		cout << "  # trigrams: " << m_iTrigrams;
		if (m_iTrigramsDropped > 0) {
			cout << " (" << m_iTrigramsDropped << " not in lexicon)"; 
		}
		cout << "\n";
	}
	cout << "-------------------------------------------------------------\n";
}

// compute a trigram language model score (log-likelihood)
float LMManager::computeTrigramScore(int iId1, int iId2, int iId3) {

	//printf("(%s,%s,%s)\n",m_lexiconManager->getStrLexUnit(iId1),m_lexiconManager->getStrLexUnit(iId2),m_lexiconManager->getStrLexUnit(iId3));

	// (1) check if p(id3|id1,id2) exists
	// check if the bigram (id1,id2) exists
	Bigram *bigram = getBigram(iId1,iId2);
	if (bigram == NULL) {
		// check if the bigram (id2,id3) exists
		Bigram *bigram2 = getBigram(iId2,iId3);
		// if the bigram does not exist, do the back-off: bo(id2) x p(id3)
		if (bigram2 == NULL) {
			return m_unigrams[iId2].fProbabilityBackoff + m_unigrams[iId3].fProbability;
		}	
		// return p(id3|id2)
		else {
			return bigram2->fProbability;
		}
	} 
	else {
		// check if the trigram (id1,id2,id3) exist
		Trigram *trigram = getTrigram(bigram,iId3);
		// if the trigram (id1,id2,id3) does not exist, do the back off: bo(id1,id2) x p(id3|id2)
		if (trigram == NULL) {
			// check if the bigram (id2,id3) exists
			Bigram *bigram2 = getBigram(iId2,iId3);
			// if the bigram does not exist, do the back-off: bo(id2) x p(id3) and return: bo(id1,id2) x bo(id2) x p(id3)
			if (bigram2 == NULL) {
				return (bigram->fProbabilityBackoff + (m_unigrams[iId2].fProbabilityBackoff + m_unigrams[iId3].fProbability));
			}
			// return bo(id1,id2) x p(id3|id2)
			else {
				return bigram->fProbabilityBackoff + bigram2->fProbability;
			}
		} else {
			// return the trigram score
			return trigram->fProbability;
		}
	}
}

// compute a bigram language model score (log-likelihood)
float LMManager::computeBigramScore(int iId1, int iId2) {

	// (1) check if p(id2|id1) exists
	Bigram *bigram = getBigram(iId1,iId2);
	// (2) if it does not, do the back off: bo(id1) x p(id2)
	if (bigram == NULL) {
		return (m_unigrams[iId2].fProbability + m_unigrams[iId1].fProbabilityBackoff);
	} else {
		return bigram->fProbability;
	}
}

// compute a unigram language model score (log-likelihood)
float LMManager::computeUnigramScore(int iId) {

	//assert(iId < m_iUnigrams);

	// p(id)
	return m_unigrams[iId].fProbability;
}

// sort the Bigrams by lexUnit-id (to prevent problems when spcial lexical units are not first in the language model)
void LMManager::sortBigrams() {

	// iterate through the unigrams
	for(int j=0 ; j < m_iUnigrams-1 ; ++j) {
		bool bSwapped = false;
		do {
			bSwapped = false;
			for(int i=0 ; i < m_unigrams[j].iBigrams-1 ; ++i) {
				if (m_unigrams[j].bigrams[i].iLexUnit > m_unigrams[j].bigrams[i+1].iLexUnit) {
					Bigram bigramAux;
					memcpy(&bigramAux,&m_unigrams[j].bigrams[i],sizeof(Bigram));
					memcpy(&m_unigrams[j].bigrams[i],&m_unigrams[j].bigrams[i+1],sizeof(Bigram));
					memcpy(&m_unigrams[j].bigrams[i+1],&bigramAux,sizeof(Bigram));
					bSwapped = true;
				}
			}
		} while(bSwapped);
	}
}


// convert n-gram descriptor to integer format
int LMManager::getNGram(const char *strNGram) {

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

// build the LM-graph
void LMManager::buildLMGraph() {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	unsigned int iStatesCreated = 0;
	unsigned int iArcsCreated = 0;
	
	LMStateTemp *states = NULL;
	LMStateTemp *stateInitial = NULL;
	LMStateTemp *stateFinal = NULL;
			
	// uniform: a single state with one transition per word, same probability
	if (m_iNGram == LM_NGRAM_ZEROGRAM) {
	
		m_iLMStateInitial = m_iLMStateFinal = 0;
		return;
	} 
	// unigram: a single state and one transition per word (no backoffs)
	else if (m_iNGram == LM_NGRAM_UNIGRAM) {
	
		m_iLMStateInitial = m_iLMStateFinal = 0;
		return;
	}
	// bigram: one state per unigram (+ unigram backoff state) and one transition per existing bigram + backoff transitions
	else if (m_iNGram == LM_NGRAM_BIGRAM) {
		
		unsigned int iStateID = 0;
	
		// get the unigrams
		int iUnigramsDiscarded = 0;
		
		// allocate memory to store the states (one state per unigram plus the bigram backoff)
		states = new LMStateTemp[m_iUnigrams+1];	
		for(int i=0 ; i < m_iUnigrams+1 ; ++i) {
			states[i].iState = i;
		}
		
		// create a map from unigram index (lexical unit index) to state
		LMStateTemp **stateMap = new LMStateTemp*[m_iUnigrams];
		for(int i=0 ; i < m_iUnigrams ; ++i) {
			stateMap[i] = NULL;
		}		
		
		// create the back-off state (to handle all the bigram back-offs)
		LMStateTemp *stateBackoff = &states[iStateID];
		stateBackoff->iState = iStateID++;
		++iStatesCreated;	
		
		// create one state per unigram
		for(int i=0 ; i < m_iUnigrams ; ++i) {
			// skip the unknown lexical unit
			if (i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
				++iUnigramsDiscarded;
				continue;
			}	
			
			// create the state
			LMStateTemp *state = &states[iStateID];
			state->iState = iStateID++;
			++iStatesCreated;
			stateMap[i] = state;
		
			// keep the final state (end of sentence)
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				stateFinal = state;
			}
		}
		
		// create an arc from the back-off state to all the unigram states
		for(int i=0 ; i < m_iUnigrams ; ++i) {
			// skip lexical units for which no state was created
			if (stateMap[i] == NULL) {
				continue;
			}
			// skip the <s>
			if (i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) {
				continue;
			}
			
			LMArcTemp *arc = NULL;
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				arc = newLMArcTemp(i,m_unigrams[i].fProbability,stateMap[i]);
			} else {
				arc = newLMArcTemp(i,m_unigrams[i].fProbability,stateMap[i]);
			}	
			++iArcsCreated;
			stateBackoff->lArc.push_back(arc);
			assert(stateBackoff != arc->stateDest);
		}
				
		// create arcs from each unigram state (unigram states are first in the vector)
		for(int i=0 ; i < m_iUnigrams ; ++i) {
			// skip lexical units for which no state was created
			if (stateMap[i] == NULL) {
				continue;
			}				
			// skip the </s> lexical unit (this is the final state)
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				continue;
			}
			
			// get the existing bigrams for this unigram
			int iBigrams = m_unigrams[i].iBigrams;
			Bigram *bigrams = m_unigrams[i].bigrams;
			int iBigramsDiscarded = 0;
			for(int j=0; j<iBigrams ; ++j) {
				// skip lexical units for which no state was created
				if (stateMap[bigrams[j].iLexUnit] == NULL) {
					++iBigramsDiscarded;
					continue;
				}
				// skip useless transition <s> <s>
				if ((i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) && (i == bigrams[j].iLexUnit)) {
					++iBigramsDiscarded;	
					continue;
				}
				
				// note that all the pronunciation variants have the same insertion penalty associated
				LMArcTemp *arc = NULL;
				if ((bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) || 	
					(bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) {	
					// no insertion penalty
					arc = newLMArcTemp(bigrams[j].iLexUnit,bigrams[j].fProbability,stateMap[bigrams[j].iLexUnit]);
				} else {
					// insertion penalty
					arc = newLMArcTemp(bigrams[j].iLexUnit,bigrams[j].fProbability,stateMap[bigrams[j].iLexUnit]);
				}
				++iArcsCreated;
				stateMap[i]->lArc.push_back(arc);
			}	
			// create a transition to the bigram backoff state (if necessary)
			if (iBigrams-iBigramsDiscarded < m_iUnigrams-iUnigramsDiscarded) {  // only if there are back-offs, this is the usual case
				// create a transition to the bigram back-off state
				LMArcTemp *arc = newLMArcTemp(BACKOFF_ARC,m_unigrams[i].fProbabilityBackoff,stateBackoff);
				stateMap[i]->lArc.push_back(arc);
				++iArcsCreated;
				assert(stateMap[i] != arc->stateDest);
			} else {
				assert(0);
			}
		}
		
		// create the initial state and connect it to the state (<s>)
		stateInitial = stateMap[m_lexiconManager->m_lexUnitBegSentence->iLexUnit];
		assert(stateInitial != NULL);
		
		delete [] stateMap;	
	}	
	// trigram: one state per existing bigram, one state per trigram backoff, one state for all the bigram backoffs
	else if (m_iNGram == LM_NGRAM_TRIGRAM) {
		
		unsigned int iStateID = 0;
	
		// get the unigrams (important: the array of unigrams has iVocabularySize index 
		// (so lexUnit Ids can be used to access it))
		// create a map from unigram index (lexical unit index) to state
		LMStateTemp **stateMapUnigram = new LMStateTemp*[m_iUnigrams];
		for(int i=0 ; i < m_iUnigrams ; ++i) {
			stateMapUnigram[i] = NULL;
		}
		
		// get the bigrams
		// create a map from bigram index (lexical unit index) to bigram state
		LMStateTemp **stateMapBigram = new LMStateTemp*[m_iBigrams];
		for(int i=0 ; i < m_iBigrams ; ++i) {
			stateMapBigram[i] = NULL;
		}		
	
		// allocate memory to store the states (one state per unigram and bigram plus the bigram backoff)
		states = new LMStateTemp[m_iUnigrams+m_iBigrams+1];	
		/*for(int i=0 ; i < m_iUnigrams+m_iBigrams+1 ; ++i) {
			states[i].iState = iStateID++;
		}*/	
		
		// map used to keep final states, no transitions will be created from these states
		//MStateBool mStateFinal;
	
		// TODO: there cannot be transitions from sentence marker states, for example </s> <s> cannot go to <s> <s> but end, the decoder is supposed to
		// be utterance based so the end of sentence will be the end of utterance, there cannot be end of sentence markers within the utterance
	
		//unsigned int iStateID = 0;
		
		// create the bigram back-off state (to handle all the bigram back-offs)
		LMStateTemp *stateBigramBackoff = &states[iStateID];
		stateBigramBackoff->iState = iStateID++;
		++iStatesCreated;
		
		// create the unigram states (one per unigram)
		LMStateTemp *stateFinalBackup = NULL;
		LMStateTemp *stateInitialBackup = NULL;
		for(int i=0 ; i < m_iUnigrams ; ++i) {	
		
			// -> skip the unknown lexical unit
			// -> the beginning of sentence bigram (</s>,<s>) or (<s>,<s>) is the initial state 
			// case 1) (</s>,<s>) -> (<s>) -> (w1) when the trigram (</s>,<s>,w1) does not exist
			// case 2) (</s>,<s>) -> (<s>,w1) when the trigram (</s>,<s>,w1) does exist
			// -> the end of sentence unigram (</s>) is the final state
			// case 1) (w1,w2) -> (w2) -> (</s>) when the trigram (w1,w2,</s>) does not exist
			// case 2) (w1,w2) -> (w2,</s>) -> (</s>) when the trigram (w1,w2,</s>) exists
			if (i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
				continue;
			}
			
			LMStateTemp *state = &states[iStateID];
			state->iState = iStateID++;
			++iStatesCreated;
			stateMapUnigram[i] = state;
			
			// keep a back-up initial state just in case there is no bigram state acting as initial
			if (i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) {
				stateInitialBackup = state;
			}
			// keep a back-up final state just in case there is no bigram state acting as final
			if (i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
				stateFinalBackup = state;
			}
			
			// create a transition to the bigram back-off (if necessary)
			if (m_unigrams[i].iBigrams < m_iUnigrams) {

				LMArcTemp *arc = newLMArcTemp(BACKOFF_ARC,m_unigrams[i].fProbabilityBackoff,stateBigramBackoff);
				++iArcsCreated;
				state->lArc.push_back(arc);
			}			
		}
			
		// create the bigram states (one per bigram)
		stateInitial = NULL;
		for(int i=0 ; i < m_iUnigrams ;  ++i) {
		
			// skip bigrams that come from the unknown lexical unit
			if (i == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
				continue;
			}
		
			for(int j=0 ; j < m_unigrams[i].iBigrams ; ++j) {
		
				// skip bigrams that go to the unknown lexical unit 
				if (m_unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit) {
					continue;
				}
				
				// skip bigrams that go to the start of sentence except the initial state (</s>,<s>) or (<s>,<s>)
				if ((m_unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) && 
					(i != m_lexiconManager->m_lexUnitEndSentence->iLexUnit) && 
					(i != m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					continue;
				}
				// skip bigrams that come from the end of sentence except the initial state: (</s> <s>)
				if ((i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) && 
					(m_unigrams[i].bigrams[j].iLexUnit != m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					continue;
				}
				// allow bigrams that come from the start of sentence except: (<s>,</s>) and (<s>,<unk>)
				if ((i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) && 
					((m_unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) ||
					(m_unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitUnknown->iLexUnit))) {
					continue;
				}
			
				// get the position of the bigram in the global array of bigrams
				int iIndex = &(m_unigrams[i].bigrams[j])-m_bigrams;
				assert((iIndex >= 0) && (iIndex < m_iBigrams));	
				
				LMStateTemp *state = &states[iStateID];
				state->iState = iStateID++;
				++iStatesCreated;
				stateMapBigram[iIndex] = state;
				
				// keep the initial state, which is (</s>,<s>) or (<s>,<s>), it can also be <s> if the others do not exist
				if (((i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) || 
					(i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit)) && 
					(m_unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					if (stateInitial != NULL) {
						BVC_ERROR << "multiple initial states found in language model!!";
					}
					stateInitial = state; 
				}
				
				// keep the final state, which is (</s>,<s>), it can also be just <s> if (</s>,<s>) does not exist
				if ((i == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) && 
					(m_unigrams[i].bigrams[j].iLexUnit == m_lexiconManager->m_lexUnitBegSentence->iLexUnit)) {
					if (stateFinal != NULL) {
						BVC_ERROR << "multiple final states found in language model!!";
					}
					stateFinal = state; 
				}
			}
		}
		// if there is no bigram acting as final state (</s>,<s>) use the unigram </s>
		if (stateFinal == NULL) {
			stateFinal = stateFinalBackup;
			if (stateFinal == NULL) {
				BVC_ERROR << "no final state </s> found in language model!!";
			}
		}	
		// if there is no bigram acting as initial state (<s>,<s>) or (<s>,<s>), use the unigram <s>
		if (stateInitial == NULL) {
			stateInitial = stateInitialBackup;
			if (stateInitial == NULL) {
				BVC_ERROR << "no initial state <s> found in language model!!";
			}
		}	
		
		// connect unigram states to brigram states
		for(int i=0 ; i < m_iUnigrams ; ++i) {	
		
			// skip non existent unigrams
			if (stateMapUnigram[i] == NULL) {
				continue;
			}
			
			// create a transition to the bigram states that come from this unigram	
			for(int j=0 ; j < m_unigrams[i].iBigrams ; ++j) {
			
				// get the position of the bigram in the global array of bigrams
				int iIndex = &(m_unigrams[i].bigrams[j])-m_bigrams;
				assert((iIndex >= 0) && (iIndex < m_iBigrams));	
				// skip non existent bigrams
				if (stateMapBigram[iIndex] == NULL) {
					continue;
				}	
				
				LMArcTemp *arc = newLMArcTemp(m_unigrams[i].bigrams[j].iLexUnit,
					m_unigrams[i].bigrams[j].fProbability,stateMapBigram[iIndex]);
				++iArcsCreated;
				stateMapUnigram[i]->lArc.push_back(arc);
			}
		}
		
		// create a transition from the bigram back-off state to all the unigram states
		for(int i=0 ; i < m_iUnigrams ; ++i) {
		
			// skip non existent unigrams
			if (stateMapUnigram[i] == NULL) {
				continue;
			}
			
			// do not create a back-off to the beginning of sentence
			if (i == m_lexiconManager->m_lexUnitBegSentence->iLexUnit) {
				continue;
			}
			
			LMArcTemp *arc = newLMArcTemp(i,m_unigrams[i].fProbability,stateMapUnigram[i]);
			++iArcsCreated;
			stateBigramBackoff->lArc.push_back(arc);
		}
		
		// create the transitions from each of the bigram states (bigram states are first in the vector)
		for(int i=0 ; i < m_iBigrams ; ++i) {
		
			// skip non-existent bigram states
			if (stateMapBigram[i] == NULL) {
				continue;
			}
			
			// connect the bigram state to the backup final state </s> if necessary
			if (stateFinal == stateFinalBackup) {
				if (m_bigrams[i].iLexUnit == m_lexiconManager->m_lexUnitEndSentence->iLexUnit) {
					LMArcTemp *arc = newLMArcTemp(m_lexiconManager->m_lexUnitBegSentence->iLexUnit,0.0,stateFinal);
					++iArcsCreated;
					stateMapBigram[i]->lArc.push_back(arc);
				}
			}
		
			// get the existing trigrams for this bigram
			int iTrigrams = m_bigrams[i].iTrigrams;
			Trigram *trigrams = m_bigrams[i].trigrams;
			for(int j=0; j<iTrigrams ; ++j) {
				// get the position of the bigram in the global array of bigrams
				Bigram *bigramAux = getBigram(m_bigrams[i].iLexUnit,trigrams[j].iLexUnit);
				if (bigramAux == NULL) {
				   const char *str1 = "x";
				   const char *str2 = m_lexiconManager->getStrLexUnit(m_bigrams[i].iLexUnit);
				   const char *str3 = m_lexiconManager->getStrLexUnit(trigrams[j].iLexUnit);
					BVC_WARNING << "unable to connect bigram (" << str1 << "," << str2 << ") to bigram (" << str2 << "," << str3 << ") through trigram (" << str1 << "," << str2 << "," << str3 << "), destination bigram was not found: trigram ignored!";
					continue;
				}	
				int iIndex = (bigramAux-m_bigrams);
				assert((iIndex >= 0) && (iIndex < m_iBigrams));
				
				// skip non existent bigrams states
				if (stateMapBigram[iIndex] == NULL) {
					continue;
				}
			
				LMArcTemp *arc = newLMArcTemp(trigrams[j].iLexUnit,trigrams[j].fProbability,stateMapBigram[iIndex]);
				++iArcsCreated;
				stateMapBigram[i]->lArc.push_back(arc);
			}
			// create a transition to the right unigram state (with the backoff link)
			if (iTrigrams < m_iUnigrams) {
			
				if (stateMapUnigram[m_bigrams[i].iLexUnit] != NULL) {
					LMArcTemp *arc = newLMArcTemp(BACKOFF_ARC,m_bigrams[i].fProbabilityBackoff,
						stateMapUnigram[m_bigrams[i].iLexUnit]);
					stateMapBigram[i]->lArc.push_back(arc);	
					++iArcsCreated;
				}	
			}
		}
		printf("states created: %d\n",iStatesCreated);
		/*for(int i=0 ; i < iStatesCreated ; ++i) {
			printStateTransitions(&states[i]);
		}*/
		
		// perform sanity checks to make sure that all the states created are connected (seen) and all the transitions too,
		bool *bStateProcessed = new bool[iStatesCreated];
		for(unsigned int i=0 ; i<iStatesCreated ; ++i) {
			bStateProcessed[i] = false;
		}
		unsigned int iStatesSeen = 0;
		unsigned int iTransitionSeen = 0;
		LLMStateTemp lState;
		lState.push_back(stateInitial);
		bStateProcessed[stateInitial->iState] = true;
		while(lState.empty() == false) {
			
			// get the next state to process
			LMStateTemp *stateFrom = lState.front();
			lState.pop_front();
			
			for(LLMArcTemp::iterator it = stateFrom->lArc.begin() ; it != stateFrom->lArc.end() ; ++it) {
			
				// insert G destination states into the queue (if not already processed)
				if (bStateProcessed[(*it)->stateDest->iState] == false) {
					lState.push_back((*it)->stateDest);
					bStateProcessed[(*it)->stateDest->iState] = true;
				}
				++iTransitionSeen;
			}
			++iStatesSeen;	
		}
		for(unsigned int i=0 ; i<iStatesCreated ; ++i) {
			if (bStateProcessed[i] == false) {
				BVC_WARNING << "lm-state not seen: " << i;
			}
			//assert(bStateProcessed[i] == true);
		}	
		assert(iStatesSeen == iStatesCreated);
		assert(iTransitionSeen == iArcsCreated);
		printf("building G: states created: %d states seen: %d\n",iStatesCreated,iStatesSeen);
		printf("building G: trans created: %d tran seen: %d\n",iArcsCreated,iTransitionSeen);
		delete [] bStateProcessed;
		
		delete [] stateMapUnigram;	
		delete [] stateMapBigram;	

		//*stateInitialG = stateInitial;
	}	
	// fourgram: one state per existing trigram + one state per existing bigram + one state per existing unigram
	else {
		assert(m_iNGram == LM_NGRAM_FOURGRAM);
		
		// not supported	
		BVC_ERROR << "fourgram language model not supported";
	}
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	
	printf("-- G graph -----------------------------\n");
	printf(" n-gram: %s\n",getStrNGram(m_iNGram));
	printf(" # states:      %12u\n",iStatesCreated);
	printf(" # transitions: %12u\n",iArcsCreated);
	printf(" building time: %.2f seconds\n",(dTimeEnd-dTimeBegin)/1000.0);
	if (stateFinal == NULL) {
		printf(" warning: no final state was found!!\n");
	} else {
		printf(" final state found\n");
	}
	printf("----------------------------------------\n");
		
	// compact the graph in order to optimize lookups
	
	// count arcs
	m_iArcs = 0;
	for(unsigned int i=0 ; i<iStatesCreated ; ++i) {
		m_iArcs += states[i].lArc.size();
	}
	m_iStates = iStatesCreated;
	
	// allocate memory
	m_states = new LMState[m_iStates+1];		// extra state is needed to mark the ending arc
	m_arcs = new LMArc[m_iArcs];
	
	// keep already created nodes (indexed by the corresponding temporal node id)
	int *iStateCreated = new int[m_iStates];
	int *iFirstArc = new int[m_iStates];
	for(int i=0 ; i<m_iStates ; ++i) {
		iStateCreated[i] = -1;
		iFirstArc[i] = -1;
	}

	// fill the structures
	int iState = 0;
	int iArc = 0;
	for(int iStateTemp = 0 ; iStateTemp < m_iStates ; ++iStateTemp) {
	
		// sort arcs by lexical unit so we can use a binary search for arc lookups
		states[iStateTemp].lArc.sort(LMManager::compareArcs);
	
		LMState *stateAux;
		int iArcBase = -1;
		if (iStateCreated[iStateTemp] == -1) {
			stateAux = &m_states[iState];
			iStateCreated[iStateTemp] = iState;
			iState++;
			iArcBase = iArc;
			iArc += states[iStateTemp].lArc.size();
		} else {		
			stateAux = &m_states[iStateCreated[iStateTemp]]; 
			iArcBase = iFirstArc[iStateCreated[iStateTemp]];
		}
		stateAux->iArcBase = iArcBase;
		// store the outgoing arcs
		int i=0;
		for(LLMArcTemp::iterator jt = states[iStateTemp].lArc.begin() ; jt != states[iStateTemp].lArc.end() ; ++jt,++i) {
			if (iStateCreated[(*jt)->stateDest->iState] == -1) {
				iStateCreated[(*jt)->stateDest->iState] = iState;
				iFirstArc[iState] = iArc;
				iArc += (*jt)->stateDest->lArc.size();
				iState++;
			}
			m_arcs[iArcBase+i].iLexUnit = (*jt)->iLexUnit;
			m_arcs[iArcBase+i].fScore = (*jt)->fScore;
			m_arcs[iArcBase+i].iStateDest = iStateCreated[(*jt)->stateDest->iState];	
		}
	}
	assert(iArc == m_iArcs);
	assert(iState == m_iStates);
	m_states[m_iStates].iArcBase = m_iArcs;
	
	m_iLMStateInitial = iStateCreated[stateInitial->iState];
	m_iLMStateFinal = iStateCreated[stateFinal->iState];
	
	// clean-up
	for(unsigned int i=0 ; i<iStatesCreated ; ++i) {
		for(LLMArcTemp::iterator it = states[i].lArc.begin() ; it != states[i].lArc.end() ; ++it) {
			delete *it;
		}
	}	
	delete [] states;	
	delete [] iStateCreated;
	delete [] iFirstArc;	
}

// get the initial state
int LMManager::getInitialState() {

	return m_iLMStateInitial;
}

// update the language model state with the given lexical unit and returns the new lm state
int LMManager::updateLMState(int iLMStatePrev, int iLexUnit, float *fScore) {

	switch(m_iNGram) {
		// zerogram (uniform distribution), only one state
		case LM_NGRAM_ZEROGRAM: {
			*fScore = 0.0;
			return 0;
		}	
		// unigram, just one state
		case LM_NGRAM_UNIGRAM: {
			*fScore = m_unigrams[iLexUnit].fProbability;
			return 0;
		}
		// higher-order n-grams
		case LM_NGRAM_BIGRAM:
		case LM_NGRAM_TRIGRAM:
		case LM_NGRAM_FOURGRAM: {
		
			//m_lexiconManager->print(m_lexiconManager->getLexUnit(iLexUnit));
		
			assert((iLMStatePrev >= 0) && (iLMStatePrev < m_iStates));
			LMState *state = &m_states[iLMStatePrev];
			*fScore = 0.0;
			int iPasses = 0;
			
			while(1) {
			
				int iFirst = state->iArcBase;
				int iLast = (state+1)->iArcBase-1;
				int iMiddle;
				while(iFirst <= iLast) {
					iMiddle = (iFirst+iLast)/2;
					if (m_arcs[iMiddle].iLexUnit == iLexUnit) {
						*fScore += m_arcs[iMiddle].fScore;
						return m_arcs[iMiddle].iStateDest;				
					} else if (m_arcs[iMiddle].iLexUnit < iLexUnit) {
						iFirst = iMiddle+1;
					} else {
						iLast = iMiddle-1;
					}
				}	
				
				LMArc *arcBackoff = &m_arcs[(state+1)->iArcBase-1];
				/*if (arcBackoff->iLexUnit != BACKOFF_ARC) {
					m_lexiconManager->print(m_lexiconManager->getLexUnit(iLexUnit));
				}*/
				assert(arcBackoff->iLexUnit == BACKOFF_ARC);
				*fScore += arcBackoff->fScore;
				state = &m_states[arcBackoff->iStateDest];
					
				++iPasses;
			}
		}
		default: {
			assert(0);
			return -1;
		}
	}
}

// return the score resulting from moving to the given lm-state to the final state
float LMManager::toFinalState(int iLMState) {

	// zerogram, uniform distribution
	if (m_iNGram == LM_NGRAM_ZEROGRAM) {	
		return 0.0;
	} 
	// unigram
	else if (m_iNGram == LM_NGRAM_UNIGRAM) {
		return m_unigrams[m_lexiconManager->m_lexUnitEndSentence->iLexUnit].fProbability;
	} 
	// bigram
	else if (m_iNGram == LM_NGRAM_BIGRAM) {
		float fScore = -FLT_MAX;
		updateLMState(iLMState,m_lexiconManager->m_lexUnitEndSentence->iLexUnit,&fScore);
		return fScore;
	} 
	// trigram
	else if (m_iNGram == LM_NGRAM_TRIGRAM) {
		float fScore1 = -FLT_MAX;
		int iLMState1 = updateLMState(iLMState,m_lexiconManager->m_lexUnitEndSentence->iLexUnit,&fScore1);
		// sometimes the final state is just </s>
		if (iLMState1 == m_iLMStateFinal) {
			return fScore1;	
		}	
		// other times it is </s> <s>
		float fScore2 = -FLT_MAX;
		int iLMState2 = updateLMState(iLMState1,m_lexiconManager->m_lexUnitBegSentence->iLexUnit,&fScore2);
		assert(iLMState2 == m_iLMStateFinal);
		
		return fScore1+fScore2;
	} 
	// not supported
	else {
		assert(0);
		return -FLT_MAX;
	}
}

// return language model scores for all words in the vocabulary for a given LM-state (word history)
void LMManager::getLMScores(int iLMState, float *fLMScores, int iVocabularySize) {

	//double dTimeBegin = TimeUtils::getTimeMilliseconds();

	LMArc **lmArcBackoff = new LMArc*[m_iNGram-1];
	int iEmpty = m_iNGram-2;
	
	// get the lm-state for each n-gram order (includes actual lm-state plus backoff lm-states)
	LMState *state = m_states+iLMState;
	LMArc *arcBackoff = m_arcs+((state+1)->iArcBase-1);
	while(arcBackoff->iLexUnit == BACKOFF_ARC) {
		assert(iEmpty >= 0);
		lmArcBackoff[iEmpty--] = arcBackoff;
		arcBackoff = m_arcs+(((m_states+arcBackoff->iStateDest)+1)->iArcBase-1);
	}
	
	// initialization (debug purposes)
	for(int i=0 ; i < iVocabularySize; ++i) {
		fLMScores[i] = FLT_MAX;
	}
	
	// (1) unobserved n-grams (backoffs)
	// deal with back-off transitions in ascending order (starting from unigrams, then bigrams, etc)
	// if a higher-order n-gram exists then its score overwrites the score of lower n-gram backoff
	int iComputed = 0;
	for(int i=iEmpty+1 ; i < m_iNGram-1 ; ++i) {
	
		// get accumulated backoff score from higher order n-grams
		float fBackoffAcc = 0.0;
		for(int j=i ; j < m_iNGram-1 ; ++j) {
			fBackoffAcc += lmArcBackoff[j]->fScore;
		}
		// compute lm-scores for all observed n-grams in this n-gram table
		LMState *state = m_states+lmArcBackoff[i]->iStateDest;
		LMArc *arcFinal = m_arcs+(state+1)->iArcBase;
		if (i != iEmpty+1) {		// backoff arc is the last arc, stop before it
			arcFinal--;
		}
		for(LMArc *arc = m_arcs+state->iArcBase; arc != arcFinal ; ++arc) {
			assert((arc->iLexUnit >= 0) && (arc->iLexUnit < iVocabularySize));
			fLMScores[arc->iLexUnit] = fBackoffAcc+arc->fScore;	
			++iComputed;
		}
		//printf("computed: %d\n",iComputed);
		iComputed = 0;
	}
	// (2) observed n-grams
	assert((m_arcs+((state+1)->iArcBase)-1)->iLexUnit == BACKOFF_ARC);
	for(LMArc *arc = m_arcs+state->iArcBase ; arc != m_arcs+((state+1)->iArcBase)-1 ; ++arc) {
		assert((arc->iLexUnit >= 0) && (arc->iLexUnit < iVocabularySize));	
		fLMScores[arc->iLexUnit] = arc->fScore;	
		++iComputed;
	}	
	//printf("computed: %d\n",iComputed);
	
	// set lm-score for filler units and sentence markers (0.0)
	// that kind of lexical units are after the lexical units for all unigrams
	for(int i=m_iUnigrams ; i < iVocabularySize ; ++i) {
		fLMScores[i] = 0.0;
		//printf("lex: %s\n",m_lexiconManager->getStrLexUnit(i));
	}
	// there may or maynot be unigrams for the unknown symbol and the sentence markers 
	fLMScores[m_lexiconManager->m_lexUnitUnknown->iLexUnit] = 0.0;
	fLMScores[m_lexiconManager->m_lexUnitBegSentence->iLexUnit] = 0.0;
	fLMScores[m_lexiconManager->m_lexUnitEndSentence->iLexUnit] = 0.0;

	// check
	for(int i=0 ; i < iVocabularySize; ++i) {
		assert(fLMScores[i] != FLT_MAX);
	}
	
	delete [] lmArcBackoff;

	//double dTimeEnd = TimeUtils::getTimeMilliseconds();
	//printf("seconds: %12.8fs\n",(dTimeEnd-dTimeBegin)/1000.0);
}

};	// end-of-namespace
