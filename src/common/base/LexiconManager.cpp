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
#include "IOBase.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "LogMessage.h"

namespace Bavieca {

// constructor
LexiconManager::LexiconManager(const char *strFile, PhoneSet *phoneSet)
{
   m_strFile = strFile;
   m_phoneSet = phoneSet;
   m_mLexUnit = new MLexUnit;
   
	// number of lexical units of each kind
	m_iLexicalUnitsStandard = 0;		
	m_iLexicalUnitsFiller = 0;		
	m_iLexicalUnitsSentenceDelimiter = 0;		
	m_iLexicalUnitsUnknown = 0;
	m_iLexicalUnitsUnique = 0;				// number of lexical units regardless alternative pronunciations	
	
	// these lexical units are allways defined (needed by the decoder)
	m_lexUnitBegSentence = NULL;
	m_lexUnitEndSentence = NULL;
	m_lexUnitUnknown = NULL;
}

// destructor
LexiconManager::~LexiconManager()
{
	destroy();
}

// load the lexicon from a file
void LexiconManager::load() {

	try {
	
		int iLexUnitId = 0;	
		
		FileInput file(m_strFile.c_str(),false);
		file.open();
		
		// (1) insert the sentence delimiters and the unknown symbol (if not in the language model they will be removed later)
		
		// unknown (before the word identity is known)
		m_lexUnitUnknown = new LexUnit;
		m_lexUnitUnknown->iType = LEX_UNIT_TYPE_UNKNOWN;
		m_lexUnitUnknown->iLexUnit = iLexUnitId;
		m_lexUnitUnknown->iPronunciation = 0;
		m_lexUnitUnknown->fProbability = 1.0;
		m_lexUnitUnknown->fInsertionPenalty = 0.0;
		LexUnitX *lexUnitX = new LexUnitX;
		lexUnitX->iLexUnit = iLexUnitId++;
		lexUnitX->strLexUnit = LEX_UNIT_UNKNOWN;
		lexUnitX->vLexUnitPronunciations.push_back(m_lexUnitUnknown);
		m_lexiconX.push_back(lexUnitX);	
		m_lexicon.push_back(m_lexUnitUnknown);
		m_mLexUnit->insert(MLexUnit::value_type(lexUnitX->strLexUnit,lexUnitX));
		++m_iLexicalUnitsUnknown;
			
		// end of sentence "</s>" (must be created before the beginning of sentence to preserve alphabetical order)
		m_lexUnitEndSentence = new LexUnit;
		m_lexUnitEndSentence->iType = LEX_UNIT_TYPE_SENTENCE_DELIMITER;
		m_lexUnitEndSentence->iLexUnit = iLexUnitId;	
		m_lexUnitEndSentence->iPronunciation = 0;
		m_lexUnitEndSentence->fProbability = 1.0;
		m_lexUnitEndSentence->fInsertionPenalty = 0.0;
		lexUnitX = new LexUnitX;
		lexUnitX->iLexUnit = iLexUnitId++;
		lexUnitX->strLexUnit = LEX_UNIT_END_SENTENCE;
		lexUnitX->vLexUnitPronunciations.push_back(m_lexUnitEndSentence);
		m_lexiconX.push_back(lexUnitX);
		m_lexicon.push_back(m_lexUnitEndSentence);
		m_mLexUnit->insert(MLexUnit::value_type(lexUnitX->strLexUnit,lexUnitX));
		++m_iLexicalUnitsSentenceDelimiter;
		
		// beginning of sentence "<s>"
		m_lexUnitBegSentence = new LexUnit;
		m_lexUnitBegSentence->iType = LEX_UNIT_TYPE_SENTENCE_DELIMITER;
		m_lexUnitBegSentence->iLexUnit = iLexUnitId;
		m_lexUnitBegSentence->iPronunciation = 0;
		m_lexUnitBegSentence->fProbability = 1.0;
		m_lexUnitBegSentence->fInsertionPenalty = 0.0;
		lexUnitX = new LexUnitX;
		lexUnitX->iLexUnit = iLexUnitId++;
		lexUnitX->strLexUnit = LEX_UNIT_BEGINNING_SENTENCE;
		lexUnitX->vLexUnitPronunciations.push_back(m_lexUnitBegSentence);
		m_lexiconX.push_back(lexUnitX);
		m_lexicon.push_back(m_lexUnitBegSentence);
		m_mLexUnit->insert(MLexUnit::value_type(lexUnitX->strLexUnit,lexUnitX));
		++m_iLexicalUnitsSentenceDelimiter;
		
		// (2) insert lexical units read from the file
		char strLexUnit[MAX_LEXUNIT_LENGTH+1];
		int iLine = 0;
		string strLine;
		while(std::getline(file.getStream(),strLine)) {
			++iLine;
			// skip comments and blank lines
			if (strLine.empty() || ((strLine.c_str()[0] == '#') && (strLine.c_str()[1] == '#')) ||
				((strLine.length() >= 3) && (strLine.c_str()[0] == ';') && (strLine.c_str()[1] == ';') && (strLine.c_str()[2] == ';'))) {
				continue;
			}
			
			LexUnit *lexUnit = processLine(strLine.c_str(),iLine,strLexUnit);
			assert(lexUnit);
			
			lexUnit->fProbability = 1.0; 
			// insert it into the lexicon map
			MLexUnit::iterator it = m_mLexUnit->find(strLexUnit);
			if (it == m_mLexUnit->end()) {
				// insert the lexical unit into the lexicon
				LexUnitX *lexUnitX = new LexUnitX;
				lexUnitX->iLexUnit = iLexUnitId++;
				lexUnitX->vLexUnitPronunciations.push_back(lexUnit);
				char *strLexUnitFit = new char[strlen(strLexUnit)+1];
				strcpy(strLexUnitFit,strLexUnit);
				pair<MLexUnit::iterator, bool> pair = m_mLexUnit->insert(MLexUnit::value_type(strLexUnitFit,lexUnitX));
				lexUnitX->strLexUnit = strLexUnitFit;
				lexUnit->iLexUnit = lexUnitX->iLexUnit;
				m_lexiconX.push_back(lexUnitX);
			} else {
				// insert the lexical unit into the lexicon 
				lexUnit->iLexUnit = it->second->iLexUnit;
				it->second->vLexUnitPronunciations.push_back(lexUnit);
			}
			// insert it into the lexicon
			m_lexicon.push_back(lexUnit);
			attachIndex(lexUnit);
		}

		file.close();
		
		// attach an index to those units that are not in the lexicon file
		attachIndex(m_lexUnitUnknown);
		attachIndex(m_lexUnitEndSentence);
		attachIndex(m_lexUnitBegSentence); 
		
		// create lexical unit ids including pronunciations
		int iLexUnitPron = 0;
		for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
			for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
				(*jt)->iLexUnitPron = iLexUnitPron;
				++iLexUnitPron;
			}
		} 
		assert((int)m_lexicon.size() == iLexUnitPron);
		
		// sanity check
		for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
			int j = 0;
			for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt, ++j) {
				assert((*jt)->iPronunciation == j);
			}
		}
   
   } catch(ExceptionBase) {
   	BVC_ERROR << "unable to load the lexicon file: " << m_strFile;
   }
}

// process a line in the lexicon and extract the lexical unit and its phonetic transcription
LexUnit *LexiconManager::processLine(const char *strLine, int iLine, char *strLexUnit) {
	
	int i = 0;
	int iCharacters = 0;
	int iTokens = 0;
	char strBuffer[MAX_LEXUNIT_LENGTH+1];
	unsigned char iType = UCHAR_MAX;
	float fPronunciationProbability = 1.0;
	
	LexUnit *lexUnit = new LexUnit;
	
	// process characters until the end of line
	while((strLine[i] != '\0') && (strLine[i] != '\r') && (strLine[i] != '\n')) {
		// if a separator is found then 
		if ((strLine[i] == ' ') || (strLine[i] == '\t')) {
			if (iCharacters > 0) {
				// the first token is the lexical unit
				if (iTokens == 0) {
					++iTokens;
					strBuffer[iCharacters] = '\0';
					// keep the lexical unit
					strcpy(strLexUnit,strBuffer);
					iCharacters = 0;
				} 
				// the rest of the tokens are the pronunciation probability (optional) or part of the phonetic transcription
				else {
					// check if it is the pronunciation probability
					if ((iTokens == 1) && (iCharacters == 6) && ((strBuffer[0] == '0') || (strBuffer[0] == '1')) 
						&& (strBuffer[1] == '.')) {
						++iTokens;
						strBuffer[iCharacters] = '\0';
						// check correctness
						for(int j=0 ; j < iCharacters ; ++j) {
							// skip the '.'
							if (j == 1) {
								continue;
							}
							if (isdigit(strBuffer[j]) == 0) {
								BVC_ERROR << "pronunciation probability " << strBuffer 
									<< " is not in the right format: line " << iLine;
							}
						}
						// convert the string to a probability
						fPronunciationProbability = (float)atof(strBuffer+1);
						if ((fPronunciationProbability > 1.0) || (fPronunciationProbability < 0.0)) {
							BVC_ERROR << "non probabilistic value found \"" << strBuffer+1 << "\": line " << iLine;
						}
						iCharacters = 0;
					}
					// read the phone
					else {
						++iTokens;
						strBuffer[iCharacters] = '\0';
						//get the phone index
						int iPhoneIndex = m_phoneSet->getPhoneIndex(strBuffer);
						if (iPhoneIndex == -1) {
							BVC_ERROR << "unknown phonetic symbol: \"" << strBuffer << 
								"\", it is not defined in the phonetic symbol set: line " << iLine;
						}
						lexUnit->vPhones.push_back(iPhoneIndex);
						iCharacters = 0;
					}
				}
			}
		} else {
			if (iTokens == 0) {
				if (iCharacters+1 >= MAX_LEXUNIT_LENGTH) {
					BVC_ERROR << "lexical unit defined in line " << iLine << " of the lexicon exceeds the maximum length (" << 
					MAX_LEXUNIT_LENGTH << " characters)";
				}
			} else {
				if (iCharacters+1 >= MAX_PHONETIC_SYMBOL_LENGTH) {
					BVC_ERROR << "phonetic symbol defined in line " << iLine << " of the lexicon exceeds the maximum length (" << 
					MAX_PHONETIC_SYMBOL_LENGTH << " characters)";
				}	
			}
			strBuffer[iCharacters++] = strLine[i];
		}
		++i;
	}
	if (iCharacters > 0) {
		++iTokens;
		strBuffer[iCharacters] = '\0';
		//get the phone index
		int iPhoneIndex = m_phoneSet->getPhoneIndex(strBuffer);
		if (iPhoneIndex == -1) {
			BVC_ERROR << "unknown phonetic symbol: \"" << strBuffer 
				<< "\", it is not defined in the phonetic symbol set: line " << iLine;
		}
		lexUnit->vPhones.push_back(iPhoneIndex);
		iCharacters = 0;
	}
	
	// check completeness
	if (iTokens == 0) {
		BVC_ERROR << "no lexical unit defined at line " << iLine << " of the lexicon";
	} else if (iTokens < 2) {
		BVC_ERROR << "no phonetic symbol defined for lexical unit " << strLexUnit << 
			" at line " << iLine << " of the lexicon";
	}
	
	// extract the pronunciation number from the lexical unit
	int iPronunciation = 0;
	char *cAlternative = strrchr(strLexUnit,'(');
	if (strLexUnit[strlen(strLexUnit)-1] != ')')
		cAlternative = NULL;
	if (cAlternative) {
		for(int i = 1 ; cAlternative[i] != ')' ; ++i) {
			// make sure it is a number (by the way this checks that the end of string is not reached before the ')')
			if ((cAlternative[i] > 57) || (cAlternative[i] < 48)) {
				BVC_ERROR << "lexicon file is not in a valid format: line " + iLine;
			}
			iPronunciation = iPronunciation*10;
			iPronunciation += cAlternative[i]-48;	
		}			
		iPronunciation--;
		*cAlternative = '\0'; 
		
		// get the lexical unit type and check its correctness
		iType = getLexicalUnitType(strLexUnit);
		if ((iType != LEX_UNIT_TYPE_STANDARD) && (iType != LEX_UNIT_TYPE_FILLER)) {
			BVC_ERROR << "lexical unit \"" << strLexUnit << "\" at line " << iLine << " cannot be included in the lexicon";
		}
		
		// make sure alternative pronunciations of this lexical unit were already seen
		bool bCorrectId = true;
		int iSize = 0;
		MLexUnit::iterator it = m_mLexUnit->find(strLexUnit);
		if (it != m_mLexUnit->end()) {
			iSize = (int)it->second->vLexUnitPronunciations.size();
			if (iSize != iPronunciation) {
				bCorrectId = false;
			}
		} else {
			bCorrectId = false;
		}
		if (bCorrectId == false) {
			int iPronunciationOld = iPronunciation+1;			
			iPronunciation = iSize;
			BVC_WARNING << "lexical unit \"" << strLexUnit << "(" << iPronunciationOld << ")\" at line " << iLine
				 << " has an incorrect pronunciation id, renamed to: \"" << strLexUnit << "(" << iPronunciation+1 << ")\"";
		}
	} else {
		// get the lexical unit type and check its correctness
		iType = getLexicalUnitType(strLexUnit);
		if ((iType != LEX_UNIT_TYPE_STANDARD) && (iType != LEX_UNIT_TYPE_FILLER)) {
			BVC_ERROR << "lexical unit \"" << strLexUnit << "\" at line " << iLine << 
				" cannot be included in the lexicon";
		}
	
		// in case the lexical unit does not have an alternative pronunciation id and was already seen in the lexicon, create a pronunciation id and report a warning	
		MLexUnit::iterator it = m_mLexUnit->find(strLexUnit);
		if (it != m_mLexUnit->end()) {
			// get the next available pronunciation id
			iPronunciation = (int)it->second->vLexUnitPronunciations.size();
			BVC_WARNING << "lexical unit \"" << strLexUnit << "\" at line " << iLine 
				<< " was previously defined in the lexicon, renamed to: \"" 
				<< strLexUnit << "(" << iPronunciation+1 << ")\"";
		}
	}
	lexUnit->iPronunciation = iPronunciation;
	lexUnit->fProbability = fPronunciationProbability;
	lexUnit->iType = iType;
	lexUnit->fInsertionPenalty = 0.0;
	
	// accumulate type statistics
	switch(lexUnit->iType) {
		case LEX_UNIT_TYPE_STANDARD: {
		 	++m_iLexicalUnitsStandard;
		 	break;
		}
		case LEX_UNIT_TYPE_FILLER: {
		 	++m_iLexicalUnitsFiller;
		 	break;
		}
		case LEX_UNIT_TYPE_SENTENCE_DELIMITER: {
		 	++m_iLexicalUnitsSentenceDelimiter;
		 	break;
		}
		case LEX_UNIT_TYPE_UNKNOWN: {
		 	++m_iLexicalUnitsUnknown;
		 	break;
		}
		default: {
			assert(0);
		}
	}
	
	// accumulate pronunciation statistics
	if (lexUnit->iPronunciation == 0) {
		m_iLexicalUnitsUnique++;
	}
	
	return lexUnit;
}

// attach insertion-penalties to lexical units in the lexicon (needed for decoding)
void LexiconManager::attachLexUnitPenalties(float fInsertionPenaltyStandard, float fInsertionPenaltyFiller) {

	for(VLexUnit::iterator it = m_lexicon.begin() ; it != m_lexicon.end() ; ++it) {
		if ((*it)->iType == LEX_UNIT_TYPE_STANDARD) {
			(*it)->fInsertionPenalty = fInsertionPenaltyStandard;
		} 
		else if ((*it)->iType == LEX_UNIT_TYPE_FILLER) {
			(*it)->fInsertionPenalty = fInsertionPenaltyFiller;
		}
		else {
			assert(((*it)->iType == LEX_UNIT_TYPE_SENTENCE_DELIMITER) || ((*it)->iType == LEX_UNIT_TYPE_UNKNOWN));
		}
	}	
}

// return the lexical-unit type
unsigned char LexiconManager::getLexicalUnitType(const char *strLexUnit) {

	unsigned int iLength = (unsigned int)strlen(strLexUnit);
	assert(iLength > 0);

	// beginning/end of sentence
	if ((strcmp(strLexUnit,LEX_UNIT_BEGINNING_SENTENCE) == 0) || (strcmp(strLexUnit,LEX_UNIT_END_SENTENCE) == 0)) {
		return LEX_UNIT_TYPE_SENTENCE_DELIMITER;
	}	
	// unknown
	else if (strcmp(strLexUnit,LEX_UNIT_UNKNOWN) == 0) {
		return LEX_UNIT_TYPE_UNKNOWN;
	}
	// filler
	else if ((strLexUnit[0] == '<') && (strLexUnit[iLength-1] == '>')) {
		return LEX_UNIT_TYPE_FILLER;	
	}
	// standard
	else {
		return LEX_UNIT_TYPE_STANDARD;
	}	
}

// clean-up
void LexiconManager::destroy() {
   
   for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
   	if ((strcmp((*it)->strLexUnit,LEX_UNIT_UNKNOWN) != 0) &&
   		(strcmp((*it)->strLexUnit,LEX_UNIT_BEGINNING_SENTENCE) != 0) && 
   		(strcmp((*it)->strLexUnit,LEX_UNIT_END_SENTENCE) != 0)) {
			delete [] (*it)->strLexUnit;
		}
		for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
			delete *jt;
		}
      delete *it;
   }
   
   m_lexicon.clear();
   m_lexiconX.clear();
   
   if (m_mLexUnit != NULL) {
		/*for(MLexUnit::iterator it = m_mLexUnit->begin() ; it != m_mLexUnit->end() ; ++it) {
			delete [] it->first;
		}*/
	   m_mLexUnit->clear();
   	delete m_mLexUnit;
   }
}

// prints the lexicon
void LexiconManager::print(bool bPrintLexUnits) {

	printf("VLexUnit:\n");
	printf(" ## lexical units (standard): %8d \n",m_iLexicalUnitsStandard);
	printf(" ## lexical units (filler):   %8d\n",m_iLexicalUnitsFiller); 
	printf(" alt. pronunciation ratio:   %8.2f\n",((float)m_iLexicalUnitsStandard)/((float)m_iLexicalUnitsUnique));
	if (bPrintLexUnits) {
		for(VLexUnit::iterator it = m_lexicon.begin() ; it != m_lexicon.end() ; ++it) {
			const char *str = LexiconManager::getStrLexUnit((*it)->iLexUnit);	
			printf("%10d[%10d]: %s(%d)",(*it)->iLexUnit,(*it)->iLexUnitPron,str,(*it)->iPronunciation);
			for(vector<int>::iterator jt = (*it)->vPhones.begin() ; jt != (*it)->vPhones.end() ; ++jt) {
				printf(" %s",m_phoneSet->getStrPhone(*jt));
			}
			printf("\n");
		}
	}
}


// return a vector containing all the filler lexical units (special-optional lexical units)
void LexiconManager::getVLexUnitFiller(VLexUnit &vLexUnit) {

	// iterate through all the lexical units
   for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
   	if (isFiller((*it)->iLexUnit)) {
   		assert((*it)->vLexUnitPronunciations.size() == 1);
			vLexUnit.push_back((*it)->vLexUnitPronunciations.front());
   	}
   }	
}

// compute the average number of alternative pronunciations in the lexicon
float LexiconManager::computeAveragePronunciationVariations() {

	return ((float)m_iLexicalUnitsStandard)/((float)m_iLexicalUnitsUnique);
}

// creates a new lexicon file from seen lexical units (the new lexicon is a subset of the original lexicon)
// all the lexical units have to correspond to the actual lexicon object
bool LexiconManager::createLexicon(const char *strFile, MLexUnitInt &mLexUnitSeen, MLexUnitLexUnit &mLexUnitLexUnit) {

	// open the file for reading
	FILE *file = fopen(strFile,"wt");
	if (file == NULL) {
		printf("Error: unable to create the lexicon file: %s\n",strFile);
		return false;
	}
	
	// (1) compute pronunciation probabilities from lexical unit frequencies
	// note that lexical units are sorted alphabetically and by pronunciation	
	LLexUnit lLexUnit;	// list of lexical units as they will be written to the file
	int iIndex = 0;	
	int iLexUnitPrev = -1;
	int iOccurrencesTotal = 0;			// within-lexical unit occurrences (all alternative pronunciations)	
	for(MLexUnitInt::iterator it = mLexUnitSeen.begin() ; it != mLexUnitSeen.end() ; ++it) {
	
		// create a copy of the lexical unit
		LexUnit *lexUnitAux = new LexUnit;
		lexUnitAux->iLexUnit = it->first->iLexUnit;
		for(vector<int>::iterator kt = it->first->vPhones.begin() ; kt != it->first->vPhones.end() ; ++kt) {
			lexUnitAux->vPhones.push_back(*kt);
		}
		lexUnitAux->iPronunciation = UCHAR_MAX;									// yet to be determined
		lexUnitAux->iType = it->first->iType;
		lexUnitAux->fProbability = -1;												// yet to be determined
		lexUnitAux->fInsertionPenalty = it->first->fInsertionPenalty;
		
		// compute the probability
		// update the denominator?
		if (it->first->iLexUnit != iLexUnitPrev) {
			iLexUnitPrev = it->first->iLexUnit;
			iOccurrencesTotal = 0;			
			for(MLexUnitInt::iterator jt = it ; ((jt != mLexUnitSeen.end()) && (jt->first->iLexUnit == it->first->iLexUnit)) ; ++jt) {
				assert(it->second > 0.0);
				iOccurrencesTotal += jt->second;
				++iIndex;
			}
		}
		lexUnitAux->fProbability = ((float)it->second)/((float)iOccurrencesTotal);
		
		// create a new entry in the map
		mLexUnitLexUnit.insert(MLexUnitLexUnit::value_type(it->first,lexUnitAux));
		lLexUnit.push_back(lexUnitAux);
	}
	
	assert(iIndex == (int)mLexUnitSeen.size());
	
	// sort the list of lexical units alphabetically and then by pronunciation probability
	lLexUnit.sort(LexiconManager::comparePronProbability);	
	
	// get the lenght of the longest lexical unit (for formatting purposes)
	int iLengthMax = -1;
	for(MLexUnitInt::iterator it = mLexUnitSeen.begin() ; it != mLexUnitSeen.end() ; ++it) {		
		int iLen = (int)strlen(getStrLexUnit(it->first->iLexUnit));
		if (it->first->iPronunciation > 0) {
			iLen += 3;
		}
		if (iLen > iLengthMax) {
			iLengthMax = iLen;
		}
	}
	
	// (2) write the lexical units to the file
	// note: pronunciation numbers might need to be recalculated because some of the original pronunciations might not be seen
	char strAux[MAX_LEXUNIT_LENGTH+1];
	iLexUnitPrev = -1;
	int iPronunciation = -1;
	int i=0;
	for(LLexUnit::iterator it = lLexUnit.begin() ; it != lLexUnit.end() ; ++it, ++i) {
		// determine the pronunciation number
		if ((*it)->iLexUnit == iLexUnitPrev) {
			iPronunciation++;
		} else {
			iLexUnitPrev = (*it)->iLexUnit;
			iPronunciation = 0;
		}
		(*it)->iPronunciation = iPronunciation;
		// print the lexical unit
		if (iPronunciation == 0) { 
			sprintf(strAux,"%s",getStrLexUnit((*it)->iLexUnit));
		} else if (iPronunciation >= 1) {
			sprintf(strAux,"%s(%d)",getStrLexUnit((*it)->iLexUnit),(*it)->iPronunciation+1);
		} else {
			assert(0);
		}
		fprintf(file,"%-*s",iLengthMax,strAux);	
		
		// print the probability
		assert(((*it)->fProbability >= 0.0) && ((*it)->fProbability <= 1.0));
		if ((*it)->fProbability == 1.0) {
			fprintf(file,"           ");
		} else {
			fprintf(file,"   %.4f  ",(*it)->fProbability);
		}
		
		// print the phonetic transcription
		for(vector<int>::iterator jt = (*it)->vPhones.begin() ; jt != (*it)->vPhones.end() ; ++jt) {
			fprintf(file," %s",m_phoneSet->getStrPhone(*jt));
		}
		fprintf(file,"\n");
	}
	
	// close the file
	if (fclose(file) == EOF) {
		return false;
	}

	return true;
}

// return the silence lexical unit
LexUnit *LexiconManager::getLexUnitSilence() {

	// find the silence lexical unit
	LexUnitX *lexUnitXSilence = getLexUnit(LEX_UNIT_SILENCE_SYMBOL);
	assert(lexUnitXSilence->vLexUnitPronunciations.size() == 1);
	
	return lexUnitXSilence->vLexUnitPronunciations.front();
}

// tranform a sequence of lexical units from text format to object format
bool LexiconManager::getLexUnits(const char *strText, VLexUnit &vLexUnit, bool &bAllKnown) {

	bAllKnown = true;
	// process the text character by character
	int iLength = (int)strlen(strText);
	char strLexUnit[MAX_LEXUNIT_LENGTH+1];
	int iCharacters = 0;
	for(int i=0 ; i<iLength ; ++i) {
		if (isspace(strText[i]) != 0) {
			if (iCharacters > 0) {
				strLexUnit[iCharacters] = 0;
				LexUnit *lexUnit = getLexUnitPronunciation(strLexUnit);
				if (lexUnit == NULL) {
					// the lexical unit is not in the lexicon, replace it by the <UNK> lexical unit
					lexUnit = m_lexUnitUnknown;
					bAllKnown = false;
				}
				vLexUnit.push_back(lexUnit);
				iCharacters = 0;	
			}
		} else {
			// check if the lexical units exceeds the maximum length
			if (iCharacters >= MAX_LEXUNIT_LENGTH) {
				return false;
			}
			strLexUnit[iCharacters] = strText[i];
			++iCharacters;	
		}
	}
	if (iCharacters > 0) {
		strLexUnit[iCharacters] = 0;
		LexUnit *lexUnit = getLexUnitPronunciation(strLexUnit);
		if (lexUnit == NULL) {
			// the lexical unit is not in the lexicon, replace it by the <UNK> lexical unit
			lexUnit = m_lexUnitUnknown;
		}
		vLexUnit.push_back(lexUnit);
	}

	return true;
}

// tranform a sequence of lexical units from text format to object format
bool LexiconManager::getLexUnits(const char *strText, VLexUnitX &vLexUnitX, bool &bAllKnown) {

	bAllKnown = true;
	// process the text character by character
	int iLength = (int)strlen(strText);
	char strLexUnit[MAX_LEXUNIT_LENGTH+1];
	int iCharacters = 0;
	for(int i=0 ; i<iLength ; ++i) {
		if (isspace(strText[i]) != 0) {
			if (iCharacters > 0) {
				strLexUnit[iCharacters] = 0;
				LexUnitX *lexUnitX = getLexUnit(strLexUnit);
				if (lexUnitX == NULL) {					
					// the lexical unit is not in the lexicon, replace it by the <UNK> lexical unit
					lexUnitX = getLexUnit(m_lexUnitUnknown->iLexUnit);
					bAllKnown = false;
				}
				vLexUnitX.push_back(lexUnitX);
				iCharacters = 0;	
			}
		} else {
			// check if the lexical units exceed the maximum length
			if (iCharacters >= MAX_LEXUNIT_LENGTH) {
				return false;
			}
			strLexUnit[iCharacters] = strText[i];
			++iCharacters;	
		}
	}
	if (iCharacters > 0) {
		strLexUnit[iCharacters] = 0;
		LexUnitX *lexUnitX = getLexUnit(strLexUnit);
		if (lexUnitX == NULL) {
			// the lexical unit is not in the lexicon, replace it by the <UNK> lexical unit
			lexUnitX = getLexUnit(m_lexUnitUnknown->iLexUnit);
		}
		vLexUnitX.push_back(lexUnitX);
	}

	return true;
}


// check that the lexical units are consistent with language model (ith lexical unit == ith unigram)
bool LexiconManager::checkConsistency(vector<string> &vStrLexUnit) {

	int i=0;
	for(vector<string>::iterator it = vStrLexUnit.begin() ; it != vStrLexUnit.end() ; ++it, ++i) {
		//printf("%s %s\n",(*it).c_str(),m_lexiconX[i]->strLexUnit);
		if ((*it).compare(m_lexiconX[i]->strLexUnit) != 0) {
			return false;
		}
	}

	return true;
}

// arrange the lexical units to match the order in the given vector of unigrams
void LexiconManager::arrangeLexUnits(vector<string> &vStrLexUnit) { 

	// initialize all the ids to -1 and move them to a temporal vector
	for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
		(*it)->iLexUnit = -1;
	}
	
	// create a vector with the lexical units ordered as in the vector of unigrams
	int iIndex = 0;
	VLexUnitX vLexUnitX;
	for(vector<string>::iterator it = vStrLexUnit.begin() ; it != vStrLexUnit.end() ; ++it, ++iIndex) {
		LexUnitX *lexUnitX = getLexUnit((*it).c_str());
		if (lexUnitX == NULL) {
			BVC_ERROR << "lexical unit \"" << (*it).c_str() << "\" is not in the lexicon (is it in the language model?)";
		}
		lexUnitX->iLexUnit = iIndex;
		for(VLexUnit::iterator jt = lexUnitX->vLexUnitPronunciations.begin() ; jt != lexUnitX->vLexUnitPronunciations.end() ; ++jt) {
			(*jt)->iLexUnit = iIndex;
		}
		vLexUnitX.push_back(lexUnitX);
	}
	
	// delete lexical units that are not in the language model (except the unknown lexical unit,
	// the filler lexical units and the sentence delimiters, which are used by the decoder)
	for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
		if ((*it)->iLexUnit == -1) {
			// IMPORTANT: sentence delimiters and unknown need to be preserved (are needed by the decoder)
			if ((strcmp((*it)->strLexUnit,LEX_UNIT_UNKNOWN) == 0) || (strcmp((*it)->strLexUnit,LEX_UNIT_BEGINNING_SENTENCE) == 0) || (strcmp((*it)->strLexUnit,LEX_UNIT_END_SENTENCE) == 0)) {
				(*it)->iLexUnit = (int)vLexUnitX.size();
				assert((*it)->vLexUnitPronunciations.size() == 1);
				(*it)->vLexUnitPronunciations.front()->iLexUnit = (*it)->iLexUnit;
				vLexUnitX.push_back(*it);
			} 
			// fillers: such as <UHM>, <BR>, etc
			else if ((*it)->vLexUnitPronunciations.front()->iType == LEX_UNIT_TYPE_FILLER) {
				(*it)->iLexUnit = (int)vLexUnitX.size();
				for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
					(*jt)->iLexUnit = (*it)->iLexUnit;
				}
				vLexUnitX.push_back(*it);
			}
			// these are words that do not have ids
			else {
				printf("lexical unit \"%s\" was removed from the lexicon\n",(*it)->strLexUnit);
				for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
					delete *jt;
				}
				delete *it;
			}
		}	
	}
	
	// clear the data structures that keep the old arrangement of lexical units
	m_lexiconX.clear();
	m_lexicon.clear();
	m_mLexUnit->clear();
		
	// regenerate the data structures and create lexical unit ids including pronunciations
   int iLexUnitPron = 0;
	for(VLexUnitX::iterator it = vLexUnitX.begin() ; it != vLexUnitX.end() ; ++it) {
		for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
   		(*jt)->iLexUnitPron = iLexUnitPron++;
			m_lexicon.push_back(*jt);
		}
		m_lexiconX.push_back(*it);
		m_mLexUnit->insert(MLexUnit::value_type((*it)->strLexUnit,*it));
	}	
   assert((int)m_lexicon.size() == iLexUnitPron);
}

// store the lexicon into a file
bool LexiconManager::store(const char *strFile) {

	char strLexUnitPronunciation[1024];

	// create the output file
	FILE *file = fopen(strFile,"wt");
	if (file == NULL) {
		return false;
	}
	
	int iTotalLexicalUnits = m_iLexicalUnitsFiller+m_iLexicalUnitsStandard;
	
	fprintf(file,"## ---------------------------------------------------\n");
	fprintf(file,"## total lexical units:  %7d\n",iTotalLexicalUnits);
	fprintf(file,"## - fillers:            %7d\n",m_iLexicalUnitsFiller);
	fprintf(file,"## - standard(unique):   %7d\n",m_iLexicalUnitsUnique);
	fprintf(file,"## - standard(+prons):   %7d\n",m_iLexicalUnitsStandard);
	fprintf(file,"## pronunciation-ratio:  %7.2f\n",((float)m_iLexicalUnitsStandard)/((float)m_iLexicalUnitsUnique));
	fprintf(file,"## ---------------------------------------------------\n");
	
   for(VLexUnitX::iterator it = m_lexiconX.begin() ; it != m_lexiconX.end() ; ++it) {
   	for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {
			if ((isSentenceDelimiter(*jt)) || (isUnknown(*jt))) {
				continue;
			}
   		getStrLexUnitPronunciation(*jt,strLexUnitPronunciation);
   		fprintf(file,"%-30s",strLexUnitPronunciation);
   		for(vector<int>::iterator kt = (*jt)->vPhones.begin() ; kt != (*jt)->vPhones.end() ; ++kt) {
				fprintf(file," %s",m_phoneSet->getStrPhone(*kt));	
   		}
   		fprintf(file,"\n");
   	}
   } 
	
	// close the file
	if (fclose(file) == EOF) {	
		return false;
	}	
	
	return true;
}


// map lexical units from pronunciation format to non-pron format
void LexiconManager::map(VLexUnit &vLexUnitInput, VLexUnitX &vLexUnitOutput) {

	for(VLexUnit::iterator it = vLexUnitInput.begin() ; it != vLexUnitInput.end() ; ++it) {
		vLexUnitOutput.push_back(getLexUnit((*it)->iLexUnit));
	}
}

// remove non-standard lexical units
void LexiconManager::removeNonStandardLexUnits(VLexUnit &vLexUnit) {

	for(VLexUnit::iterator it = vLexUnit.begin() ; it != vLexUnit.end() ; ) {
		if (isStandard(*it) == false) {
			it = vLexUnit.erase(it);
		} else {
			++it;
		}
	}
}

};	// end-of-namespace

