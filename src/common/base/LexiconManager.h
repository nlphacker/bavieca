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


#ifndef LEXICONMANAGER_H
#define LEXICONMANAGER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

using namespace std;

#include <string>
#include <vector>
#include <list>
#include <map>
#include <sstream>

#include "Global.h"
#include "PhoneSet.h"

#if defined __linux__ || defined __APPLE__ || __MINGW32__
#include <tr1/unordered_map>
#elif _MSC_VER
#include <hash_map>
#else 
	#error "unsupported platform"
#endif

namespace Bavieca {

#define NON_LEXUNIT_ID      -100

// maximum length of a lexical unit including the end-of-string character 
#define MAX_LEXUNIT_LENGTH  1024 

// maximum number of alternative pronunciations a lexical unit can have
#define MAX_LEXUNIT_PRONUNCIATIONS	1024

// lexical unit types
#define LEX_UNIT_TYPE_STANDARD 				0
#define LEX_UNIT_TYPE_FILLER					1
#define LEX_UNIT_TYPE_SENTENCE_DELIMITER	2
#define LEX_UNIT_TYPE_UNKNOWN					3

#define LEX_UNIT_SILENCE_SYMBOL				"<SIL>"
#define LEX_UNIT_BEG_SENTENCE					"<s>"
#define LEX_UNIT_END_SENTENCE					"</s>"
#define LEX_UNIT_UNKNOWN						"<UNK>"

// value for uninitialized variables that keep lexical unit indices
#define UNINITIALIZED_LEX_UNIT   LONG_MIN

// lexical unit transcription
typedef struct {
   int iLexUnit;							// lexical unit id (unique correspondence id <> lexical unit in str format)
   int iLexUnitPron;						// lexical unit id (unique correspondence id <> lexical unit + pronunciation)
   int iIndex;								// index of the lexical unit within the lexicon file
   vector<int> vPhones;					// phonetic transciption 
   unsigned char iPronunciation;		// pronunciation number
   unsigned char iType;					// standard / filler / sentence delimiter
   float fProbability;					// pronunciation probability (respect to alternative ones of the same lex unit)
   float fInsertionPenalty;			// penalty for inserting this lexical unit during decoding
} LexUnit;	

class LexiconManager;

typedef vector<LexUnit*> VLexUnit;
typedef list<LexUnit*> LLexUnit;


// ad-hoc functions to use Duple as the key in a hash_map data structure
struct MLexUnitFunctions
{

#if defined __linux__ || defined __APPLE__ || __MINGW32__

	// comparison function (used for matching, comparison for equality)
	bool operator()(const char *strLexUnit1, const char *strLexUnit2) const {
		
		return (strcmp(strLexUnit1,strLexUnit2) == 0);
	}

#elif _MSC_VER

	static const size_t bucket_size = 4;
	static const size_t min_buckets = 8;

	// comparison function (used to order elements)
	bool operator()(const char *strLexUnit1, const char *strLexUnit2) const {	
		
		return (strcmp(strLexUnit1,strLexUnit2) < 0);
	}

#endif
	
	// hash function
	size_t operator()(const char *strLexUnit) const {
	
		unsigned int iAcc = 0;
		unsigned int iAux = 0;
		for(int i=0 ; strLexUnit[i] != 0 ; ++i) {
			if (i <= 3) {	
				iAcc <<= (8*i);
				iAcc += strLexUnit[i];
			} else {
				iAux = strLexUnit[i];
				iAux <<= (8*(i%4));
				iAcc ^= iAux;
			}
		}	
	
		return iAcc;
	}
};

typedef map<LexUnit*,int,bool(*)(const LexUnit*, const LexUnit*)> MLexUnitInt;

typedef struct {
	int iLexUnit;									// lexical unit unique identifier
	const char *strLexUnit;						// lexical unit as an array of characters (a word, a syllable, etc.)
	VLexUnit vLexUnitPronunciations;			// alternative pronunciations of the lexical unit
} LexUnitX;

typedef vector<LexUnitX*> VLexUnitX;

// maps lexical units as strings of characters to their corresponding data structure
#if defined __linux__ || defined __APPLE__ || __MINGW32__
typedef std::tr1::unordered_map<const char*,LexUnitX*,MLexUnitFunctions,MLexUnitFunctions> MLexUnit;
#elif _MSC_VER
typedef hash_map<const char*,LexUnitX*,MLexUnitFunctions> MLexUnit;
#else
	#error "unsupported platform"
#endif

// maps lexical units within one lexicon to lexical units within another lexicon
typedef map<LexUnit*,LexUnit*> MLexUnitLexUnit;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class LexiconManager {
   
   private:
      
      string m_strFile;											// lexicon file name
      VLexUnit m_lexicon;										// lexicon
      VLexUnitX m_lexiconX;									// lexicon
      MLexUnit	*m_mLexUnit;									// lexical unit as a string <> lexical unit information
      PhoneSet *m_phoneSet;									// phonetic symbol set
		int m_iLexicalUnitsStandard;		
		int m_iLexicalUnitsFiller;		
		int m_iLexicalUnitsSentenceDelimiter;		
		int m_iLexicalUnitsUnknown; 
		int m_iLexicalUnitsUnique;								// number of unique lexical units (regardless of alternative pronunciations)
		
		VLexUnit m_vLexUnitLexiconFile;						// lexical units as they were read from the file
		
		// attach an index to the lexical unit
		inline void attachIndex(LexUnit *lexUnit) {
		
			lexUnit->iIndex = (int)m_vLexUnitLexiconFile.size();
			m_vLexUnitLexiconFile.push_back(lexUnit);
		}
   
		// process a line in the lexicon and extract the lexical unit and its phonetic transcription
		LexUnit *processLine(const char *strLine, int iLine, char *strLexUnit);
		
		// return whether the lexical unit is redundant
		bool isRedundant(LexUnit *lexUnit, const char *strLexUnit);	
   
   public:
   
		// these lexical units are allways defined (needed by the decoder)
		LexUnit *m_lexUnitBegSentence;
		LexUnit *m_lexUnitEndSentence;
		LexUnit *m_lexUnitUnknown;
   
   	//note: the silence and the filler may or may not be a lexical unit, only if they are included in the lexicon
      
      // constructor
      LexiconManager(const char *strFile, PhoneSet *phoneSet);

      // destructor
      ~LexiconManager();
    
      // load the lexicon from a file
      void load();
      
      // clean-up
      void destroy();
      
		// prints the lexicon information and optionally the lexical units
		void print(bool bPrintLexUnits = false);
		
		// attach insertion-penalties to lexical units in the lexicon (needed for decoding)
		void attachLexUnitPenalties(float fInsertionPenaltyStandard, float fInsertionPenaltyFiller);	
		
		// return a lexical unit by its index in the lexicon file
		inline LexUnit *getLexUnitByIndex(unsigned int iIndex) {
		
			if (iIndex >= m_vLexUnitLexiconFile.size()) {
				return NULL;
			}
			
			return m_vLexUnitLexiconFile[iIndex];
		}
      
      // return a reference to the lexicon
      inline VLexUnit *getLexiconReference() {
      
      	return &m_lexicon;
      }
      
      // return a reference to the lexicon
      inline VLexUnitX *getLexiconXReference() {
      
      	return &m_lexiconX;
      }
		
		// return the vocabulary size (excluding alternative pronunciations)
		inline unsigned int getVocabularySize() {
		
			return (unsigned int)m_lexiconX.size();
		}
		
		// return the size of the lexicon (vocabulary size including alternative pronunciations)
		inline unsigned int getLexiconSize() {

			return (unsigned int)m_lexicon.size();
		}
		
		// return the lexical unit identifier attached to a lexical unit
		inline int getLexUnitId(const string &strLexUnit) {
		
			return getLexUnitId(strLexUnit.c_str());
		}
		
		// return the lexical unit identifier attached to a lexical unit
		inline int getLexUnitId(const char *strLexUnit) {
		
			MLexUnit::iterator it = m_mLexUnit->find(strLexUnit);
			if (it == m_mLexUnit->end()) {
				return -1;
			}
			
			return it->second->iLexUnit;
		}
		
		// return a lexical unit 
		inline LexUnitX *getLexUnit(const string &strLexUnit) {
		
			MLexUnit::iterator it = m_mLexUnit->find(strLexUnit.c_str());
			if (it == m_mLexUnit->end()) {
				return NULL;
			}
			
			return it->second;
		}
		
		// return the silence lexical unit
		LexUnit *getLexUnitSilence();	
		
		// return a lexical unit given the base and the pronunciation number 
		inline LexUnit *getLexUnit(unsigned int iLexUnit, unsigned int iPronunciation) {
		
			// check
			if (iLexUnit >= m_lexiconX.size()) {
				return NULL;
			}
			
			// check
			if (iPronunciation >= m_lexiconX[iLexUnit]->vLexUnitPronunciations.size()) {
				return NULL;
			}
		
			return m_lexiconX[iLexUnit]->vLexUnitPronunciations[iPronunciation];
		}
		
		// return a lexical unit given the base and the pronunciation number 
		inline LexUnit *getLexUnit(const char *strLexUnit, unsigned int iPronunciation) {
		
			MLexUnit::iterator it = m_mLexUnit->find(strLexUnit);
			if (it == m_mLexUnit->end()) {
				return NULL;
			}
			
			// get the right pronunciation
			if (it->second->vLexUnitPronunciations.size() < iPronunciation+1) {
				return NULL;
			}
			
			return it->second->vLexUnitPronunciations[iPronunciation];
		}
		
		// return a lexical unit 
		inline LexUnitX *getLexUnit(const char *strLexUnit) {
		
			string strAux = strLexUnit;
		
			return getLexUnit(strAux);
		}
		
		// return a lexical unit 
		inline LexUnit *getLexUnitPronunciation(const char *strLexUnitPronunciation) {
		
			char strLexUnit[MAX_LEXUNIT_LENGTH];
			unsigned char iPronunciation;
			
			if (strcmp(strLexUnitPronunciation,LEX_UNIT_BEG_SENTENCE) == 0) {
				return m_lexUnitBegSentence;
			} else if (strcmp(strLexUnitPronunciation,LEX_UNIT_END_SENTENCE) == 0) {
				return m_lexUnitEndSentence;
			}
			
			if (getLexUnitAndPronunciation(strLexUnitPronunciation,strLexUnit,&iPronunciation) == false) {
				return NULL;
			}
			
			return getLexUnit(strLexUnit,iPronunciation);
		}
		
		// return a lexical unit given its index
		inline LexUnitX *getLexUnit(int iLexUnit) {
		
			assert(m_lexiconX[iLexUnit]->iLexUnit == iLexUnit);		
		
			return m_lexiconX[iLexUnit];
		}
		
		// return a lexical unit given its index
		inline LexUnit *getLexUnitPron(int iLexUnitPron) {
		
			assert(m_lexicon[iLexUnitPron]->iLexUnitPron == iLexUnitPron);
		
			return m_lexicon[iLexUnitPron];
		}
		
		// return a lexical unit given its index
		inline int getLexUnitNoPron(int iLexUnitPron) {
		
			assert(m_lexicon[iLexUnitPron]->iLexUnitPron == iLexUnitPron);
		
			return m_lexicon[iLexUnitPron]->iLexUnit;
		}
		
		// return the canonical pronunciation of a lexical unit
		inline LexUnit *getCanonicalPronunciation(LexUnitX *lexUnitX) {
		
			return lexUnitX->vLexUnitPronunciations.front();
		}
		
		// return the canonical pronunciation of a lexical unit
		inline LexUnit *getCanonicalPronunciation(int iLexUnit) {
		
			return getLexUnit(iLexUnit)->vLexUnitPronunciations.front();
		}
		
		// return a lexical unit as a string of characters
		inline const char *getStrLexUnit(unsigned int iLexUnit) {	
			
			assert(m_lexiconX.size() > iLexUnit);
			return m_lexiconX[iLexUnit]->strLexUnit;
		}
		
		// return a lexical unit as a string of characters
		inline const char *getStrLexUnitPron(unsigned int iLexUnitPron) {	
			
			assert(m_lexicon.size() > iLexUnitPron);
			return m_lexiconX[m_lexicon[iLexUnitPron]->iLexUnit]->strLexUnit;
		}
		
		// return the type of lexical unit
		unsigned char getLexicalUnitType(const char *strLexUnit);
		
		// return whether a lexical unit is standard
		inline bool isStandard(LexUnit *lexUnit) {
			
			return (lexUnit->iType == LEX_UNIT_TYPE_STANDARD);
		}
		
		// return whether a lexical unit is a filler
		inline bool isFiller(LexUnit *lexUnit) {
			
			return (lexUnit->iType == LEX_UNIT_TYPE_FILLER);
		}		
		
		// return whether a lexical unit is a sentence delimiter
		inline bool isSentenceDelimiter(LexUnit *lexUnit) {
			
			return (lexUnit->iType == LEX_UNIT_TYPE_SENTENCE_DELIMITER);
		}	
		
		// return whether a lexical unit is unknown
		inline bool isUnknown(LexUnit *lexUnit) {
			
			return (lexUnit->iType == LEX_UNIT_TYPE_UNKNOWN);
		}	
		
		// return whether a lexical unit is standard
		inline bool isStandard(int iLexUnit) {
			
			return (getLexUnit(iLexUnit)->vLexUnitPronunciations.front()->iType == LEX_UNIT_TYPE_STANDARD);
		}
		
		// return whether a lexical unit is a filler
		inline bool isFiller(int iLexUnit) {
			
			return (getLexUnit(iLexUnit)->vLexUnitPronunciations.front()->iType == LEX_UNIT_TYPE_FILLER);
		}
		
		// return whether a lexical unit is unknown
		inline bool isUnknown(int iLexUnit) {
			
			return (getLexUnit(iLexUnit)->vLexUnitPronunciations.front()->iType == LEX_UNIT_TYPE_UNKNOWN);
		}	
		
		// return whether a lexical unit is a sentence delimiter
		inline bool isSentenceDelimiter(int iLexUnit) {
			
			return (getLexUnit(iLexUnit)->vLexUnitPronunciations.front()->iType == LEX_UNIT_TYPE_SENTENCE_DELIMITER);
		}	
		
		// return a vector containing the all the filler lexical units (special-optional lexical units)
		void getVLexUnitFiller(VLexUnit &vLexUnit);
		
		// compute the average number of alternative pronunciations in the lexicon
		float computeAveragePronunciationVariations();
		
		// comparison function
		static bool compareLexUnitsOrder(const LexUnit *lexUnit1, const LexUnit *lexUnit2) {
		
			if (lexUnit1->iLexUnit == lexUnit2->iLexUnit) {
				return (lexUnit1->iPronunciation < lexUnit2->iPronunciation);
			} else {
				return (lexUnit1->iLexUnit < lexUnit2->iLexUnit);
			}
		}
		
		// comparison function
		inline static bool compareConstChar(const char *str1, const char *str2) {
			
			return strcmp(str1,str2) < 0;
		}
		
		// comparison function
		static bool comparePronProbability(const LexUnit *lexUnit1, const LexUnit *lexUnit2) {
		
			if (lexUnit1->iLexUnit == lexUnit2->iLexUnit) {
				return (lexUnit1->fProbability > lexUnit2->fProbability);
			} else {
				return (lexUnit1->iLexUnit < lexUnit2->iLexUnit);
			}
		}
		
		// return a lexical unit with its pronunciation
		inline bool getLexUnitAndPronunciation(const char *strLexUnitSrc, char *strLexUnit, unsigned char *iPronunciation) {
		
			int iPronunciationAux = 0;
			const char *cAlternative = strrchr(strLexUnitSrc,'(');
			if (strLexUnitSrc[strlen(strLexUnitSrc)-1] != ')')
				cAlternative = NULL;
			if (cAlternative) {
				for(int i = 1 ; cAlternative[i] != ')' ; ++i) {
					// make sure it is a number (btw, this checks that the end of string is not reached before the ')')
					if ((cAlternative[i] > 57) || (cAlternative[i] < 48)) {	
						return false;	
					}
					iPronunciationAux = iPronunciationAux*10;
					iPronunciationAux += cAlternative[i]-48;	
					if (iPronunciationAux >= MAX_LEXUNIT_PRONUNCIATIONS) {
						return false;
					}
				}			
				iPronunciationAux--;
			}
			if (cAlternative) {
				strncpy(strLexUnit,strLexUnitSrc,cAlternative-strLexUnitSrc);
				strLexUnit[cAlternative-strLexUnitSrc] = 0;
			} else {
				strcpy(strLexUnit,strLexUnitSrc);
			}
			*iPronunciation = iPronunciationAux;
			
			return true;
		}
				
		
		// print the a lexical unit with its pronunciation
		inline void print(LexUnit *lexUnit) {
		
			assert(lexUnit);
		
			printf("%-30s (",getStrLexUnit(lexUnit->iLexUnit));
			for(vector<int>::iterator it = lexUnit->vPhones.begin() ; it != lexUnit->vPhones.end() ; ++it) {
				printf(" %s",m_phoneSet->getStrPhone(*it));
			}
			printf(")\n");
		}
		
		// print the a lexical unit (no pronunciation)
		inline void print(LexUnitX *lexUnitX) {
		
			printf("%-30s\n",lexUnitX->strLexUnit);
		}
		
		// print the lexical units with its pronunciations
		inline void print(VLexUnit &vLexUnit) {
		
			for(VLexUnit::iterator it = vLexUnit.begin() ; it != vLexUnit.end() ; ++it) {	
				print(*it);
			}		
		}
		
		// tranform a sequence of lexical units from text format to object format
		bool getLexUnits(const char *strText, VLexUnit &vLexUnit, bool &bAllKnown);
		
		// tranform a sequence of lexical units from text format to object format
		bool getLexUnits(const char *strText, VLexUnitX &vLexUnitX, bool &bAllKnown);
		
		// writes the lexical unit along with its pronunciation number into the given buffer
		void getStrLexUnitPronunciation(LexUnit *lexUnit, string &strLexUnitPronunciation) {	
			
			ostringstream oss;
			oss << getStrLexUnit(lexUnit->iLexUnit);
			if (lexUnit->iPronunciation > 0) {
				oss << "(" << (lexUnit->iPronunciation+1) << ")";
			}
			strLexUnitPronunciation = oss.str();
		}
		
		// map lexical units from pronunciation format to non-pron format
		void map(VLexUnit &vLexUnitInput, VLexUnitX &vLexUnitOutput);		
		
		// remove non-standard lexical units
		void removeNonStandardLexUnits(VLexUnit &vLexUnit);
};

};	// end-of-namespace

#endif
