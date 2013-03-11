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


#include "ReferenceText.h"

namespace Bavieca {

// contructor
ReferenceText::ReferenceText(LexiconManager *lexiconManager, const char *strFileReference) {

	m_lexiconManager = lexiconManager;
	m_strFileReference = strFileReference;
}

// destructor
ReferenceText::~ReferenceText() {

	destroy();
}

// clean-up
void ReferenceText::destroy() {

	for(VReferenceWord::iterator it = m_vReferenceWord.begin() ; it != m_vReferenceWord.end() ; ++it) {
		delete (*it);		
	}	
	m_vReferenceWord.clear();
}

// load the reference text
bool ReferenceText::load() {

   char strWord[MAX_LEXUNIT_LENGTH];
   char strLine[1024+1];
      
   // open the text file
   FILE *file = fopen(m_strFileReference.c_str(),"r");
   if (file == NULL) {
   	printf("file: %s could not be opened\n",m_strFileReference.c_str());
      return false;
   }
   
   int iPosition = 0;
   int iIndexWord = 0;
   int iIndexSentence = 0;
	int iIndexWithinWord = 0;
   while(fgets(strLine,1024,file) != NULL) {
      int iLength = (int)strlen(strLine);		
		for(int i=0; i<iLength ; ++i , ++iPosition) {
			// regular words and initials+acronyms
			if ((isWordCharacter(strLine[i])) || ((strLine[i] == '.') && (i+1 < iLength) && (isalnum(strLine[i+1]) != 0))) {
				strWord[iIndexWithinWord++] = toupper(strLine[i]);
			}
			else {
				// end the current word
				if (iIndexWithinWord > 0) {
					strWord[iIndexWithinWord] = 0;					
					ReferenceWord *word = new ReferenceWord;
					word->iIndex = iIndexWord++;
					word->iIndexSentence = iIndexSentence;
					word->iLexUnit = m_lexiconManager->getLexUnitId(strWord);
					if (word->iLexUnit == -1) {
						return false;
					}
					//printf("%s\n",strWord);
					m_vReferenceWord.push_back(word); 
					iIndexWithinWord = 0;
				}
				// check for an end of sentence 
				// IMP: there can be initials like U.S.A, so after the period there should be a space or the end of file
				bool bEndSentence = false;
				// these characters always signal the end of sentence
				if ((strLine[i] == '?') || (strLine[i] == '!')) {
					bEndSentence = true;
				}
				if (strLine[i] == '.') {
					// no more characters after the period: there is an end of sentence
					if (i+1 >= iLength) {
						bEndSentence = true;
					} 
					// make sure it is not part of initials/acronyms
					else if (isalnum(strLine[i+1]) == 0) {
						bEndSentence = true;
					}
				}	
				if (bEndSentence) {
				   ReferenceSentence *sentence = new ReferenceSentence;
				   if (m_vReferenceSentence.empty()) {
				   	sentence->iIndexFirstWord = 0;
				   } else {
				   	sentence->iIndexFirstWord = m_vReferenceSentence.back()->iIndexLastWord + 1;
				   }
				   // there has to be at least one word
				   if (m_vReferenceWord.empty()) {
				   	return false;
				   }
				   sentence->iIndexLastWord = m_vReferenceWord.back()->iIndex;
				   m_vReferenceSentence.push_back(sentence);
				   iIndexSentence = 0;
				}
			}
		}
		// end the current word
		if (iIndexWithinWord > 0) {
			strWord[iIndexWithinWord] = 0;
			ReferenceWord *word = new ReferenceWord;
			word->iIndex = iIndexWord++;
			word->iIndexSentence = iIndexSentence;
			word->iLexUnit = m_lexiconManager->getLexUnitId(strWord);
			if (word->iLexUnit == -1) {
				return false;
			}
			m_vReferenceWord.push_back(word); 
			iIndexWithinWord = 0;
		}	
   }
   
   if (fclose(file) == EOF) {
   	return false;
   }
          
   return true;
}

// return whether a character is a "word character" (a character that is readable)
bool ReferenceText::isWordCharacter(char c) {

   if ((c >= 48) && (c <= 57))  // digits
      return true;
   if ((c >= 65) && (c <= 99))  // capital letters
      return true;
   if ((c >= 97) && (c <= 122))  // lowercase
      return true;
   if (c == 39)   // "'"
      return true;
   if (c == 45)	// "-" dash 
   	return true;
   // special symbols
   if (c == '<')	
   	return true;
   if (c == '>')	
   	return true;
      
   return false;
}

// convert a word to lower case
void ReferenceText::toLower(char *strWord) {

   for(unsigned int i = 0 ; i < strlen(strWord) ; ++i) {
      strWord[i] = tolower(strWord[i]);
   }
}

// convert a word to upper case
void ReferenceText::toUpper(char *strWord) {

   for(unsigned int i = 0 ; i < strlen(strWord) ; ++i) {
      strWord[i] = toupper(strWord[i]);
   }
}

// print the reference text
void ReferenceText::print() {
	
	printf("------------------------------------------------------\n");
	printf("Reference text:\n");
	for(VReferenceWord::iterator it = m_vReferenceWord.begin() ; it != m_vReferenceWord.end() ; ++it) {	
		printf(" %s",m_lexiconManager->getStrLexUnit((*it)->iLexUnit));
	}
	printf("\n------------------------------------------------------\n");
}

};	// end-of-namespace

