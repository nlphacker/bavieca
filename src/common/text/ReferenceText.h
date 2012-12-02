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


#ifndef REFERENCETEXT_H
#define REFERENCETEXT_H

using namespace std;

#include <string>

#include "LexiconManager.h"

namespace Bavieca {

typedef struct {
	int iIndex;					// sentence index within the whole text
	int iIndexFirstWord;		// index of the first word within the sentence
	int iIndexLastWord;		// index of the last word within the sentence
} ReferenceSentence;

typedef vector<ReferenceSentence*> VReferenceSentence;

typedef struct {
	int iIndex;             // word index in the whole text 
   int iLexUnit;				// lexical unit
   int iIndexSentence;		// sentence index
} ReferenceWord;

typedef vector<ReferenceWord*> VReferenceWord;

/**
	@author root <root@localhost.localdomain>
*/
class ReferenceText {

	private:
	
		LexiconManager *m_lexiconManager;
		string m_strFileReference;
		VReferenceWord	m_vReferenceWord;					// words in the reference text
		VReferenceSentence m_vReferenceSentence;		//	sentences in the reference text

	public:

		// contructor
		ReferenceText(LexiconManager *lexiconManager, const char *strFileReference);

		// destructor
		~ReferenceText();
		
		// load the reference text
		bool load();
		
		// return the number of elements
		inline int size() {
			
			return m_vReferenceWord.size();
		}
		
		// return whether there are elements
		inline bool empty() {
			
			return m_vReferenceWord.empty();
		}		
		
		// return the element at the given position
		inline ReferenceWord *operator[ ](int iPosition) {
			
			return m_vReferenceWord[iPosition];
		}
		
		// return the element at a given position
		inline int getLexUnitAt(unsigned int iPosition) {
		
			if (iPosition >= m_vReferenceWord.size()) {
				return m_lexiconManager->m_lexUnitUnknown->iLexUnit;
			}
			return m_vReferenceWord[iPosition]->iLexUnit;
		}
		
		// clean up
		void destroy();
		
		// return whether a character is a "word character" (a character that is readable)
		bool isWordCharacter(char c);
		
		// convert a word to lower case
		void toLower(char *strWord);
		
		// convert a word to upper case
		void toUpper(char *strWord);
		
		// print the reference text
		void print();
		
};

};	// end-of-namespace

#endif
