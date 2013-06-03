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


#ifndef TEXTALIGNMENT_H
#define TEXTALIGNMENT_H

using namespace std;
#include <vector>

#include "LexiconManager.h"

namespace Bavieca {

// alignment events
enum {TEXT_ALIGNMENT_EVENT_CORRECT=1, TEXT_ALIGNMENT_EVENT_SUBSTITUTION=2, 
	TEXT_ALIGNMENT_EVENT_DELETION=4, TEXT_ALIGNMENT_EVENT_INSERTION=8};

// # of alignment events
#define TEXT_ALIGNMENT_EVENTS_NUMBER				4

// text alignment element
typedef struct {
	int iAlignmentEvent;				// alignment event (correct, ins, del, sub)
	int iIndexReference;				// index within the reference (if applicable)
	int iLexUnitReference;			// lexical unit in the reference (if applicable)
	int iIndexHypothesis;			// index within the hypothesis (if applicable)
	int iLexUnitHypothesis;			// lexical unit in the hypothesis (if applicable)
} TAElement;

typedef vector<TAElement*> VTAElement;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class TextAlignment {

	private:
	
		LexiconManager *m_lexiconManager;										// lexicon manager
		int m_iAlignmentEvents[TEXT_ALIGNMENT_EVENTS_NUMBER];				// alignment events
		VTAElement m_vTAElement;													// alignment elements
		int m_iScore;																	// alignment score
		int m_iWordsReference;
		int m_iWordsHypothesis;

	public:

		// constructor
		TextAlignment(LexiconManager *lexiconManager);

		// destructor
		~TextAlignment();
		
		// adds a new element to the alignment
		void addElement(int iAlignmentEvent, int iIndexReference, int iLexUnitReference, 
			int iIndexHypothesis, int iLexUnitHypothesis);
			
		// return the number of correct words
		inline int getCorrect() {
		
			return m_iAlignmentEvents[getIndex(TEXT_ALIGNMENT_EVENT_CORRECT)];
		}
		
		// get the event index
		inline int getIndex(int iEvent) {
			
			switch(iEvent) {
				case 1: return 0;
				case 2: return 1;
				case 4: return 2;
				case 8: return 3;
				default: { 
					assert(0);
					return -1;
				}
			}
		}
		
		// set the alignment score
		inline void setScore(int iScore) {
		
			m_iScore = iScore;
		}
		
		// print the alignment
		void print(bool bDetailed, ostream &os);
		
		// store the alignment into the given file
		void store(const char *strFile, bool bDetailed);	
		
		// return the alignment as a string
		char *getStrAlignment();
		
		// return the alignment elements
		inline VTAElement &getAlignment() {
		
			return m_vTAElement;
		}

};

};	// end-of-namespace

#endif
