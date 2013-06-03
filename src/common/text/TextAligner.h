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


#ifndef TEXTALIGNER_H
#define TEXTALIGNER_H

#include "LexiconManager.h"
#include "TextAlignment.h"

namespace Bavieca {

// grid element for the text alignment
typedef struct _TAGridElement {
	int iEvent;							// alignment event
	int iScore;							// alignment score
	int iRef;							// reference index
	int iTra;							// transcripton index
	_TAGridElement *prev,*next;
	int iBookWord;
} TAGridElement; 

// penalty for each alignment error
#define TEXT_ALIGNMENT_PENALTY_SUBSTITUTION			4
#define TEXT_ALIGNMENT_PENALTY_DELETION				3
#define TEXT_ALIGNMENT_PENALTY_INSERTION				3

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class TextAligner {

	private:
	
		LexiconManager *m_lexiconManager;
	
	public:

		// constructor
		TextAligner(LexiconManager *lexiconManager);

		// destructor
		~TextAligner();
		
		// align two sequences of lexical units
		TextAlignment *align(VLexUnit &vLexUnitHyp, VLexUnit &vLexUnitRef, 
			bool bPenalizeTrailingDeletions = true);

};

};	// end-of-namespace

#endif
