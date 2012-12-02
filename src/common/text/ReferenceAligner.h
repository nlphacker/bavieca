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


#ifndef REFERENCEALIGNER_H
#define REFERENCEALIGNER_H

#include "BestPath.h"
#include "LexiconManager.h"
#include "ReferenceAlignment.h"
#include "ReferenceText.h"

namespace Bavieca {

typedef struct _GridElement {
   int flag;         // filled in during DP search through grid
   int score;			// alignment score
   int ref;				// reference word index
   int tra;				// transcription word index
   _GridElement *prev,*next;
   int bookword;        /* filled with correct vals for final path */
} GridElement; 

using namespace std;

// penalties used for the alignment
#define WA_INS_PEN 3    /* insertion penalty during alignment */
#define WA_DEL_PEN 3    /* deletion penalty */
#define WA_SUB_PEN 4    /* substitution penalty */

enum {WA_INS=0x1, WA_DEL=0x2, WA_SUB=0x4, WA_COR=0x8};

/**
	@author root <root@localhost.localdomain>
*/
class ReferenceAligner {

	private:
	
		LexiconManager *m_lexiconManager;

	public:

		// constructor
		ReferenceAligner(LexiconManager *lexiconManager);

		// destructor
		~ReferenceAligner();
		
		// align a recognition path against a sequence of lexical units
		ReferenceAlignment *align(BestPath *bestPath, ReferenceText *referenceText);

};

};	// end-of-namespace

#endif
