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


#ifndef FILLERMANAGER_H
#define FILLERMANAGER_H

using namespace std;

#include <string>
#include <vector>

#include "Global.h"
#include "LexiconManager.h"

namespace Bavieca {

typedef struct {
	string strLexUnit;			// lexical unit in string format
	float fInsertionPenalty;	// penalty resulting from inserting the filler during decoding
} Filler;

typedef vector<Filler*> VFiller;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FillerManager {

	private:
	
		string m_strFile;											// file containing filler penalties in case they are defined
		VFiller m_vFiller;										// array of fillers with corresponding penalties

	public:
    
    	// constructor
		FillerManager(const char *strFile);

		// destructor
		~FillerManager();
		
		// load the fillers
		void load();
		
		// attach insertion penalties to lexical units in the lexicon
		void attachInsertionPenaltyFillers(LexiconManager *lexiconManager);
};

};	// end-of-namespace

#endif
