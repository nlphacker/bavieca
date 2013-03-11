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


#ifndef LEXUNITSFILE_H
#define LEXUNITSFILE_H

using namespace std;

#include <string>

#include "LexiconManager.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class LexUnitsFile {

	private:
	
		string m_strFile;
		LexiconManager *m_lexiconManager;
		VLexUnit m_vLexUnit;
		
		// return a lexical unit
		LexUnit *getLexUnit(const char *strLexUnit, int iLine);		

	public:
    
		// constructor
		LexUnitsFile(LexiconManager *lexiconManager, const char *strFile);

		// destructor
		~LexUnitsFile();
		
		// load the lexical units from the file
		void load();
		
		// return the number of lexical units
		unsigned int size() {
			
			return (unsigned int)m_vLexUnit.size();
		}
		
		// return the lexical unit at the given position
		inline LexUnit *operator[ ](int iPosition) {
		
			return m_vLexUnit[iPosition]; 
		}
		
		// return a reference to the lexical units
		inline VLexUnit *getLexUnits() {
		
			return &m_vLexUnit;
		}
		
		// return the lexical units
		inline void getLexUnits(VLexUnit &vLexUnit) {
		
			for(VLexUnit::iterator it = m_vLexUnit.begin() ; it != m_vLexUnit.end() ; ++it) {
				vLexUnit.push_back(*it);
			}
		}
		
		// print the lexical units
		void print();	
};

};	// end-of-namespace

#endif
