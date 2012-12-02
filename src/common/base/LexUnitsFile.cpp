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


#include "LexUnitsFile.h"
#include "Global.h"
#include "LogMessage.h"
#include "FileInput.h"

namespace Bavieca {

// constructor
LexUnitsFile::LexUnitsFile(LexiconManager *lexiconManager, const char *strFile) {

	m_lexiconManager = lexiconManager;
	m_strFile = strFile;	
}

// destructor
LexUnitsFile::~LexUnitsFile() {

	m_vLexUnit.clear();
}

// load the lexical units from the file
void LexUnitsFile::load() {

	string strLine;
	string strLexUnit;
	
	FileInput file(m_strFile.c_str(),false);
	file.open();
 
	int iLine = 1;
   while(std::getline(file.getStream(),strLine)) {
      VLexUnit vLexUnitLine;
      bool bAllKnown;
      if (m_lexiconManager->getLexUnits(strLine.c_str(),vLexUnitLine,bAllKnown) == false) {
      	BVC_ERROR<< "loading the lexical units from file: " << m_strFile << " line: " << iLine;
      }
      for(VLexUnit::iterator it = vLexUnitLine.begin() ; it != vLexUnitLine.end() ; ++it) {
      	m_vLexUnit.push_back(*it);
      }
		++iLine;
   }
   
	file.close();
}

// return a lexical unit
LexUnit *LexUnitsFile::getLexUnit(const char *strLexUnit, int iLine) {

	// extract the pronunciation number from the lexical unit
	int iPronunciation = 0;	
	char *cAlternative = (char*)strrchr(strLexUnit,'(');
	if (strLexUnit[strlen(strLexUnit)-1] != ')')
		cAlternative = NULL;
	if (cAlternative) {
		for(int i = 1 ; cAlternative[i] != ')' ; ++i) {
			// make sure it is a number (btw, this checks that the end of string is not reached before the ')')
			if ((cAlternative[i] > 57) || (cAlternative[i] < 48)) {
				BVC_ERROR<< "lexical unit " << strLexUnit << " is not in correct format: line " << iLine;
			}
			iPronunciation = iPronunciation*10;
			iPronunciation += cAlternative[i]-48;	
		}			
		iPronunciation--;
		*cAlternative = '\0';
	} 
	// get the lexical unit from the lexicon
	LexUnit *lexUnit = m_lexiconManager->getLexUnit(strLexUnit,iPronunciation);
	if (lexUnit == NULL) {
		BVC_ERROR<< "lexical unit: " << strLexUnit << " was not found in the lexicon, line: " << iLine;
	}
	
	return lexUnit;	
}

// print the lexical units
void LexUnitsFile::print() {

	printf("\n");
	for(VLexUnit::iterator it = m_vLexUnit.begin() ; it != m_vLexUnit.end() ; ++it) {
		printf("%s ",m_lexiconManager->getStrLexUnit((*it)->iLexUnit));	
	}
	printf("\n");
}

};	// end-of-namespace

