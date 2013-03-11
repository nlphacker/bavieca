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


#include "FillerManager.h"
#include "BatchFile.h"
#include "LogMessage.h"

namespace Bavieca {

// constructor
FillerManager::FillerManager(const char *strFile) {
	
	m_strFile = strFile;
}

// destructor
FillerManager::~FillerManager() {

	for(VFiller::iterator it = m_vFiller.begin() ; it != m_vFiller.end() ; ++it) {
		delete *it;
	}
	m_vFiller.clear();
}

// load the fillers
void FillerManager::load() {

	BatchFile batchFile(m_strFile.c_str(),"lexUnit|insertionPenalty");
	batchFile.load();	
	
	for(unsigned int i=0 ; i < batchFile.size() ; ++i) {
		const char *strLexUnit = batchFile.getField(i,"lexUnit");
		const char *strIP = batchFile.getField(i,"insertionPenalty");
		float fIP = (float)atof(strIP);
	
		// check format
		unsigned int iLength = (unsigned int)strlen(strLexUnit);
		if ((iLength < 3) || (strLexUnit[iLength-1] != '>') || (strLexUnit[0] != '<')) {
			BVC_ERROR << "filler " << strLexUnit << " is not in the right format";
		}	
		// keep the filler
		Filler *filler = new Filler;
		filler->strLexUnit = strLexUnit;
		filler->fInsertionPenalty = fIP;
		m_vFiller.push_back(filler);
   }
}

// attach insertion penalties to filler lexical units in the lexicon
void FillerManager::attachInsertionPenaltyFillers(LexiconManager *lexiconManager) {
	
	for(VFiller::iterator it = m_vFiller.begin() ; it != m_vFiller.end() ; ++it) {

		LexUnitX *lexUnitX = lexiconManager->getLexUnit((*it)->strLexUnit.c_str());
		if (lexUnitX == NULL) {
			BVC_ERROR << "filler lexical unit " << (*it)->strLexUnit.c_str() << " is not defined in the lexicon";
		}
		for(VLexUnit::iterator jt = lexUnitX->vLexUnitPronunciations.begin() ; jt != lexUnitX->vLexUnitPronunciations.end() ; ++jt) {
			(*jt)->fInsertionPenalty = (*it)->fInsertionPenalty;
		}
	}
}

};	// end-of-namespace

