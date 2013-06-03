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


#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "LMARPA.h"
#include "LMFSM.h"
#include "LMManager.h"
#include "LogMessage.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
LMManager::LMManager(LexiconManager *lexiconManager, const char *strFile, 
	const char *strFormat, const char *strType) {

	m_lexiconManager = lexiconManager;
	m_strFile = strFile;
	m_strFormat = strFormat;
	m_strType = strType;
	m_iNGram = -1;
	
	m_bLMLoaded = false;
	m_lmARPA = NULL;
	m_lmFSM = NULL;
}

// destructor
LMManager::~LMManager() {

	if (m_bLMLoaded) {
		delete m_lmARPA;
		delete m_lmFSM;
	}
}

// load the language model from disk
void LMManager::load() {
	
	double dStartTime = TimeUtils::getTimeMilliseconds();

	// ARPA
	if (m_strFormat.compare(LM_FILE_FORMAT_ARPA) == 0) {
		m_lmARPA = new LMARPA(m_lexiconManager,m_strFile.c_str());
		m_lmARPA->load();
		m_lmFSM = new LMFSM(m_lexiconManager);
		m_lmFSM->build(m_lmARPA);
	}
	// Finite State Machine (FSM)
	else if (m_strFormat.compare(LM_FILE_FORMAT_FSM) == 0) {
		m_lmFSM = new LMFSM(m_lexiconManager);
		m_lmFSM->load(m_strFile.c_str());
	}
	// not supported
	else {
		BVC_ERROR << "language model format: " << m_strFormat << " is not supported" << endl; 
	}
	
	double dEndTime = TimeUtils::getTimeMilliseconds();
	double dSeconds = (dEndTime-dStartTime)/1000;	
	BVC_VERB << "LM loading time: " << dSeconds << " seconds";	
	
	m_bLMLoaded = true;
}



};	// end-of-namespace
