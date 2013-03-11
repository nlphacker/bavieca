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
#include "IOBase.h"
#include "TrnFile.h"
#include <string.h>

namespace Bavieca {

// constructor
TrnFile::TrnFile(const char *strFile)
{
	m_strFile = strFile;
}

// destructor
TrnFile::~TrnFile()
{
	for(VTrnEntry::iterator it = m_vTrnEntry.begin() ; it != m_vTrnEntry.end() ; ++it) {
		delete *it;
	}
	m_vTrnEntry.clear();
	m_mTrnEntry.clear();
}

// load the file
void TrnFile::load() {

	string strLine;
	int iLine = 0;
	
	FileInput file(m_strFile.c_str(),false);
	file.open();
	
	// (2) for each utterance, read the feature file name and the associated lexical units
	while(std::getline(file.getStream(),strLine)) {
		++iLine;

		size_t i = strLine.find_last_of('(');
		string strTranscription = strLine.substr(0,i);
		string strUtteranceId = strLine.substr(i);
		if ((strUtteranceId.length() < 3) || (strUtteranceId.substr(strUtteranceId.length()-1).compare(")"))) {
			BVC_ERROR << "wrong utterance-id format found int trn file, line: " << iLine;
		} 
		strUtteranceId = strUtteranceId.substr(1,strUtteranceId.length()-2);

		TrnEntry *trnEntry = new TrnEntry;
		trnEntry->strTranscription = strTranscription;
		trnEntry->strUtteranceId = strUtteranceId;
		m_vTrnEntry.push_back(trnEntry);
		m_mTrnEntry.insert(MTrnEntry::value_type(trnEntry->strUtteranceId.c_str(),trnEntry));
	}

	file.close();
}

};	// end-of-namespace

