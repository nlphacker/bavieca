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


#ifndef FILETRN_H
#define FILETRN_H

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "Global.h"

namespace Bavieca {

// entry
typedef struct {
	string strTranscription;
	string strUtteranceId;
} TrnEntry;

typedef vector<TrnEntry*> VTrnEntry;
typedef map<string,TrnEntry*> MTrnEntry;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class TrnFile {

	private:
	
		string m_strFile;
		VTrnEntry m_vTrnEntry;
		MTrnEntry m_mTrnEntry;

	public:
    
		// constructor
		TrnFile(const char *strfile);

		// destructor
		~TrnFile();
		
		// load the file
		void load();
		
		// return the size
		inline int size() {
		
			return m_vTrnEntry.size();
		}
		
		// return the given entry
		inline const char *getTranscription(const char *strUtterance) {
		
			MTrnEntry::iterator it = m_mTrnEntry.find(strUtterance);
			if (it == m_mTrnEntry.end()) {
				return NULL;
			}		
			return it->second->strTranscription.c_str();
		}

};

};	// end-of-namespace

#endif
