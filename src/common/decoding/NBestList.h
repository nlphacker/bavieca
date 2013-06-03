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

#ifndef NBESTLIST_H
#define NBESTLIST_H

using namespace std;

#include <vector>

namespace Bavieca {

class NBestListEntry;

typedef vector<NBestListEntry*> VNBestListEntry;

/**
	@author root <dani.bolanos@gmail.com>
*/
class NBestList {

	private:
	
		VNBestListEntry m_vEntries;

	public:

		// constructor
		NBestList();

		// destructor
		~NBestList();
		
		// store to disk
		void store(const char *strFile, bool bTextFormat);
		
		// load from disk
		void load(const char *strFile);
		
		// print
		void print();
		
		// add an entry to the list
		void add(NBestListEntry *entry) {
			m_vEntries.push_back(entry);
		}
};

};	// end-of-namespace

#endif
