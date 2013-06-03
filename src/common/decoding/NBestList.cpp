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

#include <iostream>

#include "FileOutput.h"
#include "IOBase.h"
#include "NBestList.h"
#include "NBestListEntry.h"

namespace Bavieca {

// constructor
NBestList::NBestList() {
}

// destructor
NBestList::~NBestList() {

	for(VNBestListEntry::iterator it = m_vEntries.begin() ; it != m_vEntries.end() ; ++it) {
		delete *it;
	}
}

// store to disk
void NBestList::store(const char *strFile, bool bTextFormat) {
	
	// text format
	if (bTextFormat) {	
		FileOutput file(strFile,false);
		file.open();
		for(VNBestListEntry::iterator it = m_vEntries.begin() ; it != m_vEntries.end() ; ++it) {
			(*it)->store(file.getStream(),bTextFormat);
		}	
		file.close();	
	}
	// binary format
	else {
		FileOutput file(strFile,true);
		file.open();
		IOBase::write(file.getStream(),(unsigned int)m_vEntries.size(),true);
		for(VNBestListEntry::iterator it = m_vEntries.begin() ; it != m_vEntries.end() ; ++it) {
			(*it)->store(file.getStream());
		}	
		file.close();	
	}
}

// load from disk
void NBestList::load(const char *strFile) {
	
	
	
}

// print
void NBestList::print() {
	
	cout << "- n-best list ------- ( " << m_vEntries.size() << " entries )--------------------------------" << endl;
	for(VNBestListEntry::iterator it = m_vEntries.begin() ; it != m_vEntries.end() ; ++it) {
		(*it)->print();
	}
	cout << "---------------------------------------------------------------------------------------------" << endl;
}


};	// end-of-namespace
