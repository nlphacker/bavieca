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
#include <iomanip>
 
#include "Global.h"
#include "IOBase.h"
#include "NBestListEntry.h"

namespace Bavieca {

// constructor
NBestListEntry::NBestListEntry(LexiconManager *lexiconManager) {
	m_lexiconManager = lexiconManager;
	m_dLikelihood = 0.0;
	m_dPP = 0.0;
}

// destructor
NBestListEntry::~NBestListEntry() {
	for(VNBestListEntryElement::iterator it = m_vElements.begin() ; it != m_vElements.end() ; ++it) {
		delete *it;
	}
}

// store to disk
void NBestListEntry::store(ostream &os, bool bTextFormat) {
	 
	if (bTextFormat) {
		ostringstream oss;
		oss << FLT(12,4) << m_dLikelihood;
		IOBase::writeString(os,oss);
		for(VNBestListEntryElement::iterator it = m_vElements.begin() ; it != m_vElements.end() ; ++it) {
			string str = " ";
			IOBase::writeString(os,str);
			(*it)->store(os,bTextFormat,m_lexiconManager);
		}
		string str = "\n";
		IOBase::writeString(os,str);
	} else {
		IOBase::write(os,(unsigned int)m_vElements.size(),true);
		for(VNBestListEntryElement::iterator it = m_vElements.begin() ; it != m_vElements.end() ; ++it) {
			(*it)->store(os);
		}
	}
}

// print the entry
void NBestListEntry::print() {
 
	cout << FLT(12,4) << m_dLikelihood << " ";
	for(VNBestListEntryElement::iterator it = m_vElements.begin() ; it != m_vElements.end() ; ++it) {
		cout << m_lexiconManager->getStrLexUnit((*it)->getLexUnit()->iLexUnit) << " ";
	}
	cout << endl;
}

};	// end-of-namespace

