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

#include "IOBase.h"
#include "LexiconManager.h"
#include "NBestListEntryElement.h"

namespace Bavieca {

// constructor
NBestListEntryElement::NBestListEntryElement(int iBegin, int iEnd, LexUnit *lexUnit, 
	double dAM, double dLM, double dIP, double dPP) {
	m_iBegin = iBegin;
	m_iEnd = iEnd;
	m_lexUnit = lexUnit;
	m_dLikelihoodAM = dAM;
	m_dLikelihoodLM = dLM;
	m_dIP = dIP;
	m_dPP = dPP;
}

// destructor
NBestListEntryElement::~NBestListEntryElement() {	
}

// store to disk
void NBestListEntryElement::store(ostream &os, bool bTextFormat, LexiconManager *lexiconManager) {

	// text format
	if (bTextFormat) {
		assert(lexiconManager);
		ostringstream oss;
		string strPron;
		lexiconManager->getStrLexUnitPronunciation(m_lexUnit,strPron);
		oss << strPron;		
		IOBase::writeString(os,oss);
	}
	// binary format
	else {
		IOBase::write(os,m_iBegin,true);	
		IOBase::write(os,m_iEnd,true);
		IOBase::write(os,m_lexUnit->iIndex,true);
		IOBase::write(os,m_dLikelihoodAM,true);
		IOBase::write(os,m_dLikelihoodLM,true);
		IOBase::write(os,m_dIP,true);
		IOBase::write(os,m_dPP,true);
	}
}

};	// end-of-namespace
