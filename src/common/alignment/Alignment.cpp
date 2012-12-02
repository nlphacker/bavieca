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


#include "Alignment.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "HMMManager.h"

namespace Bavieca {

// constructor
Alignment::Alignment(unsigned char iType)
{
	m_iType = iType;
	m_bWordLevelAlignment = false;
}

// destructor
Alignment::~Alignment()
{
	for(VFrameAlignment::iterator it = m_vFrameAlignment.begin() ; it != m_vFrameAlignment.end() ; ++it) {
		for(VStateOcc::iterator jt = (*it)->begin() ; jt != (*it)->end() ; ++jt) {
			delete *jt;
		}
		delete *it;
	}
	for(VWordAlignment::iterator it = m_vWordAlignment.begin() ; it != m_vWordAlignment.end() ; ++it) {
		delete *it;	
	}
}

// add a lex-unit alignment
void Alignment::addLexUnitAlignmentFront(int iFrameBegin, int iFrameEnd, LexUnit *lexUnit) {

	m_vWordAlignment.push_front(newWordAlignment(iFrameBegin,iFrameEnd,lexUnit));
	m_bWordLevelAlignment = true;
}

// add a lex-unit alignment
void Alignment::addLexUnitAlignmentBack(int iFrameBegin, int iFrameEnd, LexUnit *lexUnit) {

	m_vWordAlignment.push_back(newWordAlignment(iFrameBegin,iFrameEnd,lexUnit));
	m_bWordLevelAlignment = true;
}

// store to disk
void Alignment::store(const char *strFile) {

	// create the file
	FileOutput file(strFile,true);
	file.open();
			
	IOBase::write(file.getStream(),m_iType);
	IOBase::write(file.getStream(),m_bWordLevelAlignment);
	int iSize = m_vFrameAlignment.size();
	IOBase::write(file.getStream(),iSize);

	for(VFrameAlignment::iterator it = m_vFrameAlignment.begin() ; it != m_vFrameAlignment.end() ; ++it) {
		int iSize = (*it)->size();
		IOBase::write(file.getStream(),iSize);
		for(VStateOcc::iterator jt = (*it)->begin() ; jt != (*it)->end() ; ++jt) {
			IOBase::write(file.getStream(),(*jt)->iHMMState);
			IOBase::write(file.getStream(),(*jt)->dOccupation);
		}
	}
	
	// store the word-level alignment
	if (m_bWordLevelAlignment) {
		iSize = m_vWordAlignment.size();
		IOBase::write(file.getStream(),iSize);
		for(VWordAlignment::iterator it = m_vWordAlignment.begin() ; it != m_vWordAlignment.end() ; ++it) {
			IOBase::write(file.getStream(),(*it)->iFrameBegin);
			IOBase::write(file.getStream(),(*it)->iFrameEnd);
			IOBase::write(file.getStream(),(*it)->iIndex);
		}
	}
	 
   file.close();
}

// load from disk
Alignment *Alignment::load(const char *strFile, LexiconManager *lexiconManager) {

	// open the file
	FileInput file(strFile,true);
	file.open();
		
	int iFrames = -1;
	int iStates = -1;
	unsigned char iType;
	
	IOBase::read(file.getStream(),&iType);
	Alignment *alignment = new Alignment(iType);
	IOBase::read(file.getStream(),&alignment->m_bWordLevelAlignment);
	IOBase::read(file.getStream(),&iFrames);
	
	for(int i=0 ; i < iFrames ; ++i) {
		IOBase::read(file.getStream(),&iStates);
		VStateOcc *vStateOcc = new VStateOcc;
		for(int j=0 ; j < iStates ; ++j) {
			StateOcc *stateOcc = new StateOcc;
			IOBase::read(file.getStream(),&stateOcc->iHMMState);
			IOBase::read(file.getStream(),&stateOcc->dOccupation);
			vStateOcc->push_back(stateOcc);
		}
		alignment->m_vFrameAlignment.push_back(vStateOcc);
	}
	
	// load the word-level alignment
	int iWords = -1;
	IOBase::read(file.getStream(),&iWords);
	for(int i=0 ; i < iWords ; ++i) {
		WordAlignment *wordAlignment = new WordAlignment;		
		IOBase::read(file.getStream(),&wordAlignment->iFrameBegin);
		IOBase::read(file.getStream(),&wordAlignment->iFrameEnd);
		IOBase::read(file.getStream(),&wordAlignment->iIndex);
		if (lexiconManager) {
			wordAlignment->iLexUnitPron = lexiconManager->getLexUnitByIndex(wordAlignment->iIndex)->iLexUnitPron;
		} else {
			wordAlignment->iLexUnitPron = -1;
		}
		alignment->m_vWordAlignment.push_back(wordAlignment);
	}	
	 
   file.close();

	return alignment;
}

// print
void Alignment::print(LexiconManager *lexiconManager) {

	int t = 0;
	for(VFrameAlignment::iterator it = m_vFrameAlignment.begin() ; it != m_vFrameAlignment.end() ; ++it, ++t) {
		printf("t= %4d\n",t);
		for(VStateOcc::iterator jt = (*it)->begin() ; jt != (*it)->end() ; ++jt) {
			printf("state: %6d occ: %12.4f\n",(*jt)->iHMMState,(*jt)->dOccupation);
		}
	}
	for(VWordAlignment::iterator it = m_vWordAlignment.begin() ; it != m_vWordAlignment.end() ; ++it) {
		if (lexiconManager == NULL) {
			printf("%8d %8d %8d\n",(*it)->iFrameBegin,(*it)->iFrameEnd,(*it)->iLexUnitPron);
		} else if ((*it)->iLexUnitPron != -1) {
			char strLexUnitPron[MAX_LEXUNIT_LENGTH+1];
			lexiconManager->getStrLexUnitPronunciation(
				lexiconManager->getLexUnitPron((*it)->iLexUnitPron),strLexUnitPron);
			printf("%8d %8d %32s\n",(*it)->iFrameBegin,(*it)->iFrameEnd,strLexUnitPron);
		}
	}
}

// convert to phone-alignment
VPhoneAlignment *Alignment::getPhoneAlignment(LexiconManager *lexiconManager) {

	assert(m_iType == ALIGNMENT_TYPE_VITERBI);

	int t = 0;
	VPhoneAlignment *vPhoneAlignment = new VPhoneAlignment;
	for(VWordAlignment::iterator it = m_vWordAlignment.begin() ; it != m_vWordAlignment.end() ; ++it) {
		LexUnit *lexUnit = lexiconManager->getLexUnitPron((*it)->iLexUnitPron);
		if (lexUnit == NULL) {
			return NULL;	
		}
		unsigned int iPhone = 0;
		while(iPhone < lexUnit->vPhones.size()) {
			PhoneAlignment *phoneAlignment = new PhoneAlignment;
			vPhoneAlignment->push_back(phoneAlignment);
			phoneAlignment->lexUnit = lexUnit;
			phoneAlignment->iPhone = lexUnit->vPhones[iPhone];
			if (iPhone == 0) {
				if (lexUnit->vPhones.size() == 1) {
					phoneAlignment->iPosition = WITHIN_WORD_POSITION_MONOPHONE;
				} else {
					phoneAlignment->iPosition = WITHIN_WORD_POSITION_START;
				}
			} else if (iPhone == lexUnit->vPhones.size()-1) {
				phoneAlignment->iPosition = WITHIN_WORD_POSITION_END;
			} else {
				phoneAlignment->iPosition = WITHIN_WORD_POSITION_INTERNAL;
			}
			phoneAlignment->fLikelihood = 0.0;
			for(int i=0 ; i< NUMBER_HMM_STATES ; ++i) {
				phoneAlignment->fLikelihoodState[i] = 0.0;
			}
			int iHMMState = (*m_vFrameAlignment[t])[0]->iHMMState;
			int iState = 0;
			int iStatePrev = -1;
			while(1) {
				// no more frames?
				if (t == (int)m_vFrameAlignment.size()) {
					phoneAlignment->iStateEnd[iStatePrev] = t-1;
					break;
				}
				// new state?
				if ((*m_vFrameAlignment[t])[0]->iHMMState != iHMMState) {
					iHMMState = (*m_vFrameAlignment[t])[0]->iHMMState;
					iState++;
					iState %= NUMBER_HMM_STATES;
				}
				if ((iState == 0) && (iStatePrev == NUMBER_HMM_STATES-1)) {
					phoneAlignment->iStateEnd[iStatePrev] = t-1;
					break;
				}
				if (iState != iStatePrev) {
					if (iStatePrev != -1) {
						phoneAlignment->iStateEnd[iStatePrev] = t-1;
					}
					phoneAlignment->iStateBegin[iState] = t;
					iStatePrev = iState;	
				}
				++t;
			}
			++iPhone;
		}
	}

	return vPhoneAlignment;
}

};	// end-of-namespace

