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


#include "AlignmentFile.h"
#include "LogMessage.h"
#include "HMMManager.h"
#include "PhoneSet.h"

#include <iomanip>

namespace Bavieca {

// constructor
AlignmentFile::AlignmentFile(PhoneSet *phoneSet, LexiconManager *lexiconManager) {

	m_lexiconManager = lexiconManager;
	m_phoneSet = phoneSet;
}

// destructor
AlignmentFile::~AlignmentFile() {

}

// store an alignment into a file
void AlignmentFile::store(VPhoneAlignment &vPhoneAlignment, const char *strAlignmentFile) {

	assert(m_lexiconManager != NULL);
	
	FileOutput file(strAlignmentFile,false);
	file.open();

	// determine the width of each alignment number
	int iDigits = 0;
	for(int i=10 ; i<INT_MAX ; i*=10) {
		++iDigits;
		if (vPhoneAlignment.back()->iStateEnd[NUMBER_HMM_STATES-1] < i) {
			break;
		}
	}
	
	// determine the width of the phone field
	int iPhoneCharactersMax = 0;
	for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
		int iLength = strlen(m_phoneSet->getStrPhone(i));
		if (iLength > iPhoneCharactersMax) {
			iPhoneCharactersMax = iLength;
		}
	}
	
	// determine the width of the lexical unit field
	int iCharactersLexUnit = 0;
	for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {
		int iLength = strlen(m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit));
		if (iLength > iCharactersLexUnit) {
			iCharactersLexUnit = iLength;
		}
	}
	
	// print the alignment information to the file
	ostringstream oss;
	for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {
		for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			oss << setw(iDigits) << (*it)->iStateBegin[iState] << " " << setw(iDigits) << (*it)->iStateEnd[iState] 
			<< " ";
		}
		if (((*it)->iPosition == WITHIN_WORD_POSITION_START) || ((*it)->iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
			int iPronunciation = (*it)->lexUnit->iPronunciation;
			if (iPronunciation == 0) {
				oss << setw(iPhoneCharactersMax) << m_phoneSet->getStrPhone((*it)->iPhone) << " " 
				<< std::setw(12) << std::setiosflags(ios::fixed) << std::setprecision(4) << (*it)->fLikelihood 
				<< " " << m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit) << endl;
			} else {
				assert(iPronunciation > 0);
				oss << setw(iPhoneCharactersMax) << m_phoneSet->getStrPhone((*it)->iPhone) << " " 
				<< std::setw(12) << std::setiosflags(ios::fixed) << std::setprecision(4) << (*it)->fLikelihood 
				<< " " << m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit) << "(" << (iPronunciation+1) 
				<< ")" << endl;	
			}
		} else {
			oss << setw(iPhoneCharactersMax) << m_phoneSet->getStrPhone((*it)->iPhone) << " " 
			<< std::setw(12) << std::setiosflags(ios::fixed) << std::setprecision(4) << (*it)->fLikelihood << endl;
		}
	}
	IOBase::writeString(file.getStream(),oss);

	file.close();
}

// load an alignment from a file
VPhoneAlignment *AlignmentFile::load(const char *strFile) {

	FileInput file(strFile,false);
	file.open();
	
	VPhoneAlignment *vPhoneAlignment = new VPhoneAlignment; 
	string strPhone,strLexUnit;
	
	unsigned char iPosition = WITHIN_WORD_POSITION_START;
	LexUnit *lexUnit = NULL;
	int iStateBegin[NUMBER_HMM_STATES];
	int iStateEnd[NUMBER_HMM_STATES];
	float fLikelihood;
	
	int iLine = 1;
	string strLine;
   while(std::getline(file.getStream(),strLine)) {
   	std::stringstream ss(strLine);
   	// read fixed fields
		IOBase::read(ss,&iStateBegin[0],false);
		IOBase::read(ss,&iStateEnd[0],false);
		IOBase::read(ss,&iStateBegin[1],false);
		IOBase::read(ss,&iStateEnd[1],false);
		IOBase::read(ss,&iStateBegin[2],false);
		IOBase::read(ss,&iStateEnd[2],false);
		IOBase::readString(ss,strPhone);
		IOBase::read(ss,&fLikelihood,false);
		// starting phone: it is followed by the lexical unit
		if (!ss.eof()) {
			IOBase::readString(ss,strLexUnit);
			lexUnit = NULL;	
			if (m_lexiconManager != NULL) {
				lexUnit = m_lexiconManager->getLexUnitPronunciation(strLexUnit.c_str());
				if (lexUnit == NULL) {
					BVC_ERROR << "lexical unit \"" <<  strLexUnit << 
					"\" was not found int the lexicon, file: " << strFile <<", line: " << iLine;
				}	
			}
			iPosition = WITHIN_WORD_POSITION_MONOPHONE;			
		} else {
			if (iLine == 1) {
				BVC_ERROR << "incorrect format of alignment file: " << strFile << ", line: " << iLine;
			}
			iPosition = WITHIN_WORD_POSITION_INTERNAL;
		}
		PhoneAlignment *phoneAlignment = new PhoneAlignment;
		for(int i=0 ; i < NUMBER_HMM_STATES ; ++i) {
			phoneAlignment->iStateBegin[i] = iStateBegin[i];
			phoneAlignment->iStateEnd[i] = iStateEnd[i];
		}		
		int iPhone = m_phoneSet->getPhoneIndex(strPhone.c_str());
		if (iPhone == -1) {
			BVC_ERROR << "phonetic symbol " << strPhone << " was not found in the phonetic symbol set";
		}
		phoneAlignment->iPhone = iPhone;
		phoneAlignment->fLikelihood = fLikelihood;
		phoneAlignment->lexUnit = lexUnit;
		phoneAlignment->iPosition = iPosition;
		// fix the within-word position of the previous phone
		if (vPhoneAlignment->empty() == false) {
			if ((iPosition == WITHIN_WORD_POSITION_INTERNAL) && 
				(vPhoneAlignment->back()->iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
				vPhoneAlignment->back()->iPosition = WITHIN_WORD_POSITION_START;
			} else if ((iPosition == WITHIN_WORD_POSITION_MONOPHONE) && 
				(vPhoneAlignment->back()->iPosition == WITHIN_WORD_POSITION_INTERNAL)) {
				vPhoneAlignment->back()->iPosition = WITHIN_WORD_POSITION_END;
			}
		}
		vPhoneAlignment->push_back(phoneAlignment);
		++iLine;
	}
	// fix the within-word position of the previous phone
	if ((!vPhoneAlignment->empty()) && (vPhoneAlignment->back()->iPosition != WITHIN_WORD_POSITION_MONOPHONE)) {
		assert(vPhoneAlignment->back()->iPosition == WITHIN_WORD_POSITION_INTERNAL);
		vPhoneAlignment->back()->iPosition = WITHIN_WORD_POSITION_END;
	}	

	file.close();

	return vPhoneAlignment;
}

// prints the alignment
void AlignmentFile::print(VPhoneAlignment &vPhoneAlignment) {

	if (vPhoneAlignment.empty()) {
		printf("<empty alignment>\n");
		return;
	}

	// determine the width of each alignment number
	int iDigits = 0;
	for(int i=10 ; i<INT_MAX ; i*=10) {
		++iDigits;
		if (vPhoneAlignment.back()->iStateEnd[NUMBER_HMM_STATES-1] < i) {
			break;
		}
	}
	
	// determine the width of the phone field
	int iPhoneCharactersMax = 0;
	for(unsigned int i=0 ; i < m_phoneSet->size() ; ++i) {
		int iLength = strlen(m_phoneSet->getStrPhone(i));
		if (iLength > iPhoneCharactersMax) {
			iPhoneCharactersMax = iLength;
		}
	}
	
	for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {	
		for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			printf("%*d %*d ",iDigits,(*it)->iStateBegin[iState],iDigits,(*it)->iStateEnd[iState]);
		}
		char strLexUnitPronunciation[MAX_LEXUNIT_LENGTH+1];
		if (m_lexiconManager != NULL) {
			m_lexiconManager->getStrLexUnitPronunciation((*it)->lexUnit,strLexUnitPronunciation);
		} else {
			strcpy(strLexUnitPronunciation,"<unavailable>");
		}	
		printf("%*s %12.4f (%d) %s %d\n",iPhoneCharactersMax,m_phoneSet->getStrPhone((*it)->iPhone),(*it)->fLikelihood,(*it)->iPosition,strLexUnitPronunciation,(*it)->lexUnit->iLexUnit);
	}
}

// prints the lexical unit alignment
void AlignmentFile::print(VLexUnitAlignment &vLexUnitAlignment) {

	for(VLexUnitAlignment::iterator it = vLexUnitAlignment.begin() ; it != vLexUnitAlignment.end() ; ++it) {
	
		printf("%5d %5d %-20s %10.2f\n",(*it)->iBegin,(*it)->iEnd,m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit),(*it)->fLikelihood);
	}
}

// convert to an object of the class Alignment
Alignment *AlignmentFile::toAlignment(PhoneSet *phoneSet, HMMManager *hmmManager, VPhoneAlignment *vPhoneAlignment) {

	Alignment *alignment = new Alignment(ALIGNMENT_TYPE_VITERBI);
	
	// state-level alignment
	unsigned char iPhoneLeft = phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE);
	for(VPhoneAlignment::iterator it = vPhoneAlignment->begin() ; it != vPhoneAlignment->end() ; ++it) {
		// get the right context
		unsigned char iPhoneRight = UCHAR_MAX;
		VPhoneAlignment::iterator kt = it;
		++kt;
		if (kt != vPhoneAlignment->end()) {
			iPhoneRight = (*kt)->iPhone;
		} else {
			iPhoneRight = phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE);
		}
		for(unsigned char iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			// get the HMM-state 
			int iHMMState = hmmManager->getHMMStateIndex(&iPhoneLeft,(*it)->iPhone,&iPhoneRight,(*it)->iPosition,iState);
			for(int iFrame = (*it)->iStateBegin[iState] ; iFrame <= (*it)->iStateEnd[iState] ; ++iFrame) {
				VStateOcc *vStateOcc = new VStateOcc;
				vStateOcc->push_back(Alignment::newStateOcc(iHMMState,1.0));
				alignment->addFrameAlignmentBack(vStateOcc);
			}
		}	
		iPhoneLeft = (*it)->iPhone;	
	}
	
	// word alignment
	VLexUnitAlignment *vLexUnitAlignment = getVLexUnitAlignment(*vPhoneAlignment);
	for(VLexUnitAlignment::iterator it = vLexUnitAlignment->begin() ; it != vLexUnitAlignment->end() ; ++it) {
		alignment->addLexUnitAlignmentBack((*it)->iBegin,(*it)->iEnd,(*it)->lexUnit);
	}
	destroyLexUnitAlignment(vLexUnitAlignment);
	
	return alignment;
}


// return a lexical unit alignment from a phone alignment
VLexUnitAlignment *AlignmentFile::getVLexUnitAlignment(VPhoneAlignment &vPhoneAlignment) {

	VLexUnitAlignment *vLexUnitAlignment = new VLexUnitAlignment;

	LexUnitAlignment *lexUnitAlignment = NULL;
	for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {
		// there is only one starting phone per word
		if (((*it)->iPosition == WITHIN_WORD_POSITION_START) || ((*it)->iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
			lexUnitAlignment = new LexUnitAlignment();
			lexUnitAlignment->lexUnit = (*it)->lexUnit;	
			if (vLexUnitAlignment->empty() == false) {
				vLexUnitAlignment->back()->iEnd = (*it)->iStateBegin[0]-1;
			}
			lexUnitAlignment->iBegin = (*it)->iStateBegin[0];
			lexUnitAlignment->fLikelihood = (*it)->fLikelihood;
			vLexUnitAlignment->push_back(lexUnitAlignment);
		} else {
			assert(lexUnitAlignment != NULL);
			lexUnitAlignment->fLikelihood += (*it)->fLikelihood;
		}
	}
	if (vLexUnitAlignment->empty() == false) {
		vLexUnitAlignment->back()->iEnd = vPhoneAlignment.back()->iStateEnd[NUMBER_HMM_STATES-1];
	}

	return vLexUnitAlignment;
}

// return a given lexical unit alignment from a phone alignment
LexUnitAlignment *AlignmentFile::getLexUnitAlignment(VPhoneAlignment &vPhoneAlignment, unsigned int iIndex) {

	assert(vPhoneAlignment.size() >= iIndex+1);

	int iLexUnits = -1;
	LexUnitAlignment *lexUnitAlignment = NULL;
	for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {
		// there is only one starting phone per word
		if (((*it)->iPosition == WITHIN_WORD_POSITION_START) || ((*it)->iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
			++iLexUnits;
			if (iLexUnits == (int)iIndex) {
				lexUnitAlignment = new LexUnitAlignment();
				lexUnitAlignment->lexUnit = (*it)->lexUnit;	
				lexUnitAlignment->iBegin = (*it)->iStateBegin[0];
				lexUnitAlignment->fLikelihood = (*it)->fLikelihood;
			} else if (iLexUnits == (int)iIndex+1) {
				assert(lexUnitAlignment != NULL);
				lexUnitAlignment->iEnd = (*it)->iStateBegin[0]-1;
			}
		} else {
			if (iLexUnits == (int)iIndex) {
				assert(lexUnitAlignment != NULL);
				lexUnitAlignment->fLikelihood += (*it)->fLikelihood;
			}
		}
	}
	if (vPhoneAlignment.size() == iIndex+1) {
		lexUnitAlignment->iEnd = vPhoneAlignment.back()->iStateEnd[NUMBER_HMM_STATES-1];
	}

	return lexUnitAlignment;
}

// prints the alignment
void AlignmentFile::print(PhoneSet *phoneSet, LexiconManager *lexiconManager, VPhoneAlignment &vPhoneAlignment) {

	if (vPhoneAlignment.empty()) {
		printf("<empty alignment>\n");
	}

	// determine the width of each alignment number
	int iDigits = 0;
	for(int i=10 ; i<INT_MAX ; i*=10) {
		++iDigits;
		if (vPhoneAlignment.back()->iStateEnd[NUMBER_HMM_STATES-1] < i) {
			break;
		}
	}
	
	// determine the width of the phone field
	int iPhoneCharactersMax = 0;
	for(unsigned int i=0 ; i < phoneSet->size() ; ++i) {
		int iLength = strlen(phoneSet->getStrPhone(i));
		if (iLength > iPhoneCharactersMax) {
			iPhoneCharactersMax = iLength;
		}
	}
	
	for(VPhoneAlignment::iterator it = vPhoneAlignment.begin() ; it != vPhoneAlignment.end() ; ++it) {	
		for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			printf("%*d %*d ",iDigits,(*it)->iStateBegin[iState],iDigits,(*it)->iStateEnd[iState]);
		}
		string strLexUnit = lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit);
		printf("%*s %.4f (%d) %-s\n",iPhoneCharactersMax,phoneSet->getStrPhone((*it)->iPhone),(*it)->fLikelihood,(*it)->iPosition,strLexUnit.c_str());
	}
}

// prints the lexical unit alignment
void AlignmentFile::print(PhoneSet *phoneSet, LexiconManager *lexiconManager, VLexUnitAlignment &vLexUnitAlignment) {

	float fLikelihoodTotal = 0.0;

	for(VLexUnitAlignment::iterator it = vLexUnitAlignment.begin() ; it != vLexUnitAlignment.end() ; ++it) {	
		printf("%5d %5d %-20s %10.2f\n",(*it)->iBegin,(*it)->iEnd,lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit),(*it)->fLikelihood);
		fLikelihoodTotal += (*it)->fLikelihood;
	}
	printf("total likelihood: %.2f\n",fLikelihoodTotal);
}

};	// end-of-namespace


