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


#include "LexiconManager.h"
#include "BestPath.h"

#include <sstream>
#include <iomanip>

namespace Bavieca {

// constructor
BestPath::BestPath(LexiconManager *lexiconManager, float fScore)
{
	m_lexiconManager = lexiconManager;
	m_iFrameStart = 0;
	m_fScore = fScore;
}

// destructor
BestPath::~BestPath()
{
	for( list<BestPathElement*>::iterator it = m_lBestPathElements.begin() ; it != m_lBestPathElements.end() ; ++it ) {
		delete *it;
	}
	
	m_lBestPathElements.clear();
}

// print a bestPath
void BestPath::print(bool bExtended) {

	for( list<BestPathElement*>::iterator it = m_lBestPathElements.begin() ; it != m_lBestPathElements.end() ; ++it ) {	
		string strLexUnit;
		m_lexiconManager->getStrLexUnitPronunciation((*it)->lexUnit,strLexUnit);
		cout << setw(5) << (*it)->iFrameStart << setw(5) << (*it)->iFrameEnd << setw(20) << strLexUnit;
		if (bExtended) {
			cout << " (gl= " << setw(12) << std::setiosflags(ios::fixed) << std::setprecision(4) << (*it)->fScore 
			<< ") (am= " << setw(12) << std::setprecision(4) << (*it)->fScoreAcousticModel << ") (lm= " << setw(10)
			<< std::setprecision(4) << (*it)->fScoreLanguageModel << ") (ip= " << setw(4) << (*it)->fInsertionPenalty 
				<< ") (conf= " << setw(4) << (*it)->fScoreConfidence << ")" << endl;	
		}
	}
}

// write the best path to a file
void BestPath::write(ostream &os, const char *strAudioFile, const char *strSpeakerId, float fOffset, 
	bool bOutputSentenceDelimiters, bool bOutputFillers, bool bOutputConfidenceValues) {
	
	// print the lexical units one per line with starting and ending times
	for( list<BestPathElement*>::iterator it = m_lBestPathElements.begin() ; it != m_lBestPathElements.end() ; ++it ) {
	
		bool bSentenceDelimiter = m_lexiconManager->isSentenceDelimiter((*it)->lexUnit);
		
		// skip sentence delimiters if necessary
		if ((bOutputSentenceDelimiters == false) && (bSentenceDelimiter)) {
			continue;	
		}
		// skip fillers if necessary
		if ((bOutputFillers == false) && ((*it)->lexUnit->iType == LEX_UNIT_TYPE_FILLER)) {
			continue;
		}
	
		float fStart = 0;
		float fDuration = 0;
		if (!bSentenceDelimiter) {
			fStart = fOffset + ((*it)->iFrameStart + m_iFrameStart) * 0.01f; 
			fDuration = ((*it)->iFrameEnd - (*it)->iFrameStart + 1) * 0.01f;
		}
		
		ostringstream oss;
		oss << strAudioFile << " " << setw(6) << strSpeakerId << " " << setw(10) << fStart << " " << setw(10) << fDuration << " " << setw(22) << 	m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit) << " " << endl;
		
		IOBase::writeString(os,oss);
		os.flush();
	}
}

// write the best path to a file (trn format)
void BestPath::write(ostream &os, const char *strUtteranceId, 
	bool bOutputSentenceDelimiters, bool bOutputFillers, bool bOutputConfidenceValues) {
	
	ostringstream oss;
	
	// print the lexical units one per line with starting and ending times
	for( list<BestPathElement*>::iterator it = m_lBestPathElements.begin() ; it != m_lBestPathElements.end() ; ++it) {
		
		bool bSentenceDelimiter = m_lexiconManager->isSentenceDelimiter((*it)->lexUnit);
		
		// skip sentence delimiters if necessary
		if ((bOutputSentenceDelimiters == false) && (bSentenceDelimiter)) {
			continue;	
		}
		// skip fillers if necessary
		if ((bOutputFillers == false) && ((*it)->lexUnit->iType == LEX_UNIT_TYPE_FILLER)) {
			continue;
		}	
		
		// print the lexical unit
		oss << m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit) << " ";
	}
	// print the utterance id
	oss << "(" << strUtteranceId << ")" << endl;
	
	IOBase::writeString(os,oss);
	os.flush();
}


// sets the starting frame to which the word time-alignments refer to (for auto-endpointing)
void BestPath::setStartFrame(int iFrameStart) {

	m_iFrameStart = iFrameStart;
}

// compares two paths for equality 
bool BestPath::compare(BestPath *path) {

	// compare the number of elements
	if (m_lBestPathElements.size() != path->m_lBestPathElements.size()) {
		return false;
	}
	
	// compare all the elements one by one
	LBestPathElement::const_iterator it = m_lBestPathElements.begin();
	LBestPathElement::const_iterator jt = path->m_lBestPathElements.begin();
	for( ; ((it != m_lBestPathElements.end()) && (jt != path->m_lBestPathElements.end())) ; ++it, ++jt) {
		if (memcmp(*it,*jt,sizeof(BestPath)) != 0) {
			return false;
		}
	}
	
	return true;
}

// append a best path to the current one
// note: it is necessary to update the starting and ending times
// note: this function is destructive
void BestPath::append(BestPath *bestPath) {

	// add the elements to the current path
	for(LBestPathElement::iterator it = bestPath->m_lBestPathElements.begin() ; it != bestPath->m_lBestPathElements.end() ; ++it) {
		m_lBestPathElements.push_back(*it);
	}
	
	// update the score
	m_fScore += bestPath->m_fScore;	
}

// return the sequence of lexical units
void BestPath::getLexUnits(VLexUnit &vLexUnit) {

	for(LBestPathElement::iterator it = m_lBestPathElements.begin() ; it != m_lBestPathElements.end() ; ++it) {
		vLexUnit.push_back((*it)->lexUnit);	
	}
}

};	// end-of-namespace

