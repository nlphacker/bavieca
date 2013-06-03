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

#include "FileOutput.h"
#include "TextAlignment.h"

namespace Bavieca {

// constructor
TextAlignment::TextAlignment(LexiconManager *lexiconManager)
{
	m_lexiconManager = lexiconManager;
	for(int i=0 ; i < TEXT_ALIGNMENT_EVENTS_NUMBER ; ++i) {
		m_iAlignmentEvents[i] = 0;
	}
	m_iScore = -1;
	m_iWordsReference = 0;
	m_iWordsHypothesis = 0;
}

// destructor
TextAlignment::~TextAlignment()
{
	for(VTAElement::iterator it = m_vTAElement.begin() ; it != m_vTAElement.end(); ++it) {
		delete *it;
	}
	m_vTAElement.clear();
}

// adds a new element to the alignment
void TextAlignment::addElement(int iAlignmentEvent, int iIndexReference, int iLexUnitReference, 
	int iIndexHypothesis, int iLexUnitHypothesis) {
	
	m_iAlignmentEvents[getIndex(iAlignmentEvent)]++;

	TAElement *element = new TAElement;
	element->iAlignmentEvent = iAlignmentEvent;
	element->iIndexReference = iIndexReference;
	element->iLexUnitReference = iLexUnitReference;
	element->iIndexHypothesis = iIndexHypothesis;
	element->iLexUnitHypothesis = iLexUnitHypothesis;
	
	if (iIndexReference != -1) {
		++m_iWordsReference;
	}
	if (iIndexHypothesis != -1) {
		++m_iWordsHypothesis;
	}
	
	m_vTAElement.push_back(element);
}

// print the alignment
void TextAlignment::print(bool bDetailed, ostream &os) {

	os << "# align score:       " << m_iScore << endl;
	os << "# reference words:   " << m_iWordsReference << endl;
	os << "# hypothesis words:  " << m_iWordsHypothesis << endl;
	os << "# correct:           " << m_iAlignmentEvents[0] << endl;
	os << "# substitutions:     " << m_iAlignmentEvents[1] << endl;
	os << "# deletions:         " << m_iAlignmentEvents[2] << endl;
	os << "# insertions:        " << m_iAlignmentEvents[3] << endl;
	// actual alignment
	if (bDetailed) {
		os << "# event         reference            hypothesis" << endl;
		for(VTAElement::iterator it = m_vTAElement.begin() ; it != m_vTAElement.end() ; ++it) {
			switch((*it)->iAlignmentEvent) {
				case TEXT_ALIGNMENT_EVENT_CORRECT: {
					os << "correct:      " << m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference) << " " <<
						m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis) << endl;
					break;
				}
				case TEXT_ALIGNMENT_EVENT_SUBSTITUTION: {
					os << "substitution: " << m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference) << " " <<
						m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis) << endl;
					break;
				}
				case TEXT_ALIGNMENT_EVENT_DELETION: {
					os << "deletion:     " << m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference) << endl;
					break;
				}
				case TEXT_ALIGNMENT_EVENT_INSERTION: {
					os << "insertion:    " << m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis) << endl;
					break;
				}
			}
		}
	}
}

// store the alignment into the given file
void TextAlignment::store(const char *strFile, bool bDetailed) {

	FileOutput file(strFile,false);
	print(bDetailed,file.getStream());
}


// return the alignment as a stringetCorrect
char *TextAlignment::getStrAlignment() {

	ostringstream strAlignment;
	
	int i=0;
	for(VTAElement::iterator it = m_vTAElement.begin() ; it != m_vTAElement.end() ; ++it,++i) {
		if (i) {
			strAlignment << "|";
		}
		switch((*it)->iAlignmentEvent) {
			case TEXT_ALIGNMENT_EVENT_CORRECT: {
				strAlignment << "correct(" << m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference) <<
					"," << m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis) << ")";
				break;
			}
			case TEXT_ALIGNMENT_EVENT_SUBSTITUTION: {
				strAlignment << "substitution(" << m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference) <<
					"," << m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis) << ")";
				break;
			}
			case TEXT_ALIGNMENT_EVENT_DELETION: {
				strAlignment << "deletion(" << m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference) << ")";
				break;
			}
			case TEXT_ALIGNMENT_EVENT_INSERTION: {
				strAlignment << "insertion(" << m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis) << ")";
				break;
			}
		}
	}	
	
	if (strAlignment.str().compare("") != 0) {
		char *str = new char[strAlignment.str().length()+1];
		strcpy(str,strAlignment.str().c_str());
		return str;
	} else {
		return NULL;
	}
}

};	// end-of-namespace


