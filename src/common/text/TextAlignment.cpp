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
void TextAlignment::print(bool bDetailed, FILE *file) {

	fprintf(file,"# align score:       %6d\n",m_iScore);
	fprintf(file,"# reference words:   %6d\n",m_iWordsReference);
	fprintf(file,"# hypothesis words:  %6d\n",m_iWordsHypothesis);
	fprintf(file,"# correct:           %6d\n",m_iAlignmentEvents[0]);
	fprintf(file,"# substitutions:     %6d\n",m_iAlignmentEvents[1]);
	fprintf(file,"# deletions:         %6d\n",m_iAlignmentEvents[2]);
	fprintf(file,"# insertions:        %6d\n",m_iAlignmentEvents[3]);
	// actual alignment
	if (bDetailed == true) {
		fprintf(file,"# event         reference            hypothesis\n");
		for(VTAElement::iterator it = m_vTAElement.begin() ; it != m_vTAElement.end() ; ++it) {
			switch((*it)->iAlignmentEvent) {
				case TEXT_ALIGNMENT_EVENT_CORRECT: {
					fprintf(file,"correct:        %-20s %-20s\n",m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference),
					m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis));
					break;
				}
				case TEXT_ALIGNMENT_EVENT_SUBSTITUTION: {
					fprintf(file,"substitution:   %-20s %-20s\n",m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference),
					m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis));
					break;
				}
				case TEXT_ALIGNMENT_EVENT_DELETION: {
					fprintf(file,"deletion:       %-20s %-20s\n",m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference),"-");
					break;
				}
				case TEXT_ALIGNMENT_EVENT_INSERTION: {
					fprintf(file,"insertion:      %-20s %-20s\n","-",m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis));
					break;
				}
			}
		}
	}
}

// store the alignment into the given file
bool TextAlignment::store(const char *strFile, bool bDetailed) {

	// create the file
	FILE *file = fopen(strFile,"wt");
	if (file == NULL) {
		return false;
	}
	
	print(bDetailed,file);
	
	// close the file
	if (fclose(file) == EOF) {
		return false;
	}

	return true;
}


// return the alignment as a stringetCorrect
char *TextAlignment::getStrAlignment() {

	string strAlignment = "";
	char strAux[1024];
	
	int i=0;
	for(VTAElement::iterator it = m_vTAElement.begin() ; it != m_vTAElement.end() ; ++it,++i) {
		if (i) {
			strAlignment += "|";
		}
		switch((*it)->iAlignmentEvent) {
			case TEXT_ALIGNMENT_EVENT_CORRECT: {
				sprintf(strAux,"correct(%s,%s)",m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference),
				m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis));
				break;
			}
			case TEXT_ALIGNMENT_EVENT_SUBSTITUTION: {
				sprintf(strAux,"substitution(%s,%s)",m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference),
				m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis));
				break;
			}
			case TEXT_ALIGNMENT_EVENT_DELETION: {
				sprintf(strAux,"deletion(%s)",m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference));
				break;
			}
			case TEXT_ALIGNMENT_EVENT_INSERTION: {
				sprintf(strAux,"insertion(%s)",m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis));
				break;
			}
		}
		strAlignment += strAux;
	}	
	
	if (strAlignment.compare("") != 0) {
		char *str = new char[strAlignment.length()+1];
		strcpy(str,strAlignment.c_str());
		return str;
	} else {
		return NULL;
	}
}

};	// end-of-namespace


