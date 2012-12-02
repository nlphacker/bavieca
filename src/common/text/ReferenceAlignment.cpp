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


#include "ReferenceAlignment.h"

namespace Bavieca {

// constructor
ReferenceAlignment::ReferenceAlignment(LexiconManager *lexiconManager) {

	m_lexiconManager = lexiconManager;
	
	// initialization
	for(int i=0 ; i < ALIGNMENT_EVENTS_TOTAL ; ++i) {
		m_iAlignmentEvents[i] = 0;
	}
}

// destructor
ReferenceAlignment::~ReferenceAlignment() {

	// remove the elements if any
	if (m_vReferenceAlignmentElement.empty() == false) {
		for(VReferenceAlignmentElement::iterator it = m_vReferenceAlignmentElement.begin() ; it != m_vReferenceAlignmentElement.end() ; ++it) {
			delete *it;
		}
		m_vReferenceAlignmentElement.clear();
	}
}

// write the alignment information to a file
bool ReferenceAlignment::writeToFile(const char *strFile) {

	// open the file for writing
	FILE *file = fopen(strFile,"w");
	if (file == NULL) {
		return false;
	}


	// write the elements one by one
	for(VReferenceAlignmentElement::iterator it = m_vReferenceAlignmentElement.begin() ; it != m_vReferenceAlignmentElement.end() ; ++it) {
		
		assert(*it != NULL);
		writeElementToFile(file,*it);
	}	
	
	fclose(file);

	return true;
}

// write an element to the file
bool ReferenceAlignment::writeElementToFile(FILE *file, ReferenceAlignmentElement *element) {

	string strAlignmentEvent = getStrAlignmentEvent(element->iAlignmentEvent);	
	
	fprintf(file,"%16s %6d %6d %20s %20s %10d %10d %10.4f\n",strAlignmentEvent.c_str(),element->iIndexReference,element->iIndexHypothesis,m_lexiconManager->getStrLexUnit(element->iLexUnitReference),m_lexiconManager->getStrLexUnit(element->iLexUnitHypothesis),element->iFrameStart,element->iFrameEnd,element->fScoreConfidence);
	
	return true;
}


// transform an alignment event to string format
const char *ReferenceAlignment::getStrAlignmentEvent(int iAlignmentEvent) {

	switch(iAlignmentEvent) {
		case ALIGNMENT_EVENT_CORRECT: {
			return "correct";
			break;
		}
		case ALIGNMENT_EVENT_INSERTION: {
			return "insertion";
			break;
		}
		case ALIGNMENT_EVENT_DELETION: {
			return "deletion";
			break;
		}
		case ALIGNMENT_EVENT_SUBSTITUTION: {
			return "substitution";
			break;
		}
		default: {
			assert(0);
		}
	}

	return NULL;
}

// print the alignment statistics
void ReferenceAlignment::printStatistics() {

	printf("--------------------------------------------------------\n");
	printf("Alignment statistics:\n");
	printf(" Correct:       %d\n",m_iAlignmentEvents[ALIGNMENT_EVENT_CORRECT]);
	printf(" Insertions:    %d\n",m_iAlignmentEvents[ALIGNMENT_EVENT_INSERTION]);
	printf(" Deletions:     %d\n",m_iAlignmentEvents[ALIGNMENT_EVENT_DELETION]);
	printf(" Substitutions: %d\n",m_iAlignmentEvents[ALIGNMENT_EVENT_SUBSTITUTION]);
	printf("--------------------------------------------------------\n");
}

};	// end-of-namespace

