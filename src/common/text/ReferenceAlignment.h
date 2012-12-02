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


#ifndef REFERENCEALIGNMENT_H
#define REFERENCEALIGNMENT_H

// alignment events
#define ALIGNMENT_EVENT_CORRECT				0
#define ALIGNMENT_EVENT_INSERTION			1
#define ALIGNMENT_EVENT_DELETION 			2
#define ALIGNMENT_EVENT_SUBSTITUTION		3

#define ALIGNMENT_EVENTS_TOTAL				4

using namespace std;

#include <string>
#include <vector>

#include "LexiconManager.h"

namespace Bavieca {

typedef struct {
	int iAlignmentEvent;				// alignment event (correct, ins, del, sub)
	int iIndexReference;				// index within the reference (if applicable)
	int iLexUnitReference;			// lexical unit in the reference (if applicable)
	int iIndexHypothesis;			// index within the hypothesis (if applicable)
	int iLexUnitHypothesis;			// lexical unit in the hypothesis (if applicable)
	int iFrameStart;					// starting frame
	int iFrameEnd;						// ending frame
	float fScoreConfidence;			// confidence score
} ReferenceAlignmentElement;

typedef vector<ReferenceAlignmentElement*> VReferenceAlignmentElement;

/**
	@author root <root@localhost.localdomain>
*/
class ReferenceAlignment {

	private:	
	
		VReferenceAlignmentElement m_vReferenceAlignmentElement;
		LexiconManager *m_lexiconManager;
		int m_iAlignmentEvents[ALIGNMENT_EVENTS_TOTAL];
		
	private:
	
		// write an element to the file
		bool writeElementToFile(FILE *file, ReferenceAlignmentElement *element);

	public:
	
		// constructor
		ReferenceAlignment(LexiconManager *lexiconManager);

		// destructor
		~ReferenceAlignment();
		
		// adds a new element to the alignment
		inline void addElement(int iAlignmentEvent, int iIndexReference, int iLexUnitReference, int iIndexHypothesis, int iLexUnitHypothesis, int iFrameStart, int iFrameEnd, float fScoreConfidence) {
		
			assert(iAlignmentEvent < ALIGNMENT_EVENTS_TOTAL);
			m_iAlignmentEvents[iAlignmentEvent]++;
		
			ReferenceAlignmentElement *element = new ReferenceAlignmentElement;
			element->iAlignmentEvent = iAlignmentEvent;
			element->iIndexReference = iIndexReference;
			element->iLexUnitReference = iLexUnitReference;
			element->iIndexHypothesis = iIndexHypothesis;
			element->iLexUnitHypothesis = iLexUnitHypothesis;
			element->iFrameStart = iFrameStart;
			element->iFrameEnd = iFrameEnd;
			element->fScoreConfidence = fScoreConfidence;
			
			/*printf("%d\n",element->iAlignmentEvent);
			printf("%d\n",element->iIndexReference);
			printf("%d\n",element->iLexUnitReference);
			printf("%d\n",element->iIndexHypothesis);
			printf("%d\n",element->iLexUnitHypothesis);
			printf("%d\n",element->iFrameStart);
			printf("%d\n",element->iFrameEnd);
			printf("%f\n",element->fScoreConfidence);*/
			
			m_vReferenceAlignmentElement.push_back(element);
		}
		
		// print the alignment statistics
		void printStatistics();
		
		// write the alignment information to a file
		bool writeToFile(const char *strFile);
		
		// transform an alignment event to string format
		const char *getStrAlignmentEvent(int iAlignmentEvent);

};

};	// end-of-namespace

#endif
