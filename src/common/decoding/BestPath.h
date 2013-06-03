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


#ifndef BESTPATH_H
#define BESTPATH_H

#include "LexiconManager.h"
#include "IOBase.h"

using namespace std;
#include <iostream>
#include <list>
#include <vector>

namespace Bavieca {

// each of the elements of a hypothesis
typedef struct {
	int iFrameStart;					// start frame aligned to the lexical unit
	int iFrameEnd;						// end frame aligned to the lexical unit
	float fScore;						// accumulated score at the end of the element
	float fScoreAcousticModel;		// acoustic score
	float fScoreLanguageModel;		// language model score
	float fScoreConfidence;			// confidence score	
	LexUnit *lexUnit;					// lexical unit
	float fInsertionPenalty;		// insertion penalty
} BestPathElement;

typedef list<BestPathElement*> LBestPathElement;
typedef vector<BestPathElement*> VBestPathElement;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class BestPath {

	private:
	
		LexiconManager *m_lexiconManager;
		LBestPathElement m_lBestPathElements;
		int m_iFrameStart;
		float m_fScore;										// global path score

		// create a new best path element
		inline BestPathElement *newBestPathElement(int iFrameStart, int iFrameEnd, float fScore, 
			float fScoreAcousticModel, float fScoreLanguageModel,	float fScoreConfidence, 
			LexUnit *lexUnit, float fInsertionPenalty, int iIndexReference = -1) {
			
			BestPathElement *bestPathElement = new BestPathElement;
			bestPathElement->iFrameStart = iFrameStart;
			bestPathElement->iFrameEnd = iFrameEnd;
			bestPathElement->fScore = fScore;
			bestPathElement->fScoreAcousticModel = fScoreAcousticModel;
			bestPathElement->fScoreLanguageModel = fScoreLanguageModel;
			bestPathElement->fScoreConfidence = fScoreConfidence;
			bestPathElement->lexUnit = lexUnit;
			bestPathElement->fInsertionPenalty = fInsertionPenalty;
			
			return bestPathElement;
		}
	
	public:
	
		// constructor
		BestPath(LexiconManager *lexiconManager, float fScore);

		// destructor
		~BestPath();	
		
		// return the number of elements
		inline unsigned int size() {
		
			return (unsigned int)m_lBestPathElements.size();	
		}
		
		// return the path score
		inline float getPathScore() {
		
			return m_fScore;
		}
		
		// return an element
		
		// print a bestPath
		void print(bool bExtended = true);	
		
		// sets the starting frame to which the word time-alignments refer to (for auto-endpointing)
		void setStartFrame(int iFrameStart);
			
		// set the global score
		inline void setScore(float fScore) {
		
			m_fScore = fScore;
		}
		
		// write the best path to a file
		void write(ostream &os, const char *strAudioFile, const char *strSpeakerId, float fOffset,
			bool bOutputSentenceDelimiters = false, bool bOutputFillers = false, bool bOutputConfidenceValues = false);	
		
		// write the best path to a file (trn format)
		void write(ostream &os, const char *strUtteranceId, 
			bool bOutputSentenceDelimiters = false, bool bOutputFillers = false, bool bOutputConfidenceValues = false);

		// convert the best path to a string of characters
		char *toText(bool bOutputSentenceDelimiters = false, bool bOutputFillers = false, bool bOutputConfidenceValues = false);
		
		// return a reference to the list containing the best path
		inline LBestPathElement* getBestPathElements() {
		
			return &m_lBestPathElements;
		}	
		
		// add a new element to the tail of the best path
		inline BestPathElement *newElementBack(int iFrameStart, int iFrameEnd, float fScore, float fScoreAcousticModel, float fScoreLanguageModel,	float fScoreConfidence, LexUnit *lexUnit, float fInsertionPenalty, int iIndexReference = -1) {	
			
			BestPathElement *bestPathElement = newBestPathElement(iFrameStart,iFrameEnd,fScore,fScoreAcousticModel,fScoreLanguageModel,
				fScoreConfidence,lexUnit,fInsertionPenalty,iIndexReference);
			
			m_lBestPathElements.push_back(bestPathElement);
			
			return bestPathElement;
		}

		// add a new element to the head of the best path
		inline BestPathElement *newElementFront(int iFrameStart, int iFrameEnd, float fScore, float fScoreAcousticModel, float fScoreLanguageModel,	float fScoreConfidence, LexUnit *lexUnit, float fInsertionPenalty, int iIndexReference = -1) {	
			
			BestPathElement *bestPathElement = newBestPathElement(iFrameStart,iFrameEnd,fScore,fScoreAcousticModel,fScoreLanguageModel,
				fScoreConfidence,lexUnit,fInsertionPenalty,iIndexReference);
			
			m_lBestPathElements.push_front(bestPathElement);
			
			return bestPathElement;
		}
		
		// compares two paths for equality 
		bool compare(BestPath *path);
		
		// append a best path to the current one
		// note: it is necessary to update the starting and ending times
		// note: this function is destructive
		void append(BestPath *bestPath);
		
		// return the sequence of lexical units
		void getLexUnits(VLexUnit &vLexUnit);
};

};	// end-of-namespace

#endif
