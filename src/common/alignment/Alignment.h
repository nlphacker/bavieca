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


#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <stdio.h>
#include <string.h>

using namespace std;

#include <vector>
#include <deque>
#include <fstream>
#include <iostream>

#include "Global.h"
#include "LexiconManager.h"

namespace Bavieca {

// state occupation
typedef struct {
	int iHMMState;
	double dOccupation;
} StateOcc;

typedef vector<StateOcc*> VStateOcc;

// frame alignment (a feature frame can be aligned to 1 or multiple HMM-states)
typedef VStateOcc FrameAlignment;	

typedef deque<FrameAlignment*> VFrameAlignment;

// word-level alignment
typedef struct {
	int iFrameBegin;
	int iFrameEnd;
	int iLexUnitPron;				// lexical unit with pronunciation
	int iIndex;						// index of the lexical unit in the lexicon file
} WordAlignment;

typedef deque<WordAlignment*> VWordAlignment;

#ifndef PHONE_ALIGNMENT
#define PHONE_ALIGNMENT

typedef struct {
	unsigned char iPhone;									// phonetic symbol relative to the phonetic symbol set
	unsigned char iPosition;								// position of the phone within the word
	LexUnit *lexUnit;											// lexical unit (pronunciation sensitive)
	int iStateBegin[NUMBER_HMM_STATES];					// first time frame aligned to each of the HMM-states in the phone
	int iStateEnd[NUMBER_HMM_STATES];					// first time frame aligned to each of the HMM-states in the phone
	float fLikelihoodState[NUMBER_HMM_STATES];		// log-likelihood of the individual states of the phone
	float fLikelihood;										// log-likelihood of the phone
} PhoneAlignment;

typedef deque<PhoneAlignment*> VPhoneAlignment;

#endif

// alignment type
#define ALIGNMENT_TYPE_FORWARD_BACKWARD		0			// soft aassignment of time frames to HMM-states
#define ALIGNMENT_TYPE_VITERBI					1			// hard alignment (frame by frame)

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class Alignment {

	private:
	
		unsigned char m_iType;
		VFrameAlignment m_vFrameAlignment;
		VWordAlignment m_vWordAlignment;
		bool m_bWordLevelAlignment;				// whether there is word-level alignment information available

	public:
    
		// constructor
		Alignment(unsigned char iType);

		// destructor
		~Alignment();
		
		// add a lex-unit alignment
		void addLexUnitAlignmentFront(int iFrameBegin, int iFrameEnd, LexUnit *lexUnit);
		
		// add a lex-unit alignment
		void addLexUnitAlignmentBack(int iFrameBegin, int iFrameEnd, LexUnit *lexUnit);
		
		// store to disk
		void store(const char *strFile);
		
		// load from disk
		static Alignment *load(const char *strFile, LexiconManager *lexiconManager);
		
		// allocator
		static inline StateOcc *newStateOcc(int iHMMState, double dOccupation) {
		
			StateOcc *stateOcc = new StateOcc;
			stateOcc->iHMMState = iHMMState;
			stateOcc->dOccupation = dOccupation;
		
			return stateOcc;
		}
		
		// allocator
		static inline WordAlignment *newWordAlignment(int iFrameBegin, int iFrameEnd, LexUnit *lexUnit) {
		
			WordAlignment *wordAlignment = new WordAlignment;
			wordAlignment->iFrameBegin = iFrameBegin;
			wordAlignment->iFrameEnd = iFrameEnd;
			if (lexUnit) {
				wordAlignment->iLexUnitPron = lexUnit->iLexUnitPron;
				wordAlignment->iIndex = lexUnit->iIndex;
			} else {
				wordAlignment->iLexUnitPron = -1;
				wordAlignment->iIndex = -1;
			}
		
			return wordAlignment;
		}
		
		// return the alignment type
		inline unsigned char getType() {
		
			return m_iType;
		}
		
		// return the number of frames
		inline unsigned int getFrames() {
		
			return (int)m_vFrameAlignment.size();
		}
		
		// return a frame alignment
		inline FrameAlignment *getFrameAlignment(int t) {
		
			return m_vFrameAlignment[t];
		}
		
		// return a word alignment
		inline VWordAlignment *getWordAlignment() {
		
			return &m_vWordAlignment;
		}
		

		// add frame alignment
		inline void addFrameAlignmentBack(VStateOcc *vStateOcc) {
		
			assert(vStateOcc->empty() == false);
			m_vFrameAlignment.push_back(vStateOcc);
		}
		
		// add frame alignment
		inline void addFrameAlignmentFront(VStateOcc *vStateOcc) {
		
			assert(vStateOcc->empty() == false);
			m_vFrameAlignment.push_front(vStateOcc);
		}
		
		// print
		void print(LexiconManager *lexiconManager = NULL);
		
		// convert to phone-alignment
		VPhoneAlignment *getPhoneAlignment(LexiconManager *lexiconManager);
};

};	// end-of-namespace

#endif
