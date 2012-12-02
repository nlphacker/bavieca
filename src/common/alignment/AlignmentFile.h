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


#ifndef ALIGNMENTFILE_H
#define ALIGNMENTFILE_H

using namespace std;

#include <string>
#include <vector>
#include <deque>
#include <list>
#include <sstream>

#include "LexiconManager.h"
#include "HMMState.h"

namespace Bavieca {

class Alignment;
class HMMManager;
class PhoneSet;

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

typedef struct {
	LexUnit *lexUnit;			// lexical unit (pronunciation sensitive)
	int iBegin;					// first frame aligned to the lexical unit
	int iEnd;					// last frame aligned to the lexical unit
	float fLikelihood;		// acoustic log-likelihood
} LexUnitAlignment;

typedef deque<LexUnitAlignment*> VLexUnitAlignment;
typedef list<LexUnitAlignment*> LLexUnitAlignment;

/**
	@author root <root@localhost.localdomain>
*/
class AlignmentFile {

	private:
	
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		string m_strFile;

	public:
	
		// constructor
		AlignmentFile(PhoneSet *phoneSet, LexiconManager *lexiconManager = NULL);

		// destructor
		~AlignmentFile();
		
		// store an alignment into a file
		void store(VPhoneAlignment &vPhoneAlignment, const char *strAlignmentFile);
		
		// load an alignment from a file
		VPhoneAlignment *load(const char *strAlignmentFile);	
			
		// prints the phone alignment
		void print(VPhoneAlignment &vPhoneAlignment);
		
		// prints the lexical unit alignment
		void print(VLexUnitAlignment &vLexUnitAlignment);	
		
		// destroy a phone alignment
		static void destroyPhoneAlignment(VPhoneAlignment *vPhoneAlignment) {
		
			for(VPhoneAlignment::iterator it = vPhoneAlignment->begin() ; it != vPhoneAlignment->end() ; ++it) {
				delete *it;
			}	
			delete vPhoneAlignment;
		}
		
		// destroy a lexUnit alignment
		static void destroyLexUnitAlignment(VLexUnitAlignment *vLexUnitAlignment) {
		
			for(VLexUnitAlignment::iterator it = vLexUnitAlignment->begin() ; it != vLexUnitAlignment->end() ; ++it) {
				delete *it;
			}	
			delete vLexUnitAlignment;
		}
		
		// convert to a object of the class Alignment
		static Alignment *toAlignment(PhoneSet *phoneSet, HMMManager *hmmManager, VPhoneAlignment *vPhoneAlignment);
		
		// return a lexical unit alignment from a phone alignment
		static VLexUnitAlignment *getVLexUnitAlignment(VPhoneAlignment &vPhoneAlignment);
		
		// return a given lexical unit alignment from a phone alignment
		static LexUnitAlignment *getLexUnitAlignment(VPhoneAlignment &vPhoneAlignment, unsigned int iIndex);	
		
		// prints the alignment
		static void print(PhoneSet *phoneSet, LexiconManager *lexiconManager, VPhoneAlignment &vPhoneAlignment);
		
		// prints the lexical unit alignment
		static void print(PhoneSet *phoneSet, LexiconManager *lexiconManager, VLexUnitAlignment &vLexUnitAlignment);	
};

};	// end-of-namespace

#endif
