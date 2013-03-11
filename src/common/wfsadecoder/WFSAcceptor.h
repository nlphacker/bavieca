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


#ifndef WFSACCEPTOR_H
#define WFSACCEPTOR_H

#include <iomanip>

#include "Global.h"
#include "LexiconManager.h"
#include "LogMessage.h"

using namespace std;

#include <vector>
#include <map>

#include <stdio.h>

namespace Bavieca {

// input symbol values
#define EPSILON_TRANSITION		0x80000000		// epsilon transition
#define FAKE_TRANSITION			0x40000000		// fake transition (to keep final state weights) 
#define LEX_UNIT_TRANSITION	0x20000000		// transition containing a lexical unit 

// mask
#define LEX_UNIT_TRANSITION_COMPLEMENT		0xDFFFFFFF			// mask to extract the lexical unit


typedef struct _TransitionX TransitionX;

typedef TransitionX* StateX;			// a state is a pointer to its first transition in the array of transitions

typedef vector<StateX*> VStateX;
typedef map<StateX*,bool> MStateX;

typedef struct _TransitionX {			// transition data
	unsigned int iSymbol;		// symbol (either the index of a mixture of gaussians, a lexical unit or epsilon)
	float fWeight;					// weight
	StateX *state;					// destination state
} TransitionOpt;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class WFSAcceptor {

	private:
	
		unsigned int m_iStates;					// number of states
		unsigned int m_iTransitions;			// number of transitions
		TransitionX *m_transitions;			// transitions
		StateX *m_states;							// states
		StateX *m_stateInitial;					// initial state
		
		// symbols
		unsigned int m_iHMMStates;
		unsigned int m_iLexUnits;
		
		// empty contructor
		WFSAcceptor() {	
		}

	public:	
		
		// constructor
		WFSAcceptor(StateX *states, unsigned int iStates, 
			StateX *stateInitial, TransitionX *transitions, unsigned int iTransitions, 
			unsigned int iHMMStates, unsigned int iLexUnits);

		// destructor
		~WFSAcceptor();
		
		// return the initial state
		inline StateX *getInitialState() {
		
			return m_stateInitial;
		}
		
		// return the transitions
		inline TransitionX *getTransitions(unsigned int *iTransitions) {
		
			*iTransitions = m_iTransitions;
		
			return m_transitions;
		}
		
		// determines wether the sequence of input symbols is accepted by this acceptor
		bool acceptSequence(vector<unsigned int> vInputSymbols);
		
		// return the set of destination states after observing the given input symbol from the given input state
		void getNextStates(StateX *stateSource, unsigned int iSymbol, VStateX &vStateDest);	
		
		// print the acceptor information
		void print() {
		
			float fSizeMB = (m_iStates*sizeof(StateX)+m_iTransitions*sizeof(TransitionX))/(1024.0*1024.0);
		
			BVC_VERB << "--- acceptor -------------------------";
			BVC_VERB << " # states:      " << setw(10) << m_iStates;
			BVC_VERB << " # transitions: " << setw(10) << m_iTransitions;
			BVC_VERB << " size:          " << FLT(10,2) << fSizeMB << " MB";
			BVC_VERB << "--------------------------------------";
		}
		
		// store the acceptor to a file
		void store(LexiconManager *lexiconManager, const char *strFile);		
		
		// load the acceptor from a file
		static WFSAcceptor *load(LexiconManager *lexiconManager, const char *strFile);
		
		// check the transition correcness
		bool checkTransition(TransitionX *transitionX) {
			
			// make sure the transition symbol is correct
			// epsilon:
			if (transitionX->iSymbol == EPSILON_TRANSITION) {	
				return true;
			}
			// fake transition
			if (transitionX->iSymbol == FAKE_TRANSITION) {
				return true;
			}
			// lexical unit
			else if (transitionX->iSymbol & LEX_UNIT_TRANSITION) {
				int iLexUnit = transitionX->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				assert((iLexUnit >= 0) && (iLexUnit < (int)m_iLexUnits));
				if ((iLexUnit < 0) || (iLexUnit >= (int)m_iLexUnits)) {
					return false;
				}
			}
			// HMM-state
			else {
				assert(transitionX->iSymbol < LEX_UNIT_TRANSITION);
				int iHMMState = transitionX->iSymbol;				
				assert((iHMMState >= 0) && (iHMMState < (int)m_iHMMStates));
				if ((iHMMState < 0) && (iHMMState >= (int)m_iHMMStates)) {
					return false;
				}
			}
			
			return true;
		}
};

};	// end-of-namespace

#endif
