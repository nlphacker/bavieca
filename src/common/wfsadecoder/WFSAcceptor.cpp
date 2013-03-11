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


#include "WFSAcceptor.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"

namespace Bavieca {

WFSAcceptor::WFSAcceptor(StateX *states, unsigned int iStates, StateX *stateInitial, 
	TransitionX *transitions, unsigned int iTransitions, 
	unsigned int iHMMStates, unsigned int iLexUnits)
{
	m_states = states;
	m_iStates = iStates;
	m_stateInitial = stateInitial;
	m_transitions = transitions;
	m_iTransitions = iTransitions;
	m_iHMMStates = iHMMStates;
	m_iLexUnits = iLexUnits;
	
}

WFSAcceptor::~WFSAcceptor()
{
	// deallocate memory
	delete [] m_states;
	delete [] m_transitions;
}

// determines wether the sequence of input symbols is accepted by this acceptor
bool WFSAcceptor::acceptSequence(vector<unsigned int> vInputSymbols) {

	VStateX *vStates1 = new VStateX;
	VStateX *vStates2 = new VStateX;
	
	vStates1->push_back(m_stateInitial);
	
	for(unsigned int i = 0 ; i < vInputSymbols.size() ; ++i) {
		for(VStateX::iterator it = vStates1->begin() ; it != vStates1->end() ; ++it) {
			getNextStates(*it,vInputSymbols[i],*vStates2);
		}
		vStates1->clear();
		VStateX *aux = vStates1;
		vStates1 = vStates2;
		vStates2 = aux;
		if (vStates1->empty()) {
			cout << vInputSymbols[i] << "not seen!" << endl;
			return false;
		} else {
			cout << vInputSymbols[i] << " seen" << endl;
		}
	}

	return true;
}

// return the set of destination states after observing the given input symbol from the given input state
void WFSAcceptor::getNextStates(StateX *stateSource, unsigned int iSymbol, VStateX &vStateDest) {

	TransitionX *transition = *stateSource;
	TransitionX *transitionEnd = *(stateSource+1);
	while(transition != transitionEnd) {
	
		// epsilon transition
		if (transition->iSymbol & EPSILON_TRANSITION) {
			getNextStates(transition->state,iSymbol,vStateDest);
		} 
		// lex-unit transition
		else if (transition->iSymbol & LEX_UNIT_TRANSITION) {
			getNextStates(transition->state,iSymbol,vStateDest);
		} 
		// leaf-transition
		else if ((transition->iSymbol & FAKE_TRANSITION) == 0) {
			if (transition->iSymbol == iSymbol) {
				if (transition->state != NULL) {
					vStateDest.push_back(transition->state);
				}
			}
		}
		
		++transition;
	}
}

// store the acceptor to a file
void WFSAcceptor::store(LexiconManager *lexiconManager, const char *strFile) {

	FileOutput file(strFile,true);
	file.open();
	
	IOBase::write(file.getStream(),m_iHMMStates);
	IOBase::write(file.getStream(),m_iLexUnits);
	IOBase::write(file.getStream(),m_iStates);
	IOBase::write(file.getStream(),m_iTransitions);
	
	// store the transitions along with the state information
	for(unsigned int i = 0 ; i < m_iStates ; ++i) {
		// # of transitions
		unsigned int iTransitions = m_states[i+1]-m_states[i];
		IOBase::write(file.getStream(),iTransitions);
		for(TransitionX *transition = m_states[i]; transition != m_states[i+1] ; ++transition) {
			
			if (transition->iSymbol & LEX_UNIT_TRANSITION) {
				LexUnit *lexUnit = lexiconManager->getLexUnitPron(transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
				IOBase::write(file.getStream(),lexUnit->iIndex|LEX_UNIT_TRANSITION);
			} else {
				IOBase::write(file.getStream(),transition->iSymbol);
			}	
			IOBase::write(file.getStream(),transition->fWeight);
			unsigned int iStateDest = transition->state-m_states;
			IOBase::write(file.getStream(),iStateDest);
		}
	}
	
	file.close();
}

// load the acceptor from a file
WFSAcceptor *WFSAcceptor::load(LexiconManager *lexiconManager, const char *strFile) {

	WFSAcceptor *m_wfsAcceptor = new WFSAcceptor();
	
	// open the file for reading
	FileInput file(strFile,true);
	file.open();
	
	IOBase::read(file.getStream(),&m_wfsAcceptor->m_iHMMStates);
	IOBase::read(file.getStream(),&m_wfsAcceptor->m_iLexUnits);
	IOBase::read(file.getStream(),&m_wfsAcceptor->m_iStates);
	IOBase::read(file.getStream(),&m_wfsAcceptor->m_iTransitions);	

	assert((m_wfsAcceptor->m_iHMMStates > 0) && (m_wfsAcceptor->m_iLexUnits > 0));	
	
	// allocate memory for the states and transitions
	m_wfsAcceptor->m_states = new StateX[m_wfsAcceptor->m_iStates+1];
	m_wfsAcceptor->m_transitions = new TransitionX[m_wfsAcceptor->m_iTransitions];
	m_wfsAcceptor->m_states[0] = m_wfsAcceptor->m_transitions;
	
	// load the transitions along with the state information
	unsigned int iTransitionOffset = 0;
	for(unsigned int i = 0 ; i < m_wfsAcceptor->m_iStates ; ++i) {
		// # of transitions
		unsigned int iTransitions = 0;
		IOBase::read(file.getStream(),&iTransitions);	
		m_wfsAcceptor->m_states[i+1] = m_wfsAcceptor->m_transitions+iTransitionOffset+iTransitions;
		for(unsigned int j = 0 ; j < iTransitions ; ++j) {
			// transition symbol
			unsigned int iSymbol;
			IOBase::read(file.getStream(),&iSymbol);
			if (iSymbol & LEX_UNIT_TRANSITION) {
				LexUnit *lexUnit = lexiconManager->getLexUnitByIndex(iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT);
				m_wfsAcceptor->m_transitions[iTransitionOffset].iSymbol = lexUnit->iLexUnitPron|LEX_UNIT_TRANSITION;
			} else {
				m_wfsAcceptor->m_transitions[iTransitionOffset].iSymbol = iSymbol;
			}	
			IOBase::read(file.getStream(),&m_wfsAcceptor->m_transitions[iTransitionOffset].fWeight);	
			// transition destination state
			unsigned int iStateDest = 0;
			IOBase::read(file.getStream(),&iStateDest);	
			m_wfsAcceptor->m_transitions[iTransitionOffset].state = m_wfsAcceptor->m_states+iStateDest;
			// check the transition
			if (m_wfsAcceptor->checkTransition(&m_wfsAcceptor->m_transitions[iTransitionOffset]) == false) {
				return NULL;
			}
			if (m_wfsAcceptor->m_transitions[iTransitionOffset].iSymbol & FAKE_TRANSITION) {
				m_wfsAcceptor->m_transitions[iTransitionOffset].state = NULL;
			}
			++iTransitionOffset;	
		}
	}
	m_wfsAcceptor->m_states[m_wfsAcceptor->m_iStates] = m_wfsAcceptor->m_transitions+m_wfsAcceptor->m_iTransitions;
	m_wfsAcceptor->m_stateInitial = m_wfsAcceptor->m_states;
	
	// sanity check
	for(unsigned int i=0 ; i<m_wfsAcceptor->m_iStates ; ++i) {
		assert(m_wfsAcceptor->m_states[i] < m_wfsAcceptor->m_states[i+1]);
	}

	// close the file
	file.close();
	
	return m_wfsAcceptor;
}

};	// end-of-namespace
