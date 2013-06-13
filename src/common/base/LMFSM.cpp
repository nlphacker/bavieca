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


#include "LMFSM.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "LMManager.h"
#include "LogMessage.h"
#include "TimeUtils.h"

namespace Bavieca {

//constructor
LMFSM::LMFSM(LexiconManager *lexiconManager) {

	m_lexiconManager = lexiconManager;
	m_bLoaded = false;
	m_iNGrams = NULL;
	m_iNGramOrder = -1;

	m_iLMStateInitial = -1;
	m_iLMStateFinal = -1;
	m_states = NULL;
	m_arcs = NULL;	
	m_iStates = 0;
	m_iArcs = 0;
	m_iArcsStandard = 0;
	m_iArcsBackoff = 0;		
}

// destructor
LMFSM::~LMFSM() {

	if (m_bLoaded) {

		if (m_states) {
			delete [] m_states;
			m_states = NULL;
		}
		if (m_arcs) {
			delete [] m_arcs;
			m_arcs = NULL;
		}
	
		m_iLMStateInitial = -1;
		m_iLMStateFinal = -1;
		m_iStates = 0;
		m_iArcs = 0;
		m_iArcsStandard = 0;
		m_iArcsBackoff = 0;
	
		delete [] m_iNGrams;
		m_bLoaded = false;
	}	
}

// build the FSM from the n-grams
void LMFSM::build(LMARPA *lmARPA) {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	m_iNGramOrder = lmARPA->getNGramOrder();
	//m_iNGramOrder = 0;
	m_iNGrams = new int[m_iNGramOrder+1];
	for(int i=0; i<m_iNGramOrder ; ++i) {
		m_iNGrams[i] = lmARPA->getNGrams(i);
	}
	unsigned int *iIgnoredNGrams = new unsigned int[m_iNGramOrder+1];
	for(int i=0; i<=m_iNGramOrder ; ++i) {
		iIgnoredNGrams[i] = 0;
	}
	
	int iLexUnitUnknown = m_lexiconManager->m_lexUnitUnknown->iLexUnit;
	int iLexUnitBegSentence = m_lexiconManager->m_lexUnitBegSentence->iLexUnit;
	int iLexUnitEndSentence = m_lexiconManager->m_lexUnitEndSentence->iLexUnit;
	NGram *ngramFinal = lmARPA->getNGram(&iLexUnitEndSentence,1);
	assert(ngramFinal);
	
	LMStateTemp *states = NULL;
	LMStateTemp *stateInitial = NULL;
	LMStateTemp *stateFinal = NULL;
	
	// zerogram (special case)
	if (m_iNGramOrder == LM_NGRAM_ZEROGRAM) {
		states = new LMStateTemp[1];
		states[0].iState = 0;		
		m_iStates = 1;
		stateInitial = stateFinal = states;
	
		// compacting
		compact(states,m_iStates,m_iArcs,stateInitial,stateFinal);
			
		BVC_VERB << "-- FSM ------------------------------------";
		BVC_VERB << " zerogram (uniform distribution)";
		BVC_VERB << "-------------------------------------------";
		
		// clean-up
		delete [] states;
		delete [] iIgnoredNGrams;
		m_bLoaded = true;
		return;
	}
	
	// allocate memory to store the states (for example, for trigram: backoff + #unigrams + #bigrams)
	int iStatesMax = 1;
	for(int i=0 ; i < m_iNGramOrder ; ++i) {
		iStatesMax += lmARPA->getNGrams(i);
	}
	states = new LMStateTemp[iStatesMax];
	for(int i=0 ; i < iStatesMax ; ++i) {
		states[i].iState = -1;
	}
	
	// create states and arcs
	int iStateID = 0;
	
	// maps NGram <-> state in the FSM
	MNGramState mNGramState;
	
	LMStateTemp *statesPrev = states;		// points to single back-off state
	statesPrev->iState = iStateID++;			// zerogram backoff state
	LMStateTemp *stateDest = states+1;		// points to unigram states
	
	int iNGramsPrev = -1;	
	NGram *ngramsPrev = lmARPA->getNGrams(0,&iNGramsPrev);
	mNGramState.insert(MNGramState::value_type(hashKey(ngramsPrev),statesPrev));
	for(int iOrder = 1 ; iOrder <= m_iNGramOrder ; ++iOrder) {
		
		for(int i=0 ; i < iNGramsPrev ; ++i) {		// there is one FSM state for each previous n-gram
		
			// standard arcs
			//LMStateTemp *statePrev = statesPrev+i;
			NGram *prev = ngramsPrev+i;
			assert(prev->iNGrams >= 0);
			MNGramState::iterator it = mNGramState.find(hashKey(prev));
			if (it == mNGramState.end()) {
				iIgnoredNGrams[iOrder] += prev->iNGrams;
				continue;
			}
			LMStateTemp *statePrev = it->second;
			
			// skip unreachable states (also, no transitions from final state)
			assert(statePrev->iState != -1);
			if (prev->iLexUnit == iLexUnitEndSentence) {
				iIgnoredNGrams[iOrder] += prev->iNGrams;
				stateDest += prev->iNGrams;
				continue;	
			}
			for(int j=0 ; j < prev->iNGrams ; ++j, ++stateDest) {
				int iLexUnit = prev->ngrams[j].iLexUnit;
				if ((iLexUnit == iLexUnitUnknown) /*|| ((iLexUnit == iLexUnitBegSentence) && (iOrder > 1))*/) {
					//lmARPA->print(&prev->ngrams[j]);
					++iIgnoredNGrams[iOrder];
					continue;
				} 
				// connect to final state (</s>)
				else if ((iLexUnit == iLexUnitEndSentence) && (iOrder > 1)) {
					MNGramState::iterator it = mNGramState.find(hashKey(ngramFinal));
					assert(it != mNGramState.end());
					statePrev->lArc.push_back(newArc(iLexUnit,prev->ngrams[j].fProbability,it->second));
					continue;
				}
				if (iOrder < m_iNGramOrder) {
					mNGramState.insert(MNGramState::value_type(hashKey(&prev->ngrams[j]),stateDest));
					stateDest->iState = iStateID++;
					statePrev->lArc.push_back(newArc(iLexUnit,prev->ngrams[j].fProbability,stateDest));	
				} else {
					//lmARPA->print(&prev->ngrams[j]);
					assert(iOrder == m_iNGramOrder);	
					MNGramState::iterator it = mNGramState.find(hashKey(prev,true,(m_iNGramOrder > 1) ? iLexUnit : -1));	
					if (it != mNGramState.end()) {
						statePrev->lArc.push_back(newArc(iLexUnit,prev->ngrams[j].fProbability,it->second));
					} else {
						//lmARPA->print(&prev->ngrams[j]);	
						++iIgnoredNGrams[iOrder];
					}
				}
			}
			// backoff arc
			if (iOrder > 1) {
				MNGramState::iterator it = mNGramState.find(hashKey(prev,true));	
				assert(it != mNGramState.end());
				statePrev->lArc.push_back(newArc(BACKOFF_ARC,prev->fProbabilityBackoff,it->second));
			}
		}	
		if (iOrder == m_iNGramOrder)
			break;
		statesPrev += lmARPA->getNGrams(iOrder-1);
		ngramsPrev = lmARPA->getNGrams(iOrder,&iNGramsPrev);
	}
	
	m_iArcs = m_iArcsStandard+m_iArcsBackoff;
	
	// get initial and final states
	if (m_iNGramOrder <= 1) {
		stateInitial = stateFinal = states;
	} else {	
		// initial state
		MNGramState::iterator it;	
		it = mNGramState.find(hashKey(lmARPA->getNGram(&iLexUnitBegSentence,1)));	
		if (it == mNGramState.end()) {
			BVC_ERROR << "initial state (<s>) not found, make sure there is an unigram for <s>";
		}
		stateInitial = it->second;
		// final state
		it = mNGramState.find(hashKey(lmARPA->getNGram(&iLexUnitEndSentence,1)));	
		if (it == mNGramState.end()) {
			BVC_ERROR << "final state (</s>) not found, make sure there is an unigram for <s>";
		}
		stateFinal = it->second;
	}
	
	m_iStates = iStatesMax;
	
	double dTimeEndBuilding = TimeUtils::getTimeMilliseconds();
	
	// sanity check
	checkConnected(stateInitial,stateFinal,states,iStateID,m_iArcs);
	
	double dTimeEndChecks = TimeUtils::getTimeMilliseconds();
	
	// compacting
	compact(states,m_iStates,m_iArcs,stateInitial,stateFinal);
	
	double dTimeEndCompacting = TimeUtils::getTimeMilliseconds();
	
	// print summary
	BVC_VERB << "-- FSM ------------------------------------";
	BVC_VERB << " # states: " << iStateID;
	BVC_VERB << " # arcs: " << m_iArcs << " (standard: " << m_iArcsStandard << ", backoff: " << m_iArcsBackoff << ")";
	for(int i=0 ; i <= m_iNGramOrder ; ++i) {
		if (iIgnoredNGrams[i] > 0) {
			BVC_VERB << " # ignored " << LMManager::getStrNGram(i) << "s: " << iIgnoredNGrams[i];
		}
	}
	int iBytes = (m_iStates*sizeof(LMState)+m_iArcs*sizeof(LMArc));
	BVC_VERB << " size: " << iBytes << " bytes (" << ((float)iBytes)/(1024.0*1024.0) << " MBs)";
	BVC_VERB << " building time:   " << (dTimeEndBuilding-dTimeBegin)/1000.0 << " seconds";
	BVC_VERB << " check time:      " << (dTimeEndChecks-dTimeEndBuilding)/1000.0 << " seconds";
	BVC_VERB << " compacting time: " << (dTimeEndCompacting-dTimeEndChecks)/1000.0 << " seconds";
	BVC_VERB << " total time: " << (dTimeEndCompacting-dTimeBegin)/1000.0 << " seconds";
	BVC_VERB << "-------------------------------------------";
	
	// clean-up
	for(int i=0 ; i<iStatesMax ; ++i) {
		for(LLMArcTemp::iterator it = states[i].lArc.begin() ; it != states[i].lArc.end() ; ++it) {
			delete *it;
		}
	}	
	delete [] states;		
	delete [] iIgnoredNGrams;
	
	m_bLoaded = true;
}

// perform sanity checks to make sure that all the states/arcs created are connected
void LMFSM::checkConnected(LMStateTemp *stateInitial, LMStateTemp *stateFinal, 
	LMStateTemp *stateBackoffZerogram, int iStates, int iArcs) {
	
	bool *bStateProcessed = new bool[iStates];
	for(int i=0 ; i<iStates ; ++i) {
		bStateProcessed[i] = false;
	}
	int iStatesSeen = 0;
	int iArcsSeen = 0;
	LLMStateTemp lState;
	lState.push_back(stateInitial);
	bStateProcessed[stateInitial->iState] = true;
	while(lState.empty() == false) {
		
		// get the next state to process
		LMStateTemp *stateFrom = lState.front();
		lState.pop_front();
		
		bool bBackoff = false;
		for(LLMArcTemp::iterator it = stateFrom->lArc.begin() ; it != stateFrom->lArc.end() ; ++it) {
		
			if ((*it)->iLexUnit == BACKOFF_ARC) { bBackoff = true; }
		
			// insert G destination states into the queue (if not already processed)
			if (bStateProcessed[(*it)->stateDest->iState] == false) {
				lState.push_back((*it)->stateDest);
				bStateProcessed[(*it)->stateDest->iState] = true;
			}
			// make sure the transition has a valid symbol
			if ((*it)->iLexUnit != BACKOFF_ARC) {
			
				if (!m_lexiconManager->isStandard((*it)->iLexUnit) && 
					((*it)->iLexUnit != m_lexiconManager->m_lexUnitEndSentence->iLexUnit) &&
					((stateFrom != stateBackoffZerogram) || ((*it)->stateDest != stateInitial) || ((*it)->iLexUnit != m_lexiconManager->m_lexUnitBegSentence->iLexUnit))) {
					BVC_WARNING << "unexpected lexical unit: " << m_lexiconManager->getStrLexUnit((*it)->iLexUnit) << " !!";
				}
			}
			
			++iArcsSeen;
		}
		if ((bBackoff == false) && (stateFrom != stateFinal) && (stateFrom != stateBackoffZerogram)) {
			BVC_ERROR << "backoff-arc expected but not found, lm is not well formed";
		}
		
		++iStatesSeen;	
	}
	for(int i=0 ; i<iStates ; ++i) {
		if (bStateProcessed[i] == false) {
			BVC_WARNING << "lm-state not seen: " << i;
		}
	}	
	assert(iStatesSeen == iStates);
	assert(iArcsSeen == iArcs);
	delete [] bStateProcessed;	
}

// compact the FSM to use less memory and speed-up lookups (better locality)
// note: most of the cpu time goes to sorting arcs by lexical unit
void LMFSM::compact(LMStateTemp *states, int iStates, int iArcs, LMStateTemp *stateInitial, LMStateTemp *stateFinal) {
	
	// allocate memory
	m_iStates = iStates;
	m_iArcs = iArcs;
	m_states = new LMState[m_iStates+1];		// extra state is needed to mark the ending arc
	m_arcs = new LMArc[m_iArcs];
	
	// keep already created nodes (indexed by the corresponding temporal node id)
	int *iStateCreated = new int[m_iStates];
	int *iFirstArc = new int[m_iStates];
	for(int i=0 ; i<m_iStates ; ++i) {
		iStateCreated[i] = iFirstArc[i] = -1;
	}

	// fill the structures
	int iState = 0;
	int iArc = 0;
	int iStateActual = 0;
	for(int iStateTemp = 0 ; iStateTemp < m_iStates ; ++iStateTemp) {
	
		if (states[iStateTemp].iState == -1) {
			continue;
		}	
	
		// sort arcs by lexical unit so we can use a binary search for arc lookups
		states[iStateTemp].lArc.sort(LMFSM::compareArcs);
	
		LMState *stateAux;
		int iArcBase = -1;
		if (iStateCreated[iStateActual] == -1) {
			stateAux = &m_states[iState];
			iStateCreated[iStateActual] = iState++;
			iArcBase = iArc;
			iArc += (int)states[iStateTemp].lArc.size();
		} else {		
			stateAux = &m_states[iStateCreated[iStateActual]]; 
			iArcBase = iFirstArc[iStateCreated[iStateActual]];
		}
		stateAux->iArcBase = iArcBase;
		// store the outgoing arcs
		int j=0;
		for(LLMArcTemp::iterator jt = states[iStateTemp].lArc.begin() ; jt != states[iStateTemp].lArc.end() ; ++jt,++j) {
			if (iStateCreated[(*jt)->stateDest->iState] == -1) {
				iStateCreated[(*jt)->stateDest->iState] = iState;
				iFirstArc[iState] = iArc;
				iArc += (int)(*jt)->stateDest->lArc.size();
				iState++;
			}
			m_arcs[iArcBase+j].iLexUnit = (*jt)->iLexUnit;
			m_arcs[iArcBase+j].fScore = (*jt)->fScore;
			m_arcs[iArcBase+j].iStateDest = iStateCreated[(*jt)->stateDest->iState];
		}
		++iStateActual;
	}
	
	m_iStates = iStateActual;
		
	assert(iArc == m_iArcs);
	assert(iState == m_iStates);
	m_states[m_iStates].iArcBase = m_iArcs;
	
	m_iLMStateInitial = iStateCreated[stateInitial->iState];
	m_iLMStateFinal = iStateCreated[stateFinal->iState];	

	// sanity check
	for(int i = 0 ; i < m_iStates ; ++i) {	
		if (states[i].iState == -1) {
			continue;
		}	
		assert(m_states[i].iArcBase <= m_states[i+1].iArcBase);
	}	
	
	delete [] iStateCreated;
	delete [] iFirstArc;
}

// store to disk
void LMFSM::store(const char *strFile) {

	assert(m_bLoaded == true);
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	FileOutput file(strFile,true);
	file.open();
	
	// store the actual FSM	
	IOBase::write(file.getStream(),m_iStates);
	IOBase::write(file.getStream(),m_iArcs);
	
	// store the transitions along with the state information
	for(int i = 0 ; i < m_iStates ; ++i) {
		// # of arcs
		int iArcs = m_states[i+1].iArcBase-m_states[i].iArcBase;
		assert((iArcs > 0) || (i == m_iLMStateFinal));
		IOBase::write(file.getStream(),iArcs);
		// actual arcs
		for(LMArc *arc = m_arcs+m_states[i].iArcBase ; arc != m_arcs+m_states[i+1].iArcBase ; ++arc) {	
			IOBase::write(file.getStream(),arc->iLexUnit);
			IOBase::write(file.getStream(),arc->fScore);
			IOBase::write(file.getStream(),arc->iStateDest);
		}
	}
	
	file.close();
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	BVC_VERB << "language model FSM store time: " << (dTimeEnd-dTimeBegin)/1000 << " seconds";	
}

// load from disk
void LMFSM::load(const char *strFile) {

	assert(m_bLoaded == false);
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// open the file for reading
	FileInput file(strFile,true);
	file.open();
	
	IOBase::read(file.getStream(),&m_iStates);
	IOBase::read(file.getStream(),&m_iArcs);

	// minimum is 1 state and zero arcs (zerogram)
	if ((m_iStates <= 0) || (m_iArcs < 0)) {
		BVC_ERROR << "wrong number of states/arcs: states = " << m_iStates << ", arcs = " << m_iArcs;	
	}
	
	// allocate memory for the states and arcs
	m_states = new LMState[m_iStates+1];
	m_arcs = new LMArc[m_iArcs];
	m_states->iArcBase = 0;
	
	// load the transitions along with the state information
	unsigned int iArcOffset = 0;
	for(int i = 0 ; i < m_iStates ; ++i) {
		// # of arcs
		unsigned int iArcs = 0;
		IOBase::read(file.getStream(),&iArcs);	
		if (iArcs == 0) {
			m_iLMStateFinal = i;
		}
		m_states[i+1].iArcBase = iArcOffset+iArcs;
		for(unsigned int j = 0 ; j < iArcs ; ++j) {
			IOBase::read(file.getStream(),&m_arcs[iArcOffset].iLexUnit);
			IOBase::read(file.getStream(),&m_arcs[iArcOffset].fScore);
			IOBase::read(file.getStream(),&m_arcs[iArcOffset].iStateDest);
			int iStateDest = m_arcs[iArcOffset].iStateDest;
			assert((iStateDest >= 0) && (iStateDest < m_iStates));
			++iArcOffset;	
		}
	}
	m_states[m_iStates].iArcBase = m_iArcs;
	m_iLMStateInitial = 0;
	assert(m_iLMStateFinal != -1);
	
	// sanity check
	for(int i=0 ; i<m_iStates ; ++i) {
		assert((m_states[i].iArcBase < m_states[i+1].iArcBase) || (i == m_iLMStateFinal));
	}

	// close the file
	file.close();
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	BVC_VERB << "language model FSM loading time: " << (dTimeEnd-dTimeBegin)/1000 << " seconds";	

	m_bLoaded = true;
}

// get the initial state
int LMFSM::getInitialState() {

	return m_iLMStateInitial;
}

// update the language model state with the given lexical unit and return the new lm state
int LMFSM::updateLMState(int iLMStatePrev, int iLexUnit, float *fScore) {

	// zerogram (uniform distribution), only one state
	if (m_iNGramOrder == LM_NGRAM_ZEROGRAM) {
		*fScore = log10f(1.0/((float)m_lexiconManager->getVocabularySize()));
		return m_iLMStateInitial;
	}	
	// higher order n-grams
	else {
		assert(m_iNGramOrder >= LM_NGRAM_UNIGRAM);	
		assert((iLMStatePrev >= 0) && (iLMStatePrev < m_iStates));
		LMState *state = &m_states[iLMStatePrev];
		*fScore = 0.0;
		int iPasses = 0;
		
		while(1) {
		
			int iFirst = state->iArcBase;
			int iLast = (state+1)->iArcBase-1;
			int iMiddle;
			while(iFirst <= iLast) {
				iMiddle = (iFirst+iLast)/2;
				if (m_arcs[iMiddle].iLexUnit == iLexUnit) {
					//cout << m_lexiconManager->getStrLexUnit(iLexUnit) << ": " << m_arcs[iMiddle].fScore << endl;
					*fScore += m_arcs[iMiddle].fScore;
					return m_arcs[iMiddle].iStateDest;
				} else if (m_arcs[iMiddle].iLexUnit < iLexUnit) {
					iFirst = iMiddle+1;
				} else {
					iLast = iMiddle-1;
				}
			}	
			
			LMArc *arcBackoff = &m_arcs[(state+1)->iArcBase-1];
			assert(arcBackoff->iLexUnit == BACKOFF_ARC);
			*fScore += arcBackoff->fScore;
			state = &m_states[arcBackoff->iStateDest];
			//cout << "backoff: " << arcBackoff->fScore << endl;
				
			++iPasses;
		}
	}
}

// return the score resulting from moving to the given lm-state to the final state
float LMFSM::toFinalState(int iLMState) {

	// zerogram (uniform distribution)
	if (m_iNGramOrder == LM_NGRAM_ZEROGRAM) {
		return log10f(1.0/((float)m_lexiconManager->getVocabularySize()));
	} 
	// higher order n-grams
	else {
		assert(m_iNGramOrder >= LM_NGRAM_UNIGRAM);
		float fScore = -FLT_MAX;
		int iLMStateAux = updateLMState(iLMState,m_lexiconManager->m_lexUnitEndSentence->iLexUnit,&fScore);
		assert(iLMStateAux == m_iLMStateFinal);
		return fScore;
	} 
}

// return language model scores for all words in the vocabulary for a given LM-state (word history)
void LMFSM::getLMScores(int iLMState, float *fLMScores, int iVocabularySize) {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	LMArc **lmArcBackoff = new LMArc*[m_iNGramOrder-1];
	int iEmpty = m_iNGramOrder-2;
	
	// get the lm-state for each n-gram order (includes actual lm-state plus backoff lm-states)
	LMState *state = m_states+iLMState;
	LMArc *arcBackoff = m_arcs+((state+1)->iArcBase-1);
	while(arcBackoff->iLexUnit == BACKOFF_ARC) {
		assert(iEmpty >= 0);
		lmArcBackoff[iEmpty--] = arcBackoff;
		arcBackoff = m_arcs+(((m_states+arcBackoff->iStateDest)+1)->iArcBase-1);
	}
	
	// initialization (debug purposes)
	for(int i=0 ; i < iVocabularySize; ++i) {
		fLMScores[i] = FLT_MAX;
	}
	
	// (1) unobserved n-grams (backoffs)
	// deal with back-off transitions in ascending order (starting from unigrams, then bigrams, etc)
	// if a higher-order n-gram exists then its score overwrites the score of lower n-gram backoff
	int iComputed = 0;
	for(int i=iEmpty+1 ; i < m_iNGramOrder-1 ; ++i) {
	
		// get accumulated backoff score from higher order n-grams
		float fBackoffAcc = 0.0;
		for(int j=i ; j < m_iNGramOrder-1 ; ++j) {
			fBackoffAcc += lmArcBackoff[j]->fScore;
		}
		// compute lm-scores for all observed n-grams in this n-gram table
		LMState *state = m_states+lmArcBackoff[i]->iStateDest;
		LMArc *arcFinal = m_arcs+(state+1)->iArcBase;
		if (i != iEmpty+1) {		// backoff arc is the last arc, stop before it
			arcFinal--;
		}
		for(LMArc *arc = m_arcs+state->iArcBase; arc != arcFinal ; ++arc) {
			assert((arc->iLexUnit >= 0) && (arc->iLexUnit < iVocabularySize));
			fLMScores[arc->iLexUnit] = fBackoffAcc+arc->fScore;	
			++iComputed;
		}
		//printf("computed: %d\n",iComputed);
		iComputed = 0;
	}
	// (2) observed n-grams
	assert((m_arcs+((state+1)->iArcBase)-1)->iLexUnit == BACKOFF_ARC);
	for(LMArc *arc = m_arcs+state->iArcBase ; arc != m_arcs+((state+1)->iArcBase)-1 ; ++arc) {
		assert((arc->iLexUnit >= 0) && (arc->iLexUnit < iVocabularySize));	
		fLMScores[arc->iLexUnit] = arc->fScore;	
		++iComputed;
	}	
	//printf("computed: %d\n",iComputed);
	
	// set lm-score for filler units and sentence markers (0.0)
	// that kind of lexical units are after the lexical units for all unigrams
	//for(int i=m_iUnigrams ; i < iVocabularySize ; ++i) {
		//fLMScores[i] = 0.0;
		//printf("lex: %s\n",m_lexiconManager->getStrLexUnit(i));
	//}
	// there may or maynot be unigrams for the unknown symbol and the sentence markers 
	fLMScores[m_lexiconManager->m_lexUnitUnknown->iLexUnit] = 0.0;
	fLMScores[m_lexiconManager->m_lexUnitBegSentence->iLexUnit] = 0.0;
	fLMScores[m_lexiconManager->m_lexUnitEndSentence->iLexUnit] = 0.0;

	// check
	for(int i=0 ; i < iVocabularySize; ++i) {
		assert(fLMScores[i] != FLT_MAX);
	}
	
	delete [] lmArcBackoff;

	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	BVC_VERB << "seconds: " << (dTimeEnd-dTimeBegin)/1000.0 << " seconds" << endl;
}

// compute the likelihood of the given sequence of words
float LMFSM::computeLikelihood(const char *str) {

	VLexUnit vLexUnit;
	bool bAllKnown;
	m_lexiconManager->getLexUnits(str,vLexUnit,bAllKnown);
	if (vLexUnit.empty()) {
		return -1.0;
	}

	float fLikelihood = 0.0;
	float fLikelihoodAux;
	int iLMState = getInitialState();
	for(VLexUnit::iterator it = vLexUnit.begin() ; it != vLexUnit.end() ; ++it) {
		iLMState = updateLMState(iLMState,(*it)->iLexUnit,&fLikelihoodAux);
		fLikelihood += fLikelihoodAux;
	}
	fLikelihood += toFinalState(iLMState);

	return fLikelihood;
}

// print
void LMFSM::print() {

	cout << "-- language model FSM -------------------" << endl;
	cout << " n-gram order: " << LMManager::getStrNGram(m_iNGramOrder) << endl;
	cout << " # states: " << m_iStates << endl;
	cout << " # arcs:   " << m_iArcs << endl;
	cout << "-----------------------------------------" << endl;	
}

};	// end-of-namespace
