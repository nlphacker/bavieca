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


#ifndef WFSABUILDER_H
#define WFSABUILDER_H

using namespace std;

#include <deque>
#include <list>
#include <string>
#include <vector>
#include <iomanip>

#if defined __linux__ || defined __APPLE__
using namespace __gnu_cxx;
#include <ext/hash_map>
#elif _WIN32
#include <hash_map>
#else
	#error "unsupported platform"
#endif

#include "Global.h"
#include "LexiconManager.h"
#include "WFSAcceptor.h"

namespace Bavieca {

class HMMManager;
class LMManager;
class PhoneSet;

#define WEIGHT_PUSHING_SEMIRING_NONE				0			// the weights are not changed (no weight pushing)
#define WEIGHT_PUSHING_SEMIRING_TROPICAL			1			// tropical semiring
#define WEIGHT_PUSHING_SEMIRING_LOG					2			// log-semiring

typedef struct _Transition Transition;
typedef struct _Transition2 Transition2;
typedef vector<_Transition*> VTransition;

typedef struct _State {
	unsigned int iId;						// state id
	VTransition vTransition;			// transitions 
	_State *stateParent;					// parent state (NULL for root states)	
	Transition *transitionBW;			// transition from the parent state that goes to the state (NULL for root states)
} State;

typedef struct {
	State **statePredecessor;
	Transition **transitionPredecessor;
	int iPredecessors;
} StatePredecessors;

typedef struct {
	unsigned int iId;						// state id
	unsigned int iTransitions;			// # transitions
	Transition2 *transitions;			// transitions 
} State2;

typedef vector<State*> VState;
typedef list<State*> LState;

typedef struct {
	State *stateLeaf;								// leaf state (the lexical unit transition would go after this state)
	unsigned int iRightContextGroup;			// index of a right context group
} LeafRightContext;

typedef vector<LeafRightContext> VLeafRightContext;

typedef struct {
	unsigned int *iLeafRightContext;
	LeafRightContext **leafRightContext;	// one for each lexical unit (it provides an array of pairs (leaf,group of right contexts that it models)	
	State *stateRoot;
} LexiconTree;

// maps a leaf transition in the final graph (before connections are made) to the right context the leaf was built for
typedef map<Transition*,unsigned char > MLeafTransitionRightContext;
typedef map<Transition*,vector<unsigned char> > MLeafTransitionVRightContext;


//#define MAX_RIGHT_CONTEXT_GROUPS			UINT_MAX
#define MAX_RIGHT_CONTEXT_GROUPS			2000

// maps a pair (stateRootL,right context) -> state in final graph (starting state)

//typedef vector<pair<Transition*,unsigned char> > VLeafTransitionRightContext;

#define WFST_EPSILON		UINT_MAX

//#define TRANSITION_TYPE_LEXICAL_UNIT			0
//#define TRANSITION_TYPE_HMM						1
//#define TRANSITION_TYPE_EPSILON				2

typedef struct _Transition {
	unsigned int iSymbol;					// input/output symbol
	float fWeight;								// weight 
	State *state;								// destination state
	unsigned int iRightContextGroup;		// right context group
} Transition;

typedef struct _Transition2 {
	unsigned int iSymbol;					// input/output symbol
	float fWeight;								// weight 
	State *state;								// destination state
} Transition2;

typedef struct {
	LexUnit *lexUnit;					// lexical unit
	unsigned int iHMMStates;		// number of HMM-states (to index the buffer)
	State *stateLeaf;
	unsigned int iRightContextGroup;
	State *stateGTo;					// state in G
	float fWeight;
} LeafToProcess;

typedef list<LeafToProcess*> LLeafToProcess;

// ad-hoc functions to use State* as the key of a hash_map structure
struct stateHashMapFunctions
{

#if defined __linux__ || defined __APPLE__
	
	// comparison function (used for matching, comparison for equality)
	bool operator()(const State *state1, const State *state2) const {
		
		return (state1 == state2);
	}

#elif _WIN32

	static const size_t bucket_size = 4;
	static const size_t min_buckets = 8;
	
	// comparison function (used for ordering)
	bool operator()(const State *state1, const State *state2) const {
		
		return (state1 < state2);
	}

#endif
	
	// hash function (elements in the subset are not sorted)
	size_t operator()(const State *state) const {	
			
		return (size_t)state;
	}
};

#if defined __linux__ || defined __APPLE__
typedef hash_map<State*,bool,stateHashMapFunctions,stateHashMapFunctions> MStateBool;
typedef hash_map<State*,unsigned int,stateHashMapFunctions,stateHashMapFunctions> MStateSymbol;
#elif _WIN32
typedef hash_map<State*,bool,stateHashMapFunctions> MStateBool;
typedef hash_map<State*,unsigned int,stateHashMapFunctions> MStateSymbol;
#endif

typedef map<State*,float> MStateWeight;

// ad-hoc functions to use State* as the key of a hash_map structure
struct statePairHashMapFunctions
{
	// comparison function
	bool operator()(const pair<State*, unsigned char> pair1, const pair<State*, unsigned char> pair2) const {		
		
		return ((pair1.first == pair2.first) && (pair1.second == pair2.second));
	}
	
	// hash function (elements in the subset are not sorted)
	size_t operator()(const pair<State*, unsigned char> pair) const {	
			
		return (size_t)(((unsigned long)pair.first)^pair.second);
	}
};

// maps a pair (State in G, left context) to the root state in L
typedef hash_map<pair<State*, unsigned char>, State*,statePairHashMapFunctions,statePairHashMapFunctions> MGStateLRoot;

typedef struct {
	State *state;
	StateX *stateX;
} StatePair;

typedef vector<StatePair> VStatePair;


// ad-hoc functions to use State* as the key of a hash_map structure
struct stateEquivalenceHashFunctions
{
	// comparison function: two states are equal if their transitions are identical (symbol, weight, destination state)
	bool operator()(const State *state1, const State *state2) const {
	
		if (state1->vTransition.size() != state2->vTransition.size()) {
			return false;
		}
	
		bool bFound = false;
		for(VTransition::const_iterator it = state1->vTransition.begin() ; it != state1->vTransition.end() ; ++it) {
			bFound = false;
			for(VTransition::const_iterator jt = state2->vTransition.begin() ; jt != state2->vTransition.end() ; ++jt) {
				if (((*it)->iSymbol == (*jt)->iSymbol) && 
					((*it)->fWeight == (*jt)->fWeight) && 
					((*it)->iRightContextGroup == (*jt)->iRightContextGroup) && 
					((*it)->state == (*jt)->state)) {
					bFound = true;
					break;
				}
			}
			if (bFound == false) {
				break;
			}
		}
		
		return bFound;
	}
	
	// hash function (elements in the subset are not sorted)
	size_t operator()(const State *state) const {	
	
		size_t key = 0;
		for(VTransition::const_iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
			key ^= ((size_t)(*it)->iSymbol)^((size_t)(*it)->state);	
		}
			
		return key;
	}
};


// map used to perform state equivalence
typedef hash_map<State*,State*,stateEquivalenceHashFunctions,stateEquivalenceHashFunctions> MStateState;

// Dijkstra
#define STATE_STATUS_UNSEEN				0
#define STATE_STATUS_QUEUED				1	
#define STATE_STATUS_PROCESSED			2

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class WFSABuilder {

	private:
	
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManager;
		LexiconManager *m_lexiconManager;
		LMManager *m_lmManager;
		unsigned char m_iNGram;
		float m_fLMScalingFactor;
		
		// stats
		int m_iLexUnitsTotal;
		int m_iHMMStatesTotal;
		
		// build the grammar (from an n-gram) as an acceptor, it returns the initial state
		bool buildG(State **states, unsigned int *iStatesG, unsigned int *iTransitionsG, State **stateInitialG, State **stateFinalG);
		
		static inline State *newState(unsigned int iId, Transition *transitionBW, State *stateParent) {
		
			State *state = new State;
			state->iId = iId;
			state->transitionBW = transitionBW;
			state->stateParent = stateParent;
			
			return state;
		}
		
		static inline Transition *newTransition(unsigned int iSymbol, float fWeight, State *state) {
		
			Transition *transition = new Transition;
			transition->iSymbol = iSymbol;
			transition->fWeight = fWeight;
			transition->state = state;	
			transition->iRightContextGroup = UINT_MAX;	
			
			//assert(state != NULL);
			/*if (state->iId == 81) {
				bool bStop = 0;
			}*/
			
			return transition;
		}
		
		// build a word-internal context dependent decoding network or a context-independent one
		WFSAcceptor *buildWordInternal2();
		
		// build a word-internal context dependent decoding network or a context-independent one
		WFSAcceptor *buildWordInternal();
		
		// build a cross-word context dependent decoding network 
		WFSAcceptor *buildCrossWord2();
		
		// build a cross-word context dependent decoding network 
		WFSAcceptor *buildCrossWord();
		
		// create the word-internal/context-independent lexical tree
		State *buildWordInternalLexiconTree(map<LexUnit*,State*> &mLexUnitLeafState, unsigned int *iTreeDepth);	
		
		// create the left-context-dependent map of trees
		LexiconTree **buildCrossWordLexiconTrees(unsigned int *iTreeDepth, vector<unsigned char*> &vRightContextGroups);
		
		// post-order state numbering
		void preOrderNumbering(State *state, unsigned int *iId);	
		
		// return the maximum tree depth
		int getTreeDepth(State *state, int iDepth = 1);
		
		// return the number of leaves
		unsigned int getLeaves(State *state) {	
			
			if (state->vTransition.empty() == false) {
				unsigned int iAcc = 0;
				for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
					iAcc += getLeaves((*it)->state);
				}
				return iAcc;
			} 
			else {
				return 1;
			}	
		}
		
		// compare states by their number
		static inline bool compareStateNumber(const LeafToProcess *leaf1, const LeafToProcess *leaf2) {
		
			return (leaf1->stateLeaf->iId < leaf2->stateLeaf->iId);
		}
		
		// make each state in the acceptor (WFSA) accessible only by transitions with the same input symbol (for more efficient decoding)
		void equalizeInput(State *stateInitial, unsigned int iStates, unsigned int iTransitions, unsigned int *iStatesFinal, unsigned int *iTransitionsFinal);
		
		// optimize the acceptor for decoding
		WFSAcceptor *optimize(State *stateInitial, unsigned int iStates, unsigned int iTransitions, State *stateFinalG);	
		
		// detroy a lexicon tree
		void destroyTree(State *state);	
		
		// check auxiliar root node: it cannot have multiple leaves with HMMs for the same phone
		/*void checkAuxiliarRootL(State *state) {
		
			bool bPhones[m_phoneSet->size()];
			for(unsigned char i = 0 ; i < m_phoneSet->size() ; ++i) {
				bPhones[i] = false;
			}
			for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
				if ((*it)->iSymbol < LEX_UNIT_TRANSITION) {
					unsigned char iPhone = m_hmmManager->getPhoneticSymbolFromHMMStateDecoding((*it)->iSymbol);
					if (bPhones[iPhone] == false) {
						bPhones[iPhone] = true;
					} else {
						assert(0);
					}
				}
			}
		}*/	
		
		// check the correctness of a transition symbol (debugging)
		void checkTransitionSymbolCorrectness(unsigned int iSymbol) {
		
			// make sure the transition symbol is correct
			// epsilon:
			if (iSymbol == EPSILON_TRANSITION) {
			
			}
			// lexical unit
			else if (iSymbol & LEX_UNIT_TRANSITION) {
				int iLexUnit = iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				assert((iLexUnit >= 0) && (iLexUnit < m_iLexUnitsTotal));
			}
			// HMM-state
			else {
				assert(iSymbol < LEX_UNIT_TRANSITION);
				int iHMMState = iSymbol;
				assert((iHMMState >= 0) && (iHMMState < m_iHMMStatesTotal));	
			}		
		}
		
		// print the array of states
		void printStates(State *states, unsigned int iStates, State *stateInitial) {
		
			cout << "--------------------------------------------------------" << endl;
			for(unsigned int i=0 ; i<iStates ; ++i) {
				if (states+i == stateInitial) {
					cout << "initial state: " << stateInitial << " (" << stateInitial->iId << ")" << endl;
				}
				printStateTransitions(states+i);
			}	
			cout << "--------------------------------------------------------" << endl;
		}
		
		void printStateTransitions(State *state) {
		
			cout << "state id: " << state->iId << endl;
			for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
				printTransition(*it);	
			}
		}
		
		void printTransition(Transition *transition) {
		
			ostringstream oss;
			// epsilon:
			if (transition->iSymbol == EPSILON_TRANSITION) {
				oss << "eps                 ";
			}
			// lexical unit
			else if (transition->iSymbol & LEX_UNIT_TRANSITION) {
				int iLexUnit = transition->iSymbol & LEX_UNIT_TRANSITION_COMPLEMENT;
				oss << setw(20) << m_lexiconManager->getStrLexUnit(iLexUnit);
			}
			// HMM-state
			else {
				assert(transition->iSymbol < LEX_UNIT_TRANSITION);
				oss << "state:  " << setw(12) << transition->iSymbol;
			}
			oss << " " << std::setw(12) << std::setiosflags(ios::fixed) << std::setprecision(4) 
				<< transition->fWeight << " " << transition->state 
				<< " (" << transition->state->iId << ")" << endl;
			cout << oss.str();
		}
		
		// return a list containing the left context of filler lexical units (including silence)
		bool *getFillerFirstPhones(unsigned char *iElements);		
		
		// return a data structure containing the left context of G-nodes
		unsigned char **getLeftContextGNodes(State *states, unsigned int iStatesG, int iStateGFinal);
	
		// apply weight pushing to a state
		/*inline float applyWeightPushing(State *state) {
		
			float fWeightMax = -FLT_MAX;
			for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
				if ((*it)->fWeight > fWeightMax) {
					fWeightMax = (*it)->fWeight;
				}	
			}
			for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
				(*it)->fWeight -= fWeightMax;
				//printf("%x %f\n",*it,(*it)->fWeight);
			}
			
			return fWeightMax;
		}*/
		
		// apply weight pushing to a state (log semiring)
		//inline float applyWeightPushing(State *state, unsigned char iSemiring = WEIGHT_PUSHING_SEMIRING_LOG) {
		inline float applyWeightPushing(State *state, unsigned char iSemiring = WEIGHT_PUSHING_SEMIRING_TROPICAL) {
		
			// no weight pushing
			if (iSemiring == WEIGHT_PUSHING_SEMIRING_NONE) {
			
				return 0.0;
		
			// tropical semiring
			} else if (iSemiring == WEIGHT_PUSHING_SEMIRING_TROPICAL) {
		
				float fWeightMax = -FLT_MAX;
				for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
					if ((*it)->fWeight > fWeightMax) {
						fWeightMax = (*it)->fWeight;
					}	
				}
				for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
					(*it)->fWeight -= fWeightMax;
				}
				
				return fWeightMax;
			}
			// log semiring
			else { 
				assert(iSemiring == WEIGHT_PUSHING_SEMIRING_LOG);
		
				double dLogSum = 0.0;
				for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
					dLogSum += exp((*it)->fWeight);
				}
				dLogSum = log(dLogSum);
				assert(finite(dLogSum));
				for(VTransition::iterator it = state->vTransition.begin() ; it != state->vTransition.end() ; ++it) {
					(*it)->fWeight -= dLogSum;
				}
				
				return dLogSum;	
			} 
		}
		
		// sort the states in G by reverse topological order of epsilon transitions: the idea is that states with epsilon transitions
		// will be processed first so when the epsilon transitions are seen the destination state is already processed
		void sortReverseTopologicalOrderEpsilonTransitions(State *stateRootG);	
		
		// perform global weight pushing, which is the first step of weighted minimization, the second step is the actual minimization
		// note: it only works for the tropical semiring, for the log-semirint there are alternative methods
		// note: there is always only one final state which is the destination state for end-of-sentence transitions
		void globalWeightPushing(State *stateInitial, State *stateFinal, int iStates, 
			float *fWeightInitial, float *fWeightFinal);
		
	
	public:
    
		// contructor
		WFSABuilder(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager, LMManager *lmManager, unsigned char iNGram, float fLMScalingFactor);

		// destructor
		~WFSABuilder();
		
		// build the decoding network as a WFSAcceptor
		WFSAcceptor *build();	

};

};	// end-of-namespace

#endif
