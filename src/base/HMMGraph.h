/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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

#ifndef HMMGRAPH_H
#define HMMGRAPH_H

#include "HMMManager.h"
#include "LexiconManager.h"

struct _FBNodeLexUnit;
struct _FBNodePhone;
struct _FBNodeHMM;

// lex-unit edge
typedef struct _FBEdgeLexUnit {
	_FBEdgeLexUnit *edgePrev;
	_FBEdgeLexUnit *edgeNext;
	_FBNodeLexUnit *nodePrev;
	_FBNodeLexUnit *nodeNext;
	LexUnit *lexUnit;	
} FBEdgeLexUnit;

// lex-unit node
typedef struct _FBNodeLexUnit {
	_FBEdgeLexUnit *edgePrev;
	_FBEdgeLexUnit *edgeNext;
	int iNode;
} FBNodeLexUnit;

typedef vector<FBNodeLexUnit*> VFBNodeLexUnit;
typedef list<FBNodeLexUnit*> LFBNodeLexUnit;
typedef map<FBNodeLexUnit*,bool> MFBNodeLexUnit;
typedef map<FBNodeLexUnit*,_FBNodePhone*> MFBNodeLexUnitNodePhone;

// phone edge
typedef struct _FBEdgePhone {
	_FBEdgePhone *edgePrev;
	_FBEdgePhone *edgeNext;
	_FBNodePhone *nodePrev;
	_FBNodePhone *nodeNext;
	unsigned char *iContextLeft;		// left context
	unsigned char *iContextRight;		// right context
	unsigned char iPhone;				// phone
	unsigned char iPosition;			// within-word position
	LexUnit *lexUnit;						// lexical unit
} FBEdgePhone;

// phone-node
typedef struct _FBNodePhone {
	_FBEdgePhone *edgePrev;				// previous node
	_FBEdgePhone *edgeNext;				// next node
	int iNode;
} FBNodePhone;

typedef vector<FBNodePhone*> VFBNodePhone;
typedef list<FBNodePhone*> LFBNodePhone;
typedef list<FBEdgePhone*> LFBEdgePhone;
typedef map<FBNodePhone*,bool> MFBNodePhone;

// HMM-edge (represents a HMM-state)
typedef struct _FBEdgeHMM {
	_FBEdgeHMM *edgePrev;
	_FBEdgeHMM *edgeNext;
	_FBNodeHMM *nodePrev;
	_FBNodeHMM *nodeNext;	
	HMMState *hmmStateEstimation;		// HMM-state		(physical n-phones -> local accumulators)
	HMMState *hmmStateUpdate;
	Accumulator *accumulator;			// accumulator		(logical n-phones -> global accumulators)
	int iActive;
	int iEdge;
	// these two fields are not needing to generate occupation stats but for recovering the best path
	LexUnit *lexUnit;						// to recover the best sequence of lexical units
	unsigned char iPosition;			// within-word position	
	// discriminative training
	//float fScoreAM;			// acoustic model score, needed to compute forward and backward score
	//float fScoreLM;			// language model score, needed to compute forward and backward score
	//float fScoreForward;		// forward score
	//float fScoreBackward;	// backward score
} FBEdgeHMM;

// HMM-node
typedef struct _FBNodeHMM {
	_FBEdgeHMM *edgePrev;
	_FBEdgeHMM *edgeNext;
	int iDistanceStart;
	int iDistanceEnd;
	int iNode;
} FBNodeHMM;

typedef vector<FBNodeHMM*> VFBNodeHMM;
typedef list<FBNodeHMM*> LFBNodeHMM;
typedef map<FBNodeHMM*,bool> MFBNodeHMM;
typedef map<FBNodePhone*,FBNodeHMM*> MFBNodePhoneNodeHMM;

#define NODE_STATE_NOT_VISITED				0
#define NODE_STATE_QUEUED						1
#define NODE_STATE_REMOVED						2
#define NODE_STATE_PROCESSED					3

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class HMMGraph {

	private:
	
		PhoneSet *m_phoneSet;							// phonetic symbol set
		HMMManager *m_hmmManagerEstimation;			// HMMs
		HMMManager *m_hmmManagerUpdate;				// HMMs
		LexiconManager *m_lexiconManager;			// lexicon manager	
		bool m_bHeadingTrailingSilence;				// whether silence will be inserted at the beginning and end of utterance
		bool m_bMultiplePronunciations;				// whether to use multiple pronunciations for each lexical unit
		VLexUnit m_vLexUnitOptional;					// optinal symbols to insert at word boundaries	
		bool m_bMergeAcrossDiffLexUnits;				// whether equivalent hmm-arcs belonging to different lex units are allowed to be merged
																// note: if we are using multiple pronunciations and want to recover the best one
																// 		then we must prevent equivalent arcs to be merged in order to easily recover it
																//       For statistics accumulation merging should be allowed because it will be faster
		unsigned char m_iContextSizeWW;				// within-word context size
		unsigned char m_iContextSizeCW;				// cross-word context size
		unsigned char m_iPhoneContextPadding;		// word-boundary phone
		
		// create a sausage of lexical units (it uses all the alternative pronunciations)
		FBNodeLexUnit *create(VLexUnitX &vLexUnitTranscription, VLexUnit &vLexUnitOptional, int *iNodes);
		
		// create a sausage of lexical units (it uses only the given pronunciation)
		FBNodeLexUnit *create(VLexUnit &vLexUnitTranscription, VLexUnit &vLexUnitOptional, int *iNodes);	
		
		// create a sausage of phones from the sausage of lexical units
		FBNodePhone *create(FBNodeLexUnit *nodeInitial, int iNodesLexUnit, int *iNodes, FBNodePhone **nodePhoneFinal);
		
		// create a HMM-graph from a phone-sausage 
		// - expansion might be needed because of the phonetic context
		FBNodeHMM *create(FBNodePhone *nodePhone, int *iNodesHMM);
		
		// destroy the lex-unit sausage
		static void destroy(FBNodeLexUnit *nodeInitial);
		
		// destroy the lex-unit graph
		static void destroy(FBNodeLexUnit *nodeInitial, int iNodes);
			
		// destroy the phone-graph
		static void destroy(FBNodePhone *nodeInitial, int iNodes);
		
		// propagate phonetic context across all phones
		void propagateContext(FBNodePhone *nodeInitial, int iNodesOriginal, int *iNodes);
		
		// print the sausage of lexical units
		void print(FBNodeLexUnit *nodeInitial);
		
		// print the graph of phones 		
		void print(FBNodePhone *nodeInitial);
		
		// enumerate the nodes
		int renumberNodes(FBNodePhone *nodeInitial);	
		
		
		inline FBNodeLexUnit *newFBNodeLexUnit(int iNode) {
		
			FBNodeLexUnit *node = new FBNodeLexUnit;
			node->edgePrev = NULL;
			node->edgeNext = NULL;
			node->iNode = iNode;
		
			return node;
		}
		
		inline FBEdgeLexUnit *newFBEdgeLexUnit(FBNodeLexUnit *nodePrev, FBNodeLexUnit *nodeNext, LexUnit *lexUnit) {
		
			FBEdgeLexUnit *edge = new FBEdgeLexUnit;
		
			edge->nodePrev = nodePrev;
			edge->nodeNext = nodeNext;
			edge->edgePrev = nodePrev->edgeNext;
			nodePrev->edgeNext = edge;
			edge->edgeNext = nodeNext->edgePrev;
			nodeNext->edgePrev = edge;
			edge->lexUnit = lexUnit;
		
			return edge;
		}
		
		inline FBNodePhone *newFBNodePhone(int iNode) {
		
			FBNodePhone *node = new FBNodePhone;
			node->edgePrev = NULL;
			node->edgeNext = NULL;
			node->iNode = iNode;
		
			return node;
		}	
		
		inline FBEdgePhone *newFBEdgePhone(FBNodePhone *nodePrev, FBNodePhone *nodeNext) {
		
			FBEdgePhone *edge = new FBEdgePhone;
			
			edge->nodePrev = nodePrev;
			edge->nodeNext = nodeNext;
			edge->edgePrev = nodePrev->edgeNext;
			nodePrev->edgeNext = edge;
			edge->edgeNext = nodeNext->edgePrev;
			nodeNext->edgePrev = edge;
			
			edge->iContextLeft = NULL;
			edge->iContextRight = NULL;
			
			return edge;
		}
		
		inline FBNodeHMM *newFBNodeHMM(int iNode) {
		
			FBNodeHMM *node = new FBNodeHMM;
			node->edgePrev = NULL;
			node->edgeNext = NULL;
			node->iNode = iNode;
		
			return node;
		}
		
		inline FBEdgeHMM *newFBEdgeHMM(FBNodeHMM *nodePrev, FBNodeHMM *nodeNext, HMMState *hmmStateEstimation, 
			HMMState *hmmStateUpdate, Accumulator *accumulator, LexUnit *lexUnit, unsigned char iPosition) {
		
			FBEdgeHMM *edge = new FBEdgeHMM;
			
			edge->nodePrev = nodePrev;
			edge->nodeNext = nodeNext;
			edge->edgePrev = nodePrev->edgeNext;
			nodePrev->edgeNext = edge;
			edge->edgeNext = nodeNext->edgePrev;
			nodeNext->edgePrev = edge;
			
			edge->hmmStateEstimation = hmmStateEstimation;
			edge->hmmStateUpdate = hmmStateUpdate;
			edge->accumulator = accumulator;
			edge->lexUnit = lexUnit;
			edge->iPosition = iPosition;
			
			return edge;
		}		
		
		// check that for each node all the successors and predecessors have attached the same phone
		void check(FBNodePhone *nodeInitial);
		
		// compact the HMM-graph by applying forward/backward node merging
		FBNodeHMM **compact(FBNodeHMM *nodeInitial, int iNodes, int *iNodesAfterCompacting, 
			int *iEdgesAfterCompacting, FBNodeHMM **nodeHMMFinal);	
			
		// return whether the nodes have the same set of predecessors
		bool samePredecessors(FBNodeHMM *node1, FBNodeHMM *node2);	
		
		// return whether the nodes have the same set of successors
		bool sameSuccessors(FBNodeHMM *node1, FBNodeHMM *node2);
			
		// return whether two HMM-edges are equal for merging purposes
		inline bool equal(FBEdgeHMM *edge1, FBEdgeHMM *edge2) {
		
			// do we allow to merge edges beloning to different lexical units?
			if (m_bMergeAcrossDiffLexUnits == false) {
				if (edge1->lexUnit->iLexUnitPron != edge2->lexUnit->iLexUnitPron) {
					return false;
				}
			}
		
			// logical accumulators
			if (m_hmmManagerUpdate->areAccumulatorsLogical()) {
				assert(edge1->hmmStateUpdate == NULL);
				assert(edge2->hmmStateUpdate == NULL);
				return ((edge1->hmmStateEstimation == edge2->hmmStateEstimation) && 
					(edge1->accumulator == edge2->accumulator));	
			} 
			// physical accumulators
			else {
				assert(edge1->accumulator == NULL);
				assert(edge2->accumulator == NULL);
				return ((edge1->hmmStateEstimation == edge2->hmmStateEstimation) && 
					(edge1->hmmStateUpdate == edge2->hmmStateUpdate));	
			}
		}

		// compute the distance from each edge to the initial and final edges in the graph
		// - it applies Dijkstra both ways
		void computeDistances(FBNodeHMM **nodes, FBNodeHMM *nodeInitial, 
			FBNodeHMM *nodeFinal, int iNodes, int iEdges);
			
		// allocates memory for a blank context
		inline unsigned char *newBlankContext() {
			
			int iSize = max(m_iContextSizeCW,m_iContextSizeWW);
			unsigned char *iContext = new unsigned char[iSize];
			for(int i=0 ; i < iSize ; ++i) {
				iContext[i] = m_iPhoneContextPadding;
			}
		
			return iContext;
		}
		
		// propagate left-context from one edge to the next
		// note: destination context must be already filled with padding symbols
		inline void shiftLeftContext(FBEdgePhone *edge, unsigned char *iContextLeft) {
		
			// cross-word context propagation 
			if ((edge->iPosition == WITHIN_WORD_POSITION_END) || 
				(edge->iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
					
				iContextLeft[m_iContextSizeCW-1] = edge->iPhone;
				for(int i=1 ; i < min((int)m_iContextSizeCW,(int)(edge->lexUnit->vPhones.size())) ; ++i) {
					iContextLeft[m_iContextSizeCW-i-1] = edge->iContextLeft[m_iContextSizeCW-i];
				}
			}
			// within-word context propagation
			else {
				iContextLeft[m_iContextSizeWW-1] = edge->iPhone;
				// word-initial phone to second phone: do nothing else
				if (edge->iPosition == WITHIN_WORD_POSITION_START) {
					// do nothing
				} 
				// word-internal to word-internal or word-end
				else {
					for(int i=0 ; i < m_iContextSizeWW-1 ; ++i) {
						iContextLeft[i] = edge->iContextLeft[i+1]; 
					}
				}
			}
		
			return;
		}

		// propagate right-context from one edge to the previous edge
		// note: destination context must be already filled with padding symbols
		inline void shiftRightContext(FBEdgePhone *edge, unsigned char *iContextRight) {
		
			// cross-word context propagation 
			if ((edge->iPosition == WITHIN_WORD_POSITION_START) || 
				(edge->iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
					
				iContextRight[0] = edge->iPhone;
				for(int i=1 ; i < min((int)m_iContextSizeCW,(int)(edge->lexUnit->vPhones.size())) ; ++i) {
					iContextRight[i] = edge->iContextRight[i-1];
				}
			}
			// within-word context propagation
			else {
				iContextRight[0] = edge->iPhone;
				// word-final phone to phone before last phone: do nothing else
				if (edge->iPosition == WITHIN_WORD_POSITION_END) {
				} 
				// word-internal to word-internal or word-end
				else {
					for(int i=1 ; i < m_iContextSizeWW ; ++i) {
						iContextRight[i] = edge->iContextRight[i-1]; 
					}
				}
			}
		
			return;
		}
	
	public:

		// constructor
		HMMGraph(PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMManager *hmmManagerEstimation, 
			HMMManager *hmmManagerUpdate, bool bMultiplePronunciations, VLexUnit &vLexUnitOptional, bool bMergeAcrossDiffLexUnits = true);

		// destructor
		~HMMGraph();
		
		// create a graph from a sequence of lexical units
		// - handles multiple pronunciations
		// - handles optional symbols (typically silence+fillers)	
		FBNodeHMM **create(VLexUnit &vLexUnitTranscription, int *iNodes, int *iEdges, 
			FBNodeHMM **nodeHMMInitial, FBNodeHMM **nodeHMMFinal);
			
		// print the graph of HMMs
		void print(FBNodeHMM *nodeInitial);	
		
		// destroy the HMM-graph
		static void destroy(FBNodeHMM *nodeInitial, int iNodes);
		
		// destroy the HMM-graph
		static void destroy(FBNodeHMM **nodes, int iNodes);

};

#endif
