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


#ifndef HYPOTHESISLATTICE_H
#define HYPOTHESISLATTICE_H

using namespace std;

#include <vector>
#include <list>
#include <map>
#include <iomanip>

#include <float.h>

#include "Global.h"
#include "HMMManager.h"
#include "LexiconManager.h"

namespace Bavieca {

class Alignment;
class BestPath;
class ExceptionBase;
class LMManager;
class Mappings;
class PhoneSet;

// lattice properties
// mandatory properties
#define LATTICE_PROPERTY_VERSION					"[version]"
#define LATTICE_PROPERTY_NODES					"[nodes]"
#define LATTICE_PROPERTY_EDGES					"[edges]"
#define LATTICE_PROPERTY_FRAMES					"[frames]"
// optional properties
#define LATTICE_PROPERTY_AM_PROB					"[am-prob]"
#define LATTICE_PROPERTY_HMMS						"[hmms]"
#define LATTICE_PROPERTY_PHONE_ALIGN			"[ph-align]"
#define LATTICE_PROPERTY_PHONE_ACCURACY		"[phone-acc]"			// discriminative training
#define LATTICE_PROPERTY_INSERTION_PENALTY	"[ip]"			
#define LATTICE_PROPERTY_LM_PROB					"[lm-prob]"
#define LATTICE_PROPERTY_NGRAM					"[n-gram]"
#define LATTICE_PROPERTY_BEST_PATH				"[best-path]"
#define LATTICE_PROPERTY_FWD_PROB				"[fwd-prob]"
#define LATTICE_PROPERTY_BWD_PROB				"[bwd-prob]"
#define LATTICE_PROPERTY_PP						"[pp]"
#define LATTICE_PROPERTY_CONFIDENCE				"[conf]"

// to keep lattice properties
typedef map<string,string> MProperty;

#define LATTICE_NODE_STATE_UNSEEN			0
#define LATTICE_NODE_STATE_QUEUED			1
#define LATTICE_NODE_STATE_PROCESSED		2
#define LATTICE_NODE_STATE_DELETED			3

#define LATTICE_FORMAT_VERSION			"0.1"

// phone-level alignments
typedef struct {
	unsigned char iPhone;						// phonetic symbol relative to the phonetic symbol set
	unsigned char iPosition;					// position of the phone within the word
	int iStateBegin[NUMBER_HMM_STATES];		// first time frame aligned to each of the HMM-states in the phone
	int iStateEnd[NUMBER_HMM_STATES];		// first time frame aligned to each of the HMM-states in the phone
	int iHMMState[NUMBER_HMM_STATES];		// index of each of the HMM-states in the phone (clustered HMM-state)
} LPhoneAlignment;

typedef vector<LPhoneAlignment*> VLPhoneAlignment;

// the lattice structure is defined so it can be traversed forward and backward

struct _LNode;

// edge in the lattice
typedef struct _LEdge {
   int iEdge;									// edge-id
   int iFrameStart;							// start frame
   int iFrameEnd;								// end frame
	LexUnit *lexUnit;							// lexical unit (including pronunciation)
	float fScoreAM;							// acoustic score
	float fScoreLM;							// language model score
	float fInsertionPenalty;				// insertion penalty
	_LNode *nodePrev;							// incoming node
	_LNode *nodeNext;							// outgoing node
	_LEdge *edgePrev;							// list of edges going to the node
	_LEdge *edgeNext;							// list of edges coming from the node
	// phone-level alignment
	int iPhones;								// number of phones
	LPhoneAlignment *phoneAlignment;		// phone-level alignment
	// forward/backward scores
	double dScoreForward;					// forward score (different meaning for n-best and post prob)
	double dScoreBackward;					// backward score (different meaning for n-best and post prob)
	// posterior probabilities + confidence estimation
	float fPP;									// posterior probability
	float fConfidence;						// confidence score	
	// phonetic context
	unsigned char *iContextLeft;			// left phonetic context
	unsigned char *iContextRight;			// right phonetic context
	// phone accuracy (used for discriminative training)
	float *fPhoneAccuracy;
	// auxiliar fields
	bool bTouched;
	bool bBestPath;
	_LEdge *edgeNextAux;	
	int iLMState;								// language model state
	int iLMStatePrev;							// language model state of previous edge
} LEdge;

typedef vector<LEdge*> VLEdge;
typedef list<LEdge*> LLEdge;
typedef map<LEdge*,bool> MLEdge;
typedef map<int,LEdge*> MIntLEdge;

// node in the lattice
typedef struct _LNode {
	int iNode;
	int iFrame;						// ending frame of edges arriving at this node
	LEdge *edgePrev;	
	LEdge *edgeNext;
	// auxiliar fields
	bool bTouched;
} LNode;

typedef vector<LNode*> VLNode;
typedef list<LNode*> LLNode;
typedef map<LNode*,bool> MLNode;

// to keep lattice depth information
typedef struct {
	int iFramesEdges;
	int iFrames;
	float fDepth;
} LatticeDepth;

// to keep the Lattice WER information
typedef struct {
	int iWordsReference;				// # words in the reference
	int iInsertions;					// # insertions
	int iDeletions;					// # deletions
	int iSubstitutions;				// # substitutions
	int iCorrect;						// # correct words
	int iErrors;						// # errors
	int iOOV;							// # Out Of Vocabulary words
	float fWER;							// word error rate
	VLexUnit vLexUnitBest;			// best hypothesis in the lattice
} LatticeWER;

// alignment events
#define ALIGNMENT_EVENT_CORRECT				0
#define ALIGNMENT_EVENT_INSERTION			1
#define ALIGNMENT_EVENT_DELETION 			2
#define ALIGNMENT_EVENT_SUBSTITUTION		3
#define ALIGNMENT_EVENT_SKIP					4

#define ALIGNMENT_EVENTS_TOTAL				4

// alignment penalties
#define ALIGNMENT_PENALTY_INSERTION			1.0
#define ALIGNMENT_PENALTY_DELETION			1.0
#define ALIGNMENT_PENALTY_SUBSTITUTION		1.0

// beam size for Word Error Rate computation
#define LATTICE_WER_COMPUTATION_BEAM		5.0

// structure for the backtrace (Lattice WER computation)
typedef struct _LWERHistoryItem {
	int iPrev;									// previous history item
	LEdge *edge;								// edge in the lattice
	int iReference;							// reference index
	unsigned char iAlignmentEvent;		// ins|del|subs
	bool bUsed;
} LWERHistoryItem;

// token to compute the lattice WER
typedef struct _LWERToken {
	float fScore;					// alignment score
	LNode *node;					// node in the lattice
	int iReference;				// reference index
	int iHistoryItem;				// last history-item (back-trace)
	_LWERToken *prev;				// previous token in the list of active tokens for a given lattice node
} LWERToken;

// confidence measures (all of them are based on posterior probabilities)
#define CONFIDENCE_MEASURE_NONE					0
#define CONFIDENCE_MEASURE_POSTERIORS			1
#define CONFIDENCE_MEASURE_ACCUMULATED			2
#define CONFIDENCE_MEASURE_MAXIMUM				3

// confidence measures (string format)
#define CONFIDENCE_MEASURE_NONE_STR						"none"
#define CONFIDENCE_MEASURE_POSTERIORS_STR				"posteriors"
#define CONFIDENCE_MEASURE_ACCUMULATED_STR			"accumulated"
#define CONFIDENCE_MEASURE_MAXIMUM_STR					"maximum"

// lattice rescoring methods
#define RESCORING_METHOD_LIKELIHOOD						"likelihood"
#define RESCORING_METHOD_POSTERIORS						"pp"

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class HypothesisLattice {

	private:
	
		PhoneSet *m_phoneSet;						// phonetic symbol set
		LexiconManager *m_lexiconManager;		// lexicon manager
		bool m_bVerbose;								// verbose output
		int m_iNodes;									// # nodes in the lattice
		int m_iEdges;									// # edges in the lattice
		LNode *m_lnodeInitial;						// initial node in the lattice
		LNode *m_lnodeFinal;							// final node in the lattice
		LNode **m_lnodes;								// nodes in the lattice
		LEdge **m_ledges;								// edges in the lattice
		
		// lattice properties
		MProperty m_mProperties;
		int m_iFrames;									// # feature frames in the utterance
		
		// scaling factors
		float m_fAMScalingFactor;
		float m_fLMScalingFactor;
		
		// Lattice WER computation
		int m_iMaxActiveTokens;						// maximum number of active tokens
		LWERToken *m_tokenCurrent;					// current active tokens
		int m_iActiveTokensCurrent;
		LWERToken *m_tokenNext;						// next active tokens
		int m_iActiveTokensNext;
		LWERToken **m_tokenNode;					// table of nodes used to keep the active tokens at each node
		LWERHistoryItem *m_historyItems;			// history items allocated
		int m_iHistoryItems;							// # history items allocated
		int m_iHistoryItemAvailable;				// next history item that is available
		LWERToken m_tokenBest;
		
		// hmm-marking
		unsigned char m_iContextSizeWW;				// within-word context size
		unsigned char m_iContextSizeCW;				// cross-word context size
		unsigned char m_iPhoneContextPadding;		// word-boundary phone

	public:

		// constructor
		HypothesisLattice(PhoneSet *phoneSet, LexiconManager *lexiconManager, bool bVerbose = false);

		// destructor
		~HypothesisLattice();
		
		// build the lattice container
		void buildContainer(LNode *lnodeInitial, LNode *lnodeFinal);
		
		// build the lattice container (reusing both the initial and final nodes)
		void buildContainer();
		
		// destroy the lattice
		void destroy();
		
		// set a property
		void setProperty(const char *strProperty, const char *strValue) {
		
			m_mProperties[strProperty] = strValue;
		}
		
		// set a property
		void setProperty(const char *strProperty, int iValue) {
		
			stringstream strValue;
			strValue << iValue;
			setProperty(strProperty,strValue.str().c_str());
		}
		
		// remove a property
		void removeProperty(const char *strProperty) {
		
			MProperty::iterator it = m_mProperties.find(strProperty);
			if (it != m_mProperties.end()) {
				m_mProperties.erase(it);
			}
		}	
		
		// check that a property is defined
		bool isProperty(const char *strProperty) {
		
			return (m_mProperties.find(strProperty) != m_mProperties.end());
		}
		
		// check that a property has the given value
		bool checkProperty(const char *strProperty, const char *strValue) {
		
			if (!isProperty(strProperty)) {
				return false;
			}
			
			return (m_mProperties[strProperty].compare(strValue) == 0);
		}
		
		// get a property value
		const char *getPropertyValue(const char *strProperty) {
		
			if (!isProperty(strProperty)) {
				return NULL;
			}
		
			return m_mProperties[strProperty].c_str();	
		}
		
		// creates a new lattice node 
		static LNode *newNode(int iFrame) {
		
			LNode *lNode = new LNode;
			lNode->iFrame = iFrame;
			lNode->edgePrev = NULL;
			lNode->edgeNext = NULL;
			lNode->bTouched = false;
		
			return lNode;
		}
		
		// creates a new lattice edge
		static LEdge *newEdge(int iFrameStart, int iFrameEnd, LexUnit *lexUnit, float fScoreAM, float fScoreLM, float fInsertionPenalty) {
		
			assert(iFrameEnd > iFrameStart);
		
			LEdge *ledge = new LEdge;
			ledge->iFrameStart = iFrameStart;
			ledge->iFrameEnd = iFrameEnd;
			ledge->lexUnit = lexUnit;
			ledge->fScoreAM = fScoreAM;
			ledge->fScoreLM = fScoreLM;
			ledge->fInsertionPenalty = fInsertionPenalty;
			
			ledge->nodePrev = NULL;
			ledge->nodeNext = NULL;
			ledge->edgePrev = NULL;
			ledge->edgeNext = NULL;	
			
			ledge->iPhones = 0;
			ledge->phoneAlignment = NULL;
			ledge->fPhoneAccuracy = NULL;
			
			ledge->iContextLeft = NULL;
			ledge->iContextRight = NULL;
			
			ledge->bTouched = false;
			ledge->bBestPath = false;
		
			return ledge;
		}
		
		// creates a copy of a lattice edge
		static LEdge *newEdge(LEdge *edgeOriginal) {
		
			LEdge *ledge = new LEdge;
			ledge->iFrameStart = edgeOriginal->iFrameStart;
			ledge->iFrameEnd = edgeOriginal->iFrameEnd;
			ledge->lexUnit = edgeOriginal->lexUnit;
			ledge->fScoreAM = edgeOriginal->fScoreAM;
			ledge->fScoreLM = edgeOriginal->fScoreLM;
			ledge->fInsertionPenalty = edgeOriginal->fInsertionPenalty;
			
			ledge->nodePrev = NULL;
			ledge->nodeNext = NULL;
			ledge->edgePrev = NULL;
			ledge->edgeNext = NULL;	
			
			ledge->iPhones = edgeOriginal->iPhones;
			ledge->phoneAlignment = NULL;
			ledge->fPhoneAccuracy = NULL;
			
			ledge->iContextLeft = edgeOriginal->iContextLeft;
			ledge->iContextRight = edgeOriginal->iContextRight;
			
			ledge->bTouched = edgeOriginal->bTouched;
			ledge->bBestPath = false;		// IMP: the copy can't be in the best path too!!
		
			return ledge;
		}	
		
		// connect an edge
		static void connectEdge(LNode *lnode1, LEdge *ledge, LNode *lnode2) {
		
			ledge->nodePrev = lnode1;
			ledge->nodeNext = lnode2;
		
			ledge->edgePrev = lnode1->edgeNext;
			lnode1->edgeNext = ledge;
			ledge->edgeNext = lnode2->edgePrev;
			lnode2->edgePrev = ledge;	
		}
		
		// print the lattice (using the lattice container)
		void print();
		
		// print the lattice (by traversing the nodes starting from the initial node)
		void printTraverse();
		
		// print a lattice edge
		void print(LEdge *ledge);
		
		// print a list of edges
		void print(LLEdge &lLEdge);
		
		// print a lattice node
		void print(LNode *lnode);
		
		// print a phone alignment
		void print(LPhoneAlignment *phoneAlignment) {
			
			printf("%5d -> %5d\n",phoneAlignment->iStateBegin[0],phoneAlignment->iStateEnd[NUMBER_HMM_STATES-1]);	
		}
		
		// print the lattice properties
		void printProperties();
		
		// store the lattice into the given file
		void store(const char *strFile, unsigned char iFormat = FILE_FORMAT_BINARY);
		
		// store the lattice into the given file (text format)
		void storeTextFormat(const char *strFile);
		
		// store the lattice into the given file (binary format)
		void storeBinaryFormat(const char *strFile);
		
		// load the lattice from the given file
		void load(const char *strFile);	
		
		// load the lattice from the given file (binary format)
		void loadBinaryFormat(const char *strFile);
			
		// traverse the lattice from left to right and merge identical edges
		void forwardEdgeMerge();
			
		// traverse the lattice from right to left and merge identical edges
		void backwardEdgeMerge();
		
		// check whether two nodes have the same set of successor nodes reached
		// by equivalent arcs
		bool equivalentSuccessors(LNode *node1, LNode *node2) {
		
			// same number of successors?
			int iSuccessors = 0;
			for(LEdge *edge1 = node1->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
				++iSuccessors;
			}
			for(LEdge *edge2 = node2->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
				--iSuccessors;
			}
			if (iSuccessors != 0) {
				return false;
			}
		 
		 	// equivalent successors?
			for(LEdge *edge1 = node1->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
				bool bFound = false;
				for(LEdge *edge2 = node2->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
					if (equivalentEdges(edge1,edge2) && (edge1->nodeNext == edge2->nodeNext)) {
						bFound = true;
						break;
					}	
				}
				if (bFound == false) {
					return false;
				}
			}	
		
			return true;
		}

		// check whether two nodes have the same set of predecessors nodes reached
		// by equivalent arcs
		bool equivalentPredecessors(LNode *node1, LNode *node2) {
		
			// same number of predecessors?
			int iSuccessors = 0;
			for(LEdge *edge1 = node1->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
				++iSuccessors;
			}
			for(LEdge *edge2 = node2->edgePrev ; edge2 != NULL ; edge2 = edge2->edgeNext) {
				--iSuccessors;
			}
			if (iSuccessors != 0) {
				return false;
			}
		 
		 	// equivalent successors?
			for(LEdge *edge1 = node1->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
				bool bFound = false;
				for(LEdge *edge2 = node2->edgePrev ; edge2 != NULL ; edge2 = edge2->edgeNext) {
					if (equivalentEdges(edge1,edge2) && (edge1->nodePrev == edge2->nodePrev)) {
						bFound = true;
						break;
					}	
				}
				if (bFound == false) {
					return false;
				}
			}	
		
			return true;
		}

		// compare two edges for ordering
		static bool compareEdgesOrder(const LEdge *edge1, const LEdge *edge2) {
		
			if (edge1->iFrameStart == edge2->iFrameStart) {
				if (edge1->iFrameEnd == edge2->iFrameEnd) {
					return (edge1->lexUnit->iLexUnitPron <= edge2->lexUnit->iLexUnitPron);
				} else {
					return (edge1->iFrameEnd < edge2->iFrameEnd);
				}
			} else {
				return (edge1->iFrameStart < edge2->iFrameStart);
			}
		
			return false;
		}
		
		// compare two edges for equality
		static bool equivalentEdges(LEdge *ledge1, LEdge *ledge2) {
		
			if ((ledge1->iFrameStart == ledge2->iFrameStart) && 
				(ledge1->iFrameEnd == ledge2->iFrameEnd) && 
				(ledge1->lexUnit == ledge2->lexUnit)) {
				return true; 
			}
		
			return false;
		}
		
		const char *bool2str(bool b) {
			
			if (b) {
				return "yes";
			} else {
				return "no";
			}
		}
		
		const char* getStrConfidenceMeasure(unsigned char iConfidenceMeasure) {
		
			switch(iConfidenceMeasure) {
				case CONFIDENCE_MEASURE_NONE: {
					return CONFIDENCE_MEASURE_NONE_STR;
				}
				case CONFIDENCE_MEASURE_POSTERIORS: {
					return CONFIDENCE_MEASURE_POSTERIORS_STR;
				}
				case CONFIDENCE_MEASURE_ACCUMULATED: {
					return CONFIDENCE_MEASURE_ACCUMULATED_STR;
				}
				case CONFIDENCE_MEASURE_MAXIMUM: {
					return CONFIDENCE_MEASURE_MAXIMUM_STR;
				}
			}
			
			return NULL;
		}
		
		// return the edges in the lattice
		LEdge **getEdges(int *iEdges) {
			
			*iEdges = m_iEdges;
		
			return m_ledges;
		}
		
		// return the nodes in the lattice
		LNode **getNodes(int *iNodes) {
		
			*iNodes = m_iNodes;
		
			return m_lnodes;
		}
		
		// return the number of edges
		int getEdges() {
		
			return m_iEdges;
		}
		
		// return the initial node
		LNode *getNodeInitial() {
		
			return m_lnodeInitial;
		}
		
		// return the final node
		LNode *getNodeFinal() {
		
			return m_lnodeFinal;
		}		
		
		// return the number of frames
		int getFrames() {
		
			return m_iFrames;
		}
		
		// mark all the nodes as not touched
		void untouchNodes();	
		
		// mark all the edges as not touched
		void untouchEdges();
			
		// set all the acoustic-scores to the given value
		void setAMScores(float fValue);
			
		// set all the language-model-scores to the given value
		void setLMScores(float fValue);
		
		// forward-backward --------------------------------------------------------------------
		
		// compute forward/backward scores
		void computeForwardBackwardScores(float fScalingAM, float fScalingLM);
		
		// compute forward score for the edge
		void forward(LEdge *edge);
		
		// compute backward score for the edge
		void backward(LEdge *edge);
		
		// posterior-probabilities -------------------------------------------------------------
		
		// compute posterior probabilities from the hypothesis graph
		void computePosteriorProbabilities();
		
		// confidence estimation ---------------------------------------------------------------
		
		// compute posterior-probability based confidence estimates
		void computeConfidenceScore(unsigned char iConfidenceMeasure);
		
		// return whether two edges overlap in time
		bool overlap(LEdge *edge1, LEdge *edge2) {
		
			return !((edge1->iFrameStart > edge2->iFrameEnd) || (edge2->iFrameStart > edge1->iFrameEnd));
		}
		
		// return wheter a edge spans the given time frame
		bool spanTimeFrame(LEdge *edge, int iFrame) {
		
			if ((edge->iFrameStart <= iFrame) && (edge->iFrameEnd >= iFrame))  {
				return true;
			}
		
			return false;
		}		
		
		// checks ------------------------------------------------------------------------------
		
		// check lattice correctness
		void check();
		
		// lattice depth -----------------------------------------------------------------------
		
		// compute the lattice depth given the reference text
		// lattice depth is defined as the average number of lattice arcs containing words
		// that cross every time frame
		LatticeDepth *computeDepth(bool bIgnoreNonWords = false);
		
		// lattice Word Error Rate -------------------------------------------------------------
		
		// compute the Lattice Word Error Rate (also known as oracle)
		LatticeWER *computeWER(VLexUnit &vLexUnitReference, Mappings *mappings = NULL);
			
		// history item garbage collection: free-up unused history items
		void historyItemGarbageCollection();
		
		// initializes a new history item to the given values
		int newLWERHistoryItem(int iPrev, unsigned char iAlignmentEvent, LEdge *edge, int iReference) {
		
			if (m_iHistoryItemAvailable == -1) {
				historyItemGarbageCollection();
			}
			assert(m_iHistoryItemAvailable != -1);
			
			LWERHistoryItem *historyItem = m_historyItems+m_iHistoryItemAvailable;
			m_iHistoryItemAvailable = (m_historyItems+m_iHistoryItemAvailable)->iPrev;
		
			historyItem->iPrev = iPrev;
			historyItem->iAlignmentEvent = iAlignmentEvent;
			historyItem->edge = edge;	
			historyItem->iReference = iReference;
		
			return (int)(historyItem-m_historyItems);
		}
		
		// delete a history item
		void deleteLWERHistoryItem(int iHistoryItem) {
		
			(m_historyItems+iHistoryItem)->iPrev = m_iHistoryItemAvailable;
			m_iHistoryItemAvailable = iHistoryItem;
		}
		
		// activates a token
		void activateToken(LWERToken *tokenPrev, LEdge *edge, unsigned char iAlignmentEvent) {
		
			LWERToken *token = &m_tokenNext[m_iActiveTokensNext];
			
			if (edge != NULL) {
				token->node = edge->nodeNext;
			} else {
				token->node = tokenPrev->node;
			}
			token->fScore = tokenPrev->fScore;
			switch(iAlignmentEvent) {
				case ALIGNMENT_EVENT_INSERTION: {
					token->fScore += ALIGNMENT_PENALTY_INSERTION;
					token->iReference = tokenPrev->iReference;
					break;
				}
				case ALIGNMENT_EVENT_DELETION: {
					token->fScore += ALIGNMENT_PENALTY_DELETION;
					token->iReference = tokenPrev->iReference+1;
					break;
				}
				case ALIGNMENT_EVENT_SUBSTITUTION: {
					token->fScore += ALIGNMENT_PENALTY_SUBSTITUTION;
					token->iReference = tokenPrev->iReference+1;
					break;
				}
				case ALIGNMENT_EVENT_CORRECT: {
					token->iReference = tokenPrev->iReference+1;
					break;
				}
				case ALIGNMENT_EVENT_SKIP: {
					token->iReference = tokenPrev->iReference;
					break;
				}
				default: {
					assert(0);
				}
			}
			token->iHistoryItem = newLWERHistoryItem(tokenPrev->iHistoryItem,iAlignmentEvent,edge,token->iReference);
			assert(token->node != NULL);
			
			// try to recombine
			bool bRecombined = false;
			for(LWERToken *tokenAux = m_tokenNode[token->node->iNode] ; tokenAux != NULL ; tokenAux = tokenAux->prev) {
				assert(tokenAux->node != NULL);
				if ((tokenAux->node == token->node) && (tokenAux->iReference == token->iReference)) {
					if (tokenAux->fScore > token->fScore) {
						tokenAux->fScore = token->fScore;
						deleteLWERHistoryItem(tokenAux->iHistoryItem);
						tokenAux->iHistoryItem = token->iHistoryItem;
					} else {
						deleteLWERHistoryItem(token->iHistoryItem);
					}
					bRecombined = true;
					break;
				}
			}	
			if (bRecombined == false) {
				token->prev = m_tokenNode[token->node->iNode];
				m_tokenNode[token->node->iNode] = token;
				assert(token->node != NULL);
				++m_iActiveTokensNext;
				assert(m_iActiveTokensNext < m_iMaxActiveTokens);
			}
		}
		
		// reset the object values
		static void reset(LatticeDepth *latticeDepth) {
		
			latticeDepth->iFrames = 0;
			latticeDepth->iFramesEdges = 0;
			latticeDepth->fDepth = 0.0;
		}
		
		// accumulate
		static void add(LatticeDepth *latticeDepth1, LatticeDepth *latticeDepth2) {
		
			latticeDepth1->iFramesEdges += latticeDepth2->iFramesEdges;
			latticeDepth1->iFrames += latticeDepth2->iFrames;
			latticeDepth1->fDepth = ((float)latticeDepth1->iFramesEdges)/((float)latticeDepth1->iFrames);
		}
		
		// print 
		static void print(LatticeDepth *latticeDepth) {
		
			printf("-- lattice stats --------------\n");
			printf("# frames edges:       %8d\n",latticeDepth->iFramesEdges);
			printf("# frames :            %8d\n",latticeDepth->iFrames);
			printf("depth:                %8.2f\n",latticeDepth->fDepth);
			printf("-------------------------------\n");
		}	
		
		// reset the object values
		static void reset(LatticeWER *latticeWER) {
		
			latticeWER->iWordsReference = 0;
			latticeWER->iInsertions = 0;
			latticeWER->iDeletions = 0;
			latticeWER->iSubstitutions = 0;
			latticeWER->iCorrect = 0;
			latticeWER->iErrors = 0;
			latticeWER->iOOV = 0;
			latticeWER->fWER = 0.0;
		}
		
		// accumulate 
		static void add(LatticeWER *latticeWER1, LatticeWER *latticeWER2) {
		
			latticeWER1->iWordsReference += latticeWER2->iWordsReference;
			latticeWER1->iInsertions += latticeWER2->iInsertions;
			latticeWER1->iDeletions += latticeWER2->iDeletions;
			latticeWER1->iSubstitutions += latticeWER2->iSubstitutions;
			latticeWER1->iCorrect += latticeWER2->iCorrect;
			latticeWER1->iErrors += latticeWER2->iErrors;
			latticeWER1->iOOV += latticeWER2->iOOV; 
			latticeWER1->fWER = (latticeWER1->iErrors*100.0f)/((float)latticeWER1->iWordsReference);
		}
		
		// print 
		static void print(LatticeWER *latticeWER) {
		
			printf("-- wer depth stats ------------\n");
			printf("events:\n");
			printf(" -sub:      %8d\n",latticeWER->iSubstitutions);
			printf(" -del:      %8d\n",latticeWER->iDeletions);
			printf(" -ins:      %8d\n",latticeWER->iInsertions);
			printf(" -correct:  %8d\n",latticeWER->iCorrect);
			printf(" -errors:   %8d\n",latticeWER->iErrors);
			printf(" -OOV:      %8d (%8.2f%%)\n",latticeWER->iOOV,(latticeWER->iOOV*100.0)/((float)latticeWER->iWordsReference));
			printf("WER:       %8.2f%%\n",latticeWER->fWER);
			printf("Accuracy:  %8.2f%%\n",100.0-latticeWER->fWER);
			printf("-------------------------------\n");
		}
		
		// lattice rescoring -------------------------------------------------------------------
		
		// rescore the lattice using the given rescoring method (Dijkstra)
		BestPath *rescore(const char *strRescoringMethod);
		
		// return the edge loglikelihood for likelihood-based rescoring
		float edgeLogLikelihood(LEdge *edge) {
		
			float fLikelihood = (edge->lexUnit->fInsertionPenalty+edge->fScoreAM)*m_fAMScalingFactor;
			fLikelihood += edge->fScoreLM*m_fLMScalingFactor;	
			
			return fLikelihood;
		}
		
		// return the edge log-posterior probability for pp-based rescoring
		float edgeLogPP(LEdge *edge) {
		
			return log(edge->fPP);
		}
		
		// set scaling factors
		void setScalingFactors(float fScaleAM, float fScaleLM) {
		
			m_fAMScalingFactor = fScaleAM;
			m_fLMScalingFactor = fScaleLM;
		}
		
		// return the best path in the lattice (must be already marked)
		BestPath *getBestPath();
		
		// return the alignment of the best path in the lattice (must be already marked)
		VLPhoneAlignment *getBestPathAlignment();
		
		// add a path to the lattice
		void addPath(Alignment *alignment, bool bIsBest = false);
		
		// clear the best path in the lattice
		void clearBestPath();
		
		// Phone-level accuracy (neded by some discriminative training methods: MWE, MPE, bMMI, etc) ------------
		
		// compute the phone accuracy for phones within each arc given the reference word sequence
		void computePhoneAccuracy(VLPhoneAlignment &vLPhoneAlignment, 
			bool bSetSilenceToZero = true, bool bSetFillersToZero = true);	
		
		// return the proportion of the reference phone that overlaps with the given hypothesis phone
		float getOverlapProportion(LPhoneAlignment *phoneAlignmentRef, LPhoneAlignment *phoneAlignmentHyp) {
		
			int iBeg1 = phoneAlignmentRef->iStateBegin[0];
			int iEnd1 = phoneAlignmentRef->iStateEnd[NUMBER_HMM_STATES-1];
			int iBeg2 = phoneAlignmentHyp->iStateBegin[0];
			int iEnd2 = phoneAlignmentHyp->iStateEnd[NUMBER_HMM_STATES-1];
			
			// no overlap?
			if ((iBeg1 > iEnd2) || (iBeg2 > iEnd1)) {	
				return 0.0;
			}
			if (iBeg1 >= iBeg2) {
				return ((float)(min(iEnd1,iEnd2)-iBeg1+1))/((float)(iEnd1-iBeg1+1));
			} else {
				return ((float)(min(iEnd1,iEnd2)-iBeg2+1))/((float)(iEnd1-iBeg1+1));
			}
		}
		
		// HMM-id marking ----------------------------------------------------------------------
		
		// mark each of the edges in the lattice with HMM-information
		void hmmMarking(HMMManager *hmmManager);
		
		// allocates memory for a blank context
		unsigned char *newBlankContext() {
			
			int iSize = max(m_iContextSizeCW,m_iContextSizeWW);
			unsigned char *iContext = new unsigned char[iSize];
			for(int i=0 ; i < iSize ; ++i) {
				iContext[i] = m_iPhoneContextPadding;
			}
		
			return iContext;
		}
		
		// propagate left-context from one edge to the next
		// note: destination context must be already filled with padding symbols
		void shiftLeftContext(unsigned char *iContextLeftSrc, unsigned char iPosition, 
			unsigned char iPhone, LexUnit *lexUnit, unsigned char *iContextLeftDest) {
		
			// cross-word context propagation 
			if ((iPosition == WITHIN_WORD_POSITION_END) || 
				(iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
					
				iContextLeftDest[m_iContextSizeCW-1] = iPhone;
				for(int i=1 ; i < min((int)m_iContextSizeCW,(int)(lexUnit->vPhones.size())) ; ++i) {
					iContextLeftDest[m_iContextSizeCW-i-1] = iContextLeftSrc[m_iContextSizeCW-i];
				}
			}
			// within-word context propagation
			else {
				iContextLeftDest[m_iContextSizeWW-1] = iPhone;
				// word-initial phone to second phone: do nothing else
				if (iPosition == WITHIN_WORD_POSITION_START) {
					// do nothing
				} 
				// word-internal to word-internal or word-end
				else {
					for(int i=0 ; i < m_iContextSizeWW-1 ; ++i) {
						iContextLeftDest[i] = iContextLeftSrc[i+1]; 
					}
				}
			}
		}

		// propagate right-context from one edge to the previous edge
		// note: destination context must be already filled with padding symbols
		void shiftRightContext(unsigned char *iContextRightSrc, unsigned char iPosition, 
			unsigned char iPhone, LexUnit *lexUnit, unsigned char *iContextRightDest) {
		
			// cross-word context propagation 
			if ((iPosition == WITHIN_WORD_POSITION_START) || 
				(iPosition == WITHIN_WORD_POSITION_MONOPHONE)) {
					
				iContextRightDest[0] = iPhone;
				for(int i=1 ; i < min((int)m_iContextSizeCW,(int)(lexUnit->vPhones.size())) ; ++i) {
					iContextRightDest[i] = iContextRightSrc[i-1];
				}
			}
			// within-word context propagation
			else {
				iContextRightDest[0] = iPhone;
				// word-final phone to phone before last phone: do nothing else
				if (iPosition == WITHIN_WORD_POSITION_END) {
				} 
				// word-internal to word-internal or word-end
				else {
					for(int i=1 ; i < m_iContextSizeWW ; ++i) {
						iContextRightDest[i] = iContextRightSrc[i-1]; 
					}
				}
			}
		}
				
		// copy context
		unsigned char *copyContext(unsigned char *iContext) {
		
			int iSize = max(m_iContextSizeCW,m_iContextSizeWW);
			unsigned char *iContextCopy = new unsigned char[iSize];
			for(int i=0 ; i < iSize ; ++i) {
				iContextCopy[i] = iContext[i];
			}
			
			return iContextCopy;
		}
		
		// print a context (debugging)
		void printContext(unsigned char *iContext) {
		
			int iSize = max(m_iContextSizeCW,m_iContextSizeWW);
			for(int i=0 ; i < iSize ; ++i) {
				if (iContext[i] != m_iPhoneContextPadding) {
					printf("%s ",m_phoneSet->getStrPhone(iContext[i]));
				} else {
					printf("<> ");
				}
			}
			printf("\n");
		}
		
		// attach lm-probabilities to the edges of the lattice
		// note: lattice expansion may be needed in order for each edge to have a unique
		//       word context
		void attachLMProbabilities(LMManager *lmManager);
		
		// attach insertion penalties
		void attachInsertionPenalty(LexiconManager *lexiconManager);	
		
		// compute lattice likelihood
		double getLikelihood();
		
};

};	// end-of-namespace

#endif
