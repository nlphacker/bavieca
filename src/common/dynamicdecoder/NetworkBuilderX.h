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


#ifndef NETWORKBUILDERX_H
#define NETWORKBUILDERX_H

#include "DynamicNetworkX.h"
#include "Global.h"
#include "HMMManager.h"
#include "LexiconManager.h"

#if defined __linux__ || defined __APPLE__ || __MINGW32__
#include <tr1/unordered_map>
#elif _WIN32
#include <hash_map>
#else 
	#error "unsupported platform"
#endif

using namespace std;

#include <vector>
#include <deque>
#include <list>

namespace Bavieca {

class PhoneSet;

#define TEMPNODE_STATE_QUEUED			0
#define TEMPNODE_STATE_PROCESSED		1
#define TEMPNODE_STATE_UNSEEN			2

struct _ArcTempX;

typedef struct _NodeTempX {
   unsigned char iType;					// this is meaningful in the case of dummy nodes (field *model != NULL) 
   unsigned char iDepth;				// node depth in the network 
	bool bWordEnd;							// whether the node corresponds to a word-end
	char iIPIndex;							// insertion penalty (index within the array)	
   int iNode;								// node id
   HMMStateDecoding *state;			// HMM-state
   list<_ArcTempX*> lArcNext;			// pointers to the next nodes in the network
   list<_ArcTempX*> lArcPrev;			// pointers to the previous nodes in the network	
   int iFI;									// index in the table of FI-nodes
} NodeTempX;

typedef struct _ArcTempX {
	unsigned char iType;					// HMM-state or Word-Identity
	HMMStateDecoding *state;			// HMM-state
	LexUnit *lexUnit;						// lexical unit
	NodeTempX *nodeSource;				// prev node
	NodeTempX *nodeDest;					// next node
	int iLANode;							// index within the look-ahead tree (which is stored as an array)
} ArcTempX;

typedef vector<ArcTempX*> VArcTempX; 
typedef deque<ArcTempX*> DArcTempX; 
typedef list<ArcTempX*> LArcTempX; 
typedef map<ArcTempX*,bool> MArcTempX;

typedef vector<NodeTempX*> VNodeTempX; 
typedef list<NodeTempX*> LNodeTempX; 
typedef map<NodeTempX*,bool> MNodeTempX;

typedef struct {
	int iMerged;				// entry to with the current entry was merged
	NodeTempX *node;  		// FI-node
	unsigned char *iKey;		// FI-key
} NodeMergeInfo;


// ad-hoc functions to use a phonetic context as the key in a hash_map data structure
struct MContextHashFunctions {

#if defined __linux__ || defined __APPLE__ || defined __MINGW32__

	// comparison function (used for matching, comparison for equality)
	bool operator()(const unsigned char *contextUnit1, const unsigned char *contextUnit2) const {
	
		int i=0;
		for( ; contextUnit1[i] != UCHAR_MAX ; ++i) {
			assert(contextUnit2[i] != UCHAR_MAX);
			if (contextUnit1[i] != contextUnit2[i]) {
				return false;
			}
		}
		assert(contextUnit2[i] == UCHAR_MAX);
		
		return true;
	}

#elif _MSC_VER

	static const size_t bucket_size = 4;
	static const size_t min_buckets = 8;

	// comparison function (used to order elements)
	bool operator()(const unsigned char *contextUnit1, const unsigned char *contextUnit2) const {
	
		int i=0;
		for( ; contextUnit1[i] != UCHAR_MAX ; ++i) {
			assert(contextUnit2[i] != UCHAR_MAX);
			if (contextUnit1[i] == contextUnit2[i]) {
				continue;
			}
			return (contextUnit1[i] < contextUnit2[i]);
		}
		assert(contextUnit2[i] == UCHAR_MAX);
		// a == b, then (a < b) and (b < a) must be false
		return false;
	}	

#endif
	
	// hash function
	size_t operator()(const unsigned char *contextUnit) const {
	
		unsigned int iAcc = 0;
		unsigned int iAux = 0;
		for(int i=0 ; contextUnit[i] != UCHAR_MAX ; ++i) {
			if (i <= 3) {	
				iAcc <<= (8*i);
				iAcc += contextUnit[i];
			} else {
				iAux = contextUnit[i];
				iAux <<= (8*(i%4));
				iAcc ^= iAux;
			}
		}
	
		return iAcc;
	}
};

// structure for easy management of different contexts (keeps existing contexts)
#if defined __linux__ || defined __APPLE__ || __MINGW32__
typedef std::tr1::unordered_map<unsigned char*,bool,MContextHashFunctions,MContextHashFunctions> MContextBool;
typedef std::tr1::unordered_map<unsigned char*,NodeTempX*,MContextHashFunctions,MContextHashFunctions> MContextNode;
typedef std::tr1::unordered_map<unsigned char*,vector< pair<unsigned char*,NodeTempX*> >,MContextHashFunctions,MContextHashFunctions> MContextV;
typedef std::tr1::unordered_map<int,NodeTempX*> MIntNodeTempX;
#elif _WIN32
typedef hash_map<unsigned char*,bool,MContextHashFunctions> MContextBool;
typedef hash_map<unsigned char*,NodeTempX*,MContextHashFunctions> MContextNode;
typedef hash_map<unsigned char*,vector< pair<unsigned char*,NodeTempX*> >,MContextHashFunctions> MContextV;
typedef hash_map<int,NodeTempX*> MIntNodeTempX;
#endif

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class NetworkBuilderX {

	private:
	
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManager;
		LexiconManager *m_lexiconManager;
		
		// network properties
   	int m_iFIDepth;
   	int m_iWIDepth;
   	int m_iFODepth;
   	
   	// node depth
   	int m_iDepthFI;
   	int m_iDepthMI;
   	int m_iDepthFO;
   	int m_iDepthWI;
   	
   	// temporal network building
   	NodeTempX *m_nodeTempRoot;
   	
   	int m_iNodesTempID;
   	
   	int m_iNodesTemp;									// total number of temporal nodes created
   	int m_iNodesTempRoot;
   	int m_iNodesTempFI;	
   	int m_iNodesTempMI;
   	int m_iNodesTempFO;	
   	int m_iNodesTempWord;
   	int m_iNodesTempHMM;
   	int m_iNodesTempNull;
   	int *m_iNodesDepth;	
   		
   	int m_iArcsTemp;
   	int m_iArcsNull;
   	int m_iArcsHMM;
   	int m_iArcsWord;
   	
   	// language model look-ahead
	   int m_iLANodes;
   	int *m_iLATree;
   	
   		   	
   	// return a copy of the given context
   	unsigned char *copyContext(unsigned char *iContext, int iSize) {
   	
   		unsigned char *iContextCopy = new unsigned char[iSize+1];
   		int i=0;
   		for( ; i < iSize+1 ; ++i) {
   			iContextCopy[i] = iContext[i];
   		}	
   		
   		return iContextCopy;
   	}
	
		// auxiliar function to create an arc		
		ArcTempX *newArcTempX(NodeTempX *nodeSource, NodeTempX *nodeDest, 
			unsigned char iType, HMMStateDecoding *state, LexUnit *lexUnit) {
			
			ArcTempX *arc = new ArcTempX;
			++m_iArcsTemp;
			
			assert(nodeSource != NULL);
			assert(nodeDest != NULL);
			
			switch(iType) {
				case ARC_TYPE_NULL: {
					++m_iArcsNull;
					arc->lexUnit = NULL;
					break;
				}
				case ARC_TYPE_WORD: {
					assert(lexUnit != NULL);
					arc->lexUnit = lexUnit;
					++m_iArcsWord;
					break;
				}
				default: {
					assert(0);
				}
			}
		
			arc->iType = iType;
			arc->nodeSource = nodeSource;	
			arc->nodeDest = nodeDest;	
			
			// make connections
			nodeSource->lArcNext.push_back(arc);
			nodeDest->lArcPrev.push_back(arc);
			
			return arc;
		}
		
		// auxiliar function to delete an arc		
		void deleteArcTempX(ArcTempX *arcTemp) {
			
			--m_iArcsTemp;
			
			switch(arcTemp->iType) {
				case ARC_TYPE_NULL: {
					--m_iArcsNull;
					break;
				}
				case ARC_TYPE_WORD: {
					--m_iArcsWord;
					break;
				}
				default: {
					assert(0);
				}
			}
		
			delete arcTemp;
		}
		
		// return whether two nodes are connected
		bool areConnected(NodeTempX *nodeSource, NodeTempX *nodeDest) {
		
			for(LArcTempX::iterator it = nodeSource->lArcNext.begin() ; it != nodeSource->lArcNext.end() ; ++it) {
				if ((*it)->nodeDest == nodeDest) {
					return true;
				}
			}
				
			return false;
		}
		
		// auxiliar function to create a node
		NodeTempX *newNodeTempX(unsigned char iType, unsigned char iDepth, HMMStateDecoding *state = NULL, char iIPIndex = -1) {
			
			NodeTempX *nodeTemp = new NodeTempX;
			nodeTemp->iType = iType;
			nodeTemp->iDepth = iDepth;	
			nodeTemp->iNode = m_iNodesTempID++;
			nodeTemp->state = state;
			nodeTemp->bWordEnd = false;
			nodeTemp->iIPIndex = iIPIndex;
			nodeTemp->iFI = -1;
			++m_iNodesTemp;
			
			// stats
			switch(iType) {
				case NODE_TYPE_ROOT: {
					m_iNodesTempRoot++;
					break;
				}
				case NODE_TYPE_FI: {
					m_iNodesTempFI++;
					break;
				}
				case NODE_TYPE_MI: {
					m_iNodesTempMI++;
					break;
				}
				case NODE_TYPE_FO: {
					m_iNodesTempFO++;
					break;
				}
				case NODE_TYPE_HMM: {
					m_iNodesTempHMM++;
					break;
				}
				case NODE_TYPE_WORD: {
					m_iNodesTempWord++;
					break;
				}
				case NODE_TYPE_NULL: {
					m_iNodesTempNull++;
					break;
				}
				default: {
					assert(0);
				}
			}
		
			return nodeTemp;
		}
		
		// auxiliar function to create a node
		void deleteNodeTempX(NodeTempX *nodeTemp) {
			
			--m_iNodesTemp;
			
			// stats
			switch(nodeTemp->iType) {
				case NODE_TYPE_ROOT: {
					m_iNodesTempRoot--;
					break;
				}
				case NODE_TYPE_FI: {
					m_iNodesTempFI--;
					break;
				}
				case NODE_TYPE_MI: {
					m_iNodesTempMI--;
					break;
				}
				case NODE_TYPE_FO: {
					m_iNodesTempFO--;
					break;
				}
				case NODE_TYPE_HMM: {
					m_iNodesTempHMM--;
					break;
				}
				case NODE_TYPE_WORD: {
					m_iNodesTempWord--;
					break;
				}
				case NODE_TYPE_NULL: {
					m_iNodesTempNull--;
					break;
				}
				default: {
					assert(0);
				}
			}
		
			delete nodeTemp;
		}
		
		// push the word-labels to the beginning of the network
		void pushWordLabels();	
		
		// compact the network
		DynamicNetworkX *compact();
		
		// merge nodes forward
		void mergeForward();
		
		// merge nodes forward starting from the given group of FO-nodes
		void mergeForwardFO(VNodeTempX &vNodeTempFO);
		
		// check that forward-compacting was done right (debugging)
		void checkForward();
		
		// merge nodes backward
		void mergeBackward();
		
		// merge nodes backward starting from the given group of MI-nodes
		void mergeBackwardMI(VNodeTempX &vNodeTempMI, NodeMergeInfo *fiNodeMergeInfo);
		
		// merge nodes backward starting from the given group of FI-nodes
		void mergeBackwardFI(VNodeTempX &vNodesFI, NodeMergeInfo *foNodeMergeInfo);
			
		// merge nodes backward starting from the given group of MI-nodes
		void mergeBackwardWordNodes(VNodeTempX &vNodeWord);
		
		// remove FI-nodes used to build the network
		void removeFINodes();
		
		// remove MI nodes used to build the network
		void removeMINodes(MContextNode &mContextBeg);		
		
		// print network stats
		void print();
		
		// print a context (auxiliar function)
		void print(unsigned char *iContext) {
		
			for(int i=0 ; iContext[i] != UCHAR_MAX ; ++i) {
				if (iContext[i] == m_phoneSet->size()) {
				   printf("<-> ");
				} else {
					printf("%3s ",m_phoneSet->getStrPhone(iContext[i]));
				}
			}
			printf("\n");
		}
		
		// print a context (auxiliar function)
		void print(unsigned char *iContext, int iSize) {
		
			for(int i=0 ; i < iSize ; ++i) {				
				if (iContext[i] == m_phoneSet->size()) {
				   printf("<-> ");
				} else {
					printf("%3s ",m_phoneSet->getStrPhone(iContext[i]));
				}
			}
			printf("\n");
		}
		
		// return whether two nodes are equivalent
		inline bool equivalentNodes(NodeTempX *node1, NodeTempX *node2) {
		
			if (node1->iType != node2->iType) {
				return false;
			}
			switch(node1->iType) {
				case NODE_TYPE_HMM: {
					return (node1->state == node2->state);
				}
				case NODE_TYPE_FI: {		
					return true;
				}
				case NODE_TYPE_FO: {
					return true;
				}
				case NODE_TYPE_MI: {
					return true;
				}
				case NODE_TYPE_ROOT: {
					return true;
				}
				case NODE_TYPE_WORD: {
					return true;
				}
				case NODE_TYPE_NULL: {
					return true;
				}
				default: {
					assert(0);
				}
			}
			
			return true;
		}
		
		// auxiliar function to delete an arc		
		inline bool equivalentArcs(ArcTempX *arc1, ArcTempX *arc2) {
		
			if (arc1->iType != arc2->iType) {
				return false;
			}
			
			if (arc1->iType == ARC_TYPE_NULL) {
				return true;
			} else {
				assert(arc1->iType == ARC_TYPE_WORD);
				return (arc1->lexUnit == arc2->lexUnit);	
			}
		}		
		
		// return wether two nodes have the same set of predecessors
		inline bool samePredecessors(NodeTempX *node1, NodeTempX *node2) {
		
			if (node1->lArcPrev.size() != node2->lArcPrev.size()) {
				return false;
			}
			
			for(LArcTempX::iterator it = node1->lArcPrev.begin() ; it != node1->lArcPrev.end() ; ++it) {
				bool bFound = false;
				for(LArcTempX::iterator jt = node2->lArcPrev.begin() ; jt != node2->lArcPrev.end() ; ++jt) {
					if ((*jt)->nodeSource == (*it)->nodeSource) {
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

		// return wether two nodes have the same set of successors
		inline bool sameSuccessors(NodeTempX *node1, NodeTempX *node2) {
		
			if (node1->lArcNext.size() != node2->lArcNext.size()) {
				return false;
			}
			
			for(LArcTempX::iterator it = node1->lArcNext.begin() ; it != node1->lArcNext.end() ; ++it) {
				bool bFound = false;
				for(LArcTempX::iterator jt = node2->lArcNext.begin() ; jt != node2->lArcNext.end() ; ++jt) {
					if ((*jt)->nodeDest == (*it)->nodeDest) {
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
		
		// build the look-ahead tree from the network and transform it to a topologically sorted array
		// - leaves correspond to words in the vocabulary
		// - positions in the array keep the index of the parent node
		int *buildLMLookAheadTree(int *iNodes);
	
	public:

		// constructor
		NetworkBuilderX(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager);

		// destructor
		~NetworkBuilderX();
		
		// build the network (search)
		DynamicNetworkX *build();

};

};	// end-of-namespace

#endif
