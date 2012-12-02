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


#ifndef DYNAMICNETWORKX_H
#define DYNAMICNETWORKX_H

#include "Global.h"
#include "HMMManager.h"
#include "LexiconManager.h"

using namespace std;

#include <list>
#include <map>
#include <vector>

namespace Bavieca {

class TimeUtils;

// node type
#define NODE_TYPE_ROOT		0	
#define NODE_TYPE_FI			1			// temporal network only
#define NODE_TYPE_MI			2			// temporal network only
#define NODE_TYPE_FO			3			// temporal network only
#define NODE_TYPE_HMM		4			// temporal network only
#define NODE_TYPE_WORD		5
#define NODE_TYPE_NULL		6

// arc type
#define ARC_TYPE_NULL		0
#define ARC_TYPE_HMM			1
#define ARC_TYPE_WORD		2

// network node: final node
typedef struct _DNode {
	unsigned char iType;						// node type (hmm|null)
	unsigned char iDepth;					// node depth
	bool bWordEnd;								// whether the node can act as an end-of word
	char iIPIndex;								// insertion penalty (index within the array)
	int iArcNext;								// first index of successor acrs
	// token activation fields
	int iActiveTokensCurrent;				// # active tokens
	int iActiveTokensCurrentBase;			// array of active tokens (current time frame)
	int iActiveTokensNext;					// # active tokens
	int iActiveTokensNextBase;				// array of active tokens (next time frame)	
} DNode;

// network arc:
typedef struct _DArc {
	unsigned char iType;						// arc type
	union {
		HMMStateDecoding *state;			// hmm-state
		LexUnit *lexUnit;						// lexical-unit
	};
	int iNodeDest;								// index of destination node
	int iLANode;								// index within the look-ahead tree (which is stored as an array)
} DArc;


/**
	@author daniel <dani.bolanos@gmail.com>
*/
class DynamicNetworkX {

	private:

		DNode *m_nodes;
		int m_iNodes;
		DArc *m_arcs;
		int m_iArcs;
		
		// look ahead
		int *m_iLATree;				// look ahead tree
		int m_iLANodes;				// # look ahead nodes
		
		// insertion-penalty values attached to entry arcs
		float *m_fIP;
		int m_iIPSize;

	public:

		// constructor
		DynamicNetworkX(DNode *nodes, int iNodes, DArc *arcs, int iArcs, int *iLATree, int iLANodes);

		// destructor
		~DynamicNetworkX();
		
		// return the root node
		inline DNode *getRootNode() {
		
			return &m_nodes[0];
		}
		
		// return the array of arcs
		inline DArc *getArcs(int *iArcs) {
		
			*iArcs = m_iArcs;
		
			return m_arcs;
		}

		// return the array of nodes
		inline DNode *getNodes(int *iNodes) {
		
			*iNodes = m_iNodes;
		
			return m_nodes;
		}
		
		// set the array of insertion-penalties for entry arcs in the network
		inline void setIP(float *fIP, int iSize) {
			
			m_fIP = fIP;
			m_iIPSize = iSize; 
		}
		
		// return the insertion-penalty for the given index
		inline float getIP(char iIndex) {
		
			assert((iIndex >= 0) && (iIndex < m_iIPSize));
			return m_fIP[(unsigned char)iIndex];
		}
		
		// print the network
		void print(LexiconManager *lexiconManager);
		
		// get the look ahead tree
		inline int *getLMLookAheadTree(int *iNodes) {
		
			*iNodes = m_iLANodes;
			return m_iLATree;
		}

};

};	// end-of-namespace

#endif
