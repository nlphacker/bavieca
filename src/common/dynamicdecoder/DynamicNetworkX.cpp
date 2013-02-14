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


#include "DynamicNetworkX.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
DynamicNetworkX::DynamicNetworkX(DNode *nodes, int iNodes, DArc *arcs, int iArcs, int *iLATree, int iLANodes)
{
	m_nodes = nodes;
	m_iNodes = iNodes;
	m_arcs = arcs;
	m_iArcs = iArcs;
	
	// look-ahead
	m_iLATree = iLATree;
	m_iLANodes = iLANodes;
	
	m_fIP = NULL;
	m_iIPSize = -1;
}

// destructor
DynamicNetworkX::~DynamicNetworkX()
{
	delete [] m_arcs;
	delete [] m_nodes;
		
	if (m_iLATree) {
		delete [] m_iLATree;
	}	
	if (m_fIP != NULL) {
		delete [] m_fIP;
	}
}

// print the network
void DynamicNetworkX::print(LexiconManager *lexiconManager) {

	printf("# nodes: %12d\n",m_iNodes);
	printf("# arcs:  %12d\n",m_iArcs);
	printf("nodes ----------------------------------\n");
	for(int i=0 ; i < m_iNodes ; ++i) {
		if (m_nodes[i].bWordEnd) {
			printf("node [%d] %d (word-end)\n",i,m_nodes[i].iDepth);
		} else {
			printf("node [%d] %d\n",i,m_nodes[i].iDepth);
		}
		DArc *arcEnd = m_arcs+(m_nodes[i+1].iArcNext);
		for(DArc *arc = m_arcs+m_nodes[i].iArcNext ; arc != arcEnd ; ++arc) {
			if (arc->iType == ARC_TYPE_NULL) {
				printf("null: -> [%d] %d\n",arc->iNodeDest,m_nodes[arc->iNodeDest].iDepth);
			} else if (arc->iType == ARC_TYPE_HMM) {
				printf("hmm: %8d -> [%d] %d\n",
					arc->state->getId(),arc->iNodeDest,m_nodes[arc->iNodeDest].iDepth);
			} else if (arc->iType == ARC_TYPE_WORD) {
				printf("word: %s -> [%d] %d\n",
					lexiconManager->getStrLexUnitPron(arc->lexUnit->iLexUnitPron),
					arc->iNodeDest,m_nodes[arc->iNodeDest].iDepth);
			}	
		}
	}	
}

};	// end-of-namespace

