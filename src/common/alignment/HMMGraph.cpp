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


#include "HMMGraph.h"
#include "TimeUtils.h"

#include <iomanip>

namespace Bavieca {

// constructor
HMMGraph::HMMGraph(PhoneSet *phoneSet, LexiconManager *lexiconManager, HMMManager *hmmManagerEstimation,
	HMMManager *hmmManagerUpdate, bool bMultiplePronunciations, VLexUnit &vLexUnitOptional, bool bMergeAcrossDiffLexUnits)
{
	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
	m_hmmManagerEstimation = hmmManagerEstimation;
	m_hmmManagerUpdate = hmmManagerUpdate;
	m_bMultiplePronunciations = bMultiplePronunciations;
	m_vLexUnitOptional = vLexUnitOptional;
	m_bMergeAcrossDiffLexUnits = bMergeAcrossDiffLexUnits; 
}

// destructor
HMMGraph::~HMMGraph()
{
}

// create a graph from a sequence of lexical units
// - handles multiple pronunciations
// - handles optional symbols (typically silence+fillers)
FBNodeHMM **HMMGraph::create(VLexUnit &vLexUnitTranscription, int *iNodes, int *iEdges, 
			FBNodeHMM **nodeHMMInitial, FBNodeHMM **nodeHMMFinal) {
	
	//double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	if (vLexUnitTranscription.empty() == true) {
		return NULL;
	}
	
	// (1) create the lex-unit graph
	FBNodeLexUnit *nodeLexUnitInitial = NULL;
	int iNodesLexUnit = -1;
	// multiple pronunciations
	if (m_bMultiplePronunciations == true) {
		VLexUnitX vLexUnitX;
		m_lexiconManager->map(vLexUnitTranscription,vLexUnitX);
		nodeLexUnitInitial = create(vLexUnitX,m_vLexUnitOptional,&iNodesLexUnit);
	}
	// given pronunciation 
	else {
		nodeLexUnitInitial = create(vLexUnitTranscription,m_vLexUnitOptional,&iNodesLexUnit);
	}
	//print(nodeLexUnitInitial);

	// (2) create the phone-graph (including alternative pronunciations and optional symbols)
	int iNodesPhone = -1;
	FBNodePhone *nodePhoneFinal = NULL;
	FBNodePhone *nodePhoneInitial = create(nodeLexUnitInitial,iNodesLexUnit,&iNodesPhone,&nodePhoneFinal);
	destroy(nodeLexUnitInitial,iNodesLexUnit);
	//print(nodePhoneInitial);
	//destroy(nodePhoneInitial,iNodesPhone);
	
	// (3) propagate context across the graph (if needed)
	if (m_hmmManagerUpdate->areAccumulatorsLogical() == true) {
		m_iContextSizeWW = m_hmmManagerUpdate->getContextSizeAccumulators();
		m_iContextSizeCW = m_hmmManagerUpdate->getContextSizeAccumulatorsCW();
		m_iPhoneContextPadding = m_phoneSet->size();
	} else {		
		m_iContextSizeWW = m_hmmManagerUpdate->getContextSizeHMM();
		m_iContextSizeCW = m_hmmManagerUpdate->getContextSizeHMMCW();
		m_iPhoneContextPadding = m_phoneSet->size();
	}
	assert(m_iContextSizeWW >= m_iContextSizeCW);
	int iNodesPhoneFinal = iNodesPhone;
	if ((m_iContextSizeWW > 0) || (m_iContextSizeCW > 0)) {
		propagateContext(nodePhoneInitial,iNodesPhone,&iNodesPhoneFinal);
	}
	//print(nodePhoneInitial);
	// (4) create the HMM-graph from the phone-graph
	int iNodesHMM = -1;
	*nodeHMMInitial = create(nodePhoneInitial,&iNodesHMM);
	destroy(nodePhoneInitial,iNodesPhoneFinal);
	//print(*nodeHMMInitial);	
	
	// (5) compact the HMM-graph to reduce its size and renumber the remaining nodes
	int iNodesHMMCompacted = -1;
	int iEdgesHMMCompacted = -1;
	FBNodeHMM **nodesHMM = compact(*nodeHMMInitial,iNodesHMM,&iNodesHMMCompacted,&iEdgesHMMCompacted,nodeHMMFinal);
	*iNodes = iNodesHMMCompacted;
	*iEdges = iEdgesHMMCompacted;
	//compact(nodeHMMInitial,iNodesHMM);
	//print(*nodeHMMInitial);
	
	// (6) compute the distance from each edge to the initial and final edges in the graph
	computeDistances(nodesHMM,*nodeHMMInitial,*nodeHMMFinal,iNodesHMMCompacted,iEdgesHMMCompacted);	
	
	//double dTimeEnd = TimeUtils::getTimeMilliseconds();
	//double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	//printf("time: %.8f\n",dTimeSeconds);

	return nodesHMM;
}

// create a sausage of lexical units (it uses all the alternative pronunciations)
/*FBNodeLexUnit *HMMGraph::create(VLexUnitX &vLexUnitTranscription, VLexUnit &vLexUnitOptional, int *iNodes) {

	int iNode = 0;
	FBNodeLexUnit *nodeInitial = newFBNodeLexUnit(iNode++);
	FBNodeLexUnit *nodePrev = nodeInitial;	
	FBNodeLexUnit *nodeNext = NULL;	
	
	int iBoundary = 0;
	for(VLexUnitX::iterator it = vLexUnitTranscription.begin() ; it != vLexUnitTranscription.end() ; ++it, ++iBoundary) {	
		if ((iBoundary > 1) && (iBoundary < vLexUnitTranscription.size()-1) && (vLexUnitOptional.empty() == false)) {	
			nodeNext = newFBNodeLexUnit(iNode++);
			for(VLexUnit::iterator jt = vLexUnitOptional.begin() ; jt != vLexUnitOptional.end() ; ++jt) {
				newFBEdgeLexUnit(nodePrev,nodeNext,*jt);
			}
			nodePrev = nodeNext;
		}
		nodeNext = newFBNodeLexUnit(iNode++);
		for(VLexUnit::iterator jt = (*it)->vLexUnitPronunciations.begin() ; jt != (*it)->vLexUnitPronunciations.end() ; ++jt) {	
			newFBEdgeLexUnit(nodePrev,nodeNext,*jt);
		}
		nodePrev = nodeNext;
	}
	
	*iNodes = iNode;
	
	return nodeInitial;
}*/

// create a sausage of lexical units (it uses only the given pronunciation)
// - optional symbols are inserted at word-boundaries only if they are not redundant (i.e. if they
//   do not appear already on the transcription)
FBNodeLexUnit *HMMGraph::create(VLexUnit &vLexUnitTranscription, VLexUnit &vLexUnitOptional, int *iNodes) {

	if (vLexUnitTranscription.empty() == true) {
		return NULL;
	}

	int iNode = 0;
	FBNodeLexUnit *nodeInitial = newFBNodeLexUnit(iNode++);
	FBNodeLexUnit *nodePrev = nodeInitial;	
	FBNodeLexUnit *nodeNext = NULL;	
	
	int iBoundary = 0;
	LexUnit *lexUnitPrev = NULL;
	for(VLexUnit::iterator it = vLexUnitTranscription.begin() ; it != vLexUnitTranscription.end() ; ++it, ++iBoundary) {
		nodeNext = newFBNodeLexUnit(iNode++);
		if (vLexUnitOptional.empty() == false) {
			FBNodeLexUnit *nodeMiddle = NULL;
			FBNodeLexUnit *nodeMiddle2 = NULL;
			for(VLexUnit::iterator jt = vLexUnitOptional.begin() ; jt != vLexUnitOptional.end() ; ++jt) {
				if ((*jt != lexUnitPrev) && (*jt != *it)) {
					if (nodeMiddle == NULL) {
						nodeMiddle = newFBNodeLexUnit(iNode++);
					}
					newFBEdgeLexUnit(nodePrev,nodeMiddle,*jt);
				}
				VLexUnit::iterator kt = it;
				// final lexical unit: add trailing optional symbol
				if ((++kt == vLexUnitTranscription.end()) && (*jt != *it)) {
					assert(nodeMiddle != NULL);
					if (nodeMiddle2 == NULL) {
						nodeMiddle2 = newFBNodeLexUnit(iNode++);
					}
					newFBEdgeLexUnit(nodeMiddle2,nodeNext,*jt);
				}
			}			
			if (nodeMiddle != NULL) {
				newFBEdgeLexUnit(nodeMiddle,nodeNext,*it);
			}
			if (nodeMiddle2 != NULL) {
				newFBEdgeLexUnit(nodeMiddle,nodeMiddle2,*it);
			}
		}
		newFBEdgeLexUnit(nodePrev,nodeNext,*it);
		nodePrev = nodeNext;
		lexUnitPrev = *it;
	}
	
	*iNodes = iNode;
	
	return nodeInitial;
}

// create a sausage of lexical units (it uses only the given pronunciation)
// - optional symbols are inserted at word-boundaries only if they are not redundant (i.e. if they
//   do not appear already on the transcription)
FBNodeLexUnit *HMMGraph::create(VLexUnitX &vLexUnitTranscription, VLexUnit &vLexUnitOptional, int *iNodes) {

	if (vLexUnitTranscription.empty() == true) {
		return NULL;
	}

	int iNode = 0;
	FBNodeLexUnit *nodeInitial = newFBNodeLexUnit(iNode++);
	FBNodeLexUnit *nodePrev = nodeInitial;	
	FBNodeLexUnit *nodeNext = NULL;	
	
	int iBoundary = 0;
	LexUnitX *lexUnitPrev = NULL;
	for(VLexUnitX::iterator it = vLexUnitTranscription.begin() ; it != vLexUnitTranscription.end() ; ++it, ++iBoundary) {
		nodeNext = newFBNodeLexUnit(iNode++);
		if (vLexUnitOptional.empty() == false) {
			FBNodeLexUnit *nodeMiddle = NULL;
			FBNodeLexUnit *nodeMiddle2 = NULL;			
			for(VLexUnit::iterator jt = vLexUnitOptional.begin() ; jt != vLexUnitOptional.end() ; ++jt) {
				// check whether the optional symbol is the same than the lexical unit
				bool bSame = false;
				for(VLexUnit::iterator lt = (*it)->vLexUnitPronunciations.begin() ; lt != (*it)->vLexUnitPronunciations.end() ; ++lt) {
					if (*lt == *jt) {
						bSame = true;
						break;
					}
				}	
				// check whether the optional symbol is the same than the previous lexical unit
				bool bSamePrev = false;
				if (lexUnitPrev != NULL) {
					for(VLexUnit::iterator lt = lexUnitPrev->vLexUnitPronunciations.begin() ; lt != lexUnitPrev->vLexUnitPronunciations.end() ; ++lt) {
						if (*lt == *jt) {
							bSamePrev = true;
							break;
						}
					}	
				}
				// create the destination state for optional symbols?
				if ((bSamePrev == false) && (bSame == false)) {
					if (nodeMiddle == NULL) {
						nodeMiddle = newFBNodeLexUnit(iNode++);
					}
					newFBEdgeLexUnit(nodePrev,nodeMiddle,*jt);
				}
				VLexUnitX::iterator kt = it;
				// final lexical unit: add trailing optional symbol
				if ((++kt == vLexUnitTranscription.end()) && (bSame == false)) {
					assert(nodeMiddle != NULL);
					if (nodeMiddle2 == NULL) {
						nodeMiddle2 = newFBNodeLexUnit(iNode++);
					}
					newFBEdgeLexUnit(nodeMiddle2,nodeNext,*jt);
				}
			}			
			if (nodeMiddle != NULL) {
				for(VLexUnit::iterator lt = (*it)->vLexUnitPronunciations.begin() ; lt != (*it)->vLexUnitPronunciations.end() ; ++lt) {
					newFBEdgeLexUnit(nodeMiddle,nodeNext,*lt);
				}
			}
			if (nodeMiddle2 != NULL) {
				for(VLexUnit::iterator lt = (*it)->vLexUnitPronunciations.begin() ; lt != (*it)->vLexUnitPronunciations.end() ; ++lt) {
					newFBEdgeLexUnit(nodeMiddle,nodeMiddle2,*lt);
				}
			}
		}
		for(VLexUnit::iterator lt = (*it)->vLexUnitPronunciations.begin() ; lt != (*it)->vLexUnitPronunciations.end() ; ++lt) {
			newFBEdgeLexUnit(nodePrev,nodeNext,*lt);
		}
		nodePrev = nodeNext;
		lexUnitPrev = *it;
	}
	
	*iNodes = iNode;
	
	return nodeInitial;
}

// create a sausage of phones from the sausage of lexical units
FBNodePhone *HMMGraph::create(FBNodeLexUnit *nodeInitial, int iNodesLexUnit, int *iNodes, FBNodePhone **nodePhoneFinal) {

	int iNode = 0;
	//FBNodeLexUnit *nodeNext = NULL;
	FBNodePhone *nodePhoneInitial = newFBNodePhone(iNode++);
	*nodePhoneFinal = NULL;
	//FBNodePhone *nodePhonePrev = nodePhoneInitial;
	//FBNodePhone *nodePhoneNext = NULL;
	unsigned char iPosition; 
	
	MFBNodeLexUnitNodePhone mNodePair;
	vector<pair<FBNodeLexUnit*,FBNodePhone*> > vNodePair;
	mNodePair.insert(MFBNodeLexUnitNodePhone::value_type(nodeInitial,nodePhoneInitial));
	vNodePair.push_back(pair<FBNodeLexUnit*,FBNodePhone*>(nodeInitial,nodePhoneInitial));
	
	while(vNodePair.empty() == false) {
	
		pair<FBNodeLexUnit*,FBNodePhone*> nodePair = vNodePair.back();
		vNodePair.pop_back();	
	
		for(FBEdgeLexUnit *edge = nodePair.first->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			
			FBNodePhone *nodePhoneNext = NULL;
			MFBNodeLexUnitNodePhone::iterator it = mNodePair.find(edge->nodeNext);
			if (it == mNodePair.end()) {
				nodePhoneNext = newFBNodePhone(iNode++);
				mNodePair.insert(MFBNodeLexUnitNodePhone::value_type(edge->nodeNext,nodePhoneNext));
				vNodePair.push_back(pair<FBNodeLexUnit*,FBNodePhone*>(edge->nodeNext,nodePhoneNext));	
				if (edge->nodeNext->edgeNext == NULL) {
					assert(*nodePhoneFinal == NULL);
					*nodePhoneFinal = nodePhoneNext;
				}
			} else {
				nodePhoneNext = it->second;
			}	
			
			FBNodePhone *nodePhonePrev = nodePair.second;
		
			// for each phoneme in the lexical unit
			int iPhoneIndex = 0;
			FBNodePhone *nodePhonePrevWithin = nodePhonePrev;
			for(vector<int>::iterator it = edge->lexUnit->vPhones.begin() ; it != edge->lexUnit->vPhones.end() ; ++it, ++iPhoneIndex) {
				// update the within-word position
				if (edge->lexUnit->vPhones.size() == 1) {
					iPosition = WITHIN_WORD_POSITION_MONOPHONE;
				} 
				else if (iPhoneIndex == 0) {
					iPosition = WITHIN_WORD_POSITION_START;
				}
				else {
					assert(iPhoneIndex > 0);
					if (iPhoneIndex+1 == (int)edge->lexUnit->vPhones.size()) {
						iPosition = WITHIN_WORD_POSITION_END;
					} else {
						iPosition = WITHIN_WORD_POSITION_INTERNAL;
					}
				}
				FBEdgePhone *edgePhone = new FBEdgePhone;
				edgePhone->lexUnit = edge->lexUnit;
				edgePhone->iContextLeft = NULL;
				edgePhone->iContextRight = NULL;
				edgePhone->iPhone = *it;
				edgePhone->iPosition = iPosition;
				edgePhone->nodePrev = nodePhonePrevWithin;
				edgePhone->edgePrev = nodePhonePrevWithin->edgeNext;
				//assert((edgePhone->edgePrev != NULL) || (edgePhone->edgeNext != NULL));
				nodePhonePrevWithin->edgeNext = edgePhone;
				if ((iPhoneIndex+1 == (int)edge->lexUnit->vPhones.size())) {
					if (edge->nodePrev != edge->nodeNext) {
						edgePhone->nodeNext = nodePhoneNext;
					} else {
						edgePhone->nodeNext = edgePhone->nodePrev;
					}
				} else {
					edgePhone->nodeNext = newFBNodePhone(iNode++);
				}
				edgePhone->edgeNext = edgePhone->nodeNext->edgePrev;
				edgePhone->nodeNext->edgePrev = edgePhone;
				nodePhonePrevWithin = edgePhone->nodeNext;
			}
		}
	}	
	
	assert(*nodePhoneFinal != NULL);
	
	*iNodes = iNode;

	return nodePhoneInitial;
}


// propagate context across all phones (expansion might be needed)
// - all edges coming from a node have the same phone attached
// - all edges going to a node have the same phone attached
// context-dependency issues:
// - cross-word context is limited to adjacent words
// - only first and last phone within the word are cross-word context dependent, thus internal
//   phones are across-word context dependent but cross-word context independent
void HMMGraph::propagateContext(FBNodePhone *nodeInitial, int iNodesOriginal, int *iNodesFinal) {	

	assert(m_iContextSizeCW > 0);
	assert(m_iContextSizeWW > 0);
	
	int iContextSizeMax = max(m_iContextSizeCW,m_iContextSizeWW);

	FBNodePhone *nodeFinal = NULL;
	int iNode = iNodesOriginal;
	int iEdgesAdded = 0;
	int iNodesAdded = 0;	
	
	// (1) left context

	LFBEdgePhone lEdge;
	
 	// create a fake edge to make the algorithm simpler
	FBEdgePhone *edgeFake = new FBEdgePhone;
	edgeFake->nodeNext = nodeInitial;
	if (iContextSizeMax > 1) {
		edgeFake->iContextLeft = new unsigned char[iContextSizeMax];
		for(int i=0 ; i < iContextSizeMax ; ++i) {
			edgeFake->iContextLeft[i] = m_phoneSet->getPhoneIndexSilence();
		}			
	}
	edgeFake->lexUnit = m_lexiconManager->getLexUnitSilence();
	edgeFake->iPhone = m_phoneSet->getPhoneIndexSilence();
	edgeFake->iPosition = WITHIN_WORD_POSITION_MONOPHONE;
	lEdge.push_back(edgeFake);
	
	// process the edges one by one
	while(lEdge.empty() == false) {
		FBEdgePhone *edgeAux = lEdge.front();
		lEdge.pop_front();
		
		// does the edge go to the final node?
		if (edgeAux->nodeNext->edgeNext == NULL) {
			// keep the final node
			if (nodeFinal == NULL) {
				nodeFinal = edgeAux->nodeNext;
			} else {
				assert(nodeFinal == edgeAux->nodeNext);
			}
			continue;	
		}
		
		// successor edges do not have a context yet
		bool bContextLeft = false;
		for(FBEdgePhone *edge = edgeAux->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (edge == edgeAux) {
				continue;
			}
			if (edge->iContextLeft != NULL) {
				bContextLeft = true;
				break;
			}
		}	
		if (bContextLeft == false) {
			// propagate context and put the edges in the queue
			for(FBEdgePhone *edge = edgeAux->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
				edge->iContextLeft = newBlankContext();	
				shiftLeftContext(edgeAux,edge->iContextLeft);
				lEdge.push_back(edge);
			}
		} 
		// successor edges already have a context
		else {
			unsigned char *iContextLeftAux = newBlankContext();
			shiftLeftContext(edgeAux,iContextLeftAux);
			// different context: duplicate node and edges and put the edges in the queue
			if (memcmp(iContextLeftAux,edgeAux->nodeNext->edgeNext->iContextLeft,iContextSizeMax*sizeof(unsigned char)) != 0) {	
				// extract the auxiliar edge from the original destination node
				FBEdgePhone **edge = &(edgeAux->nodeNext->edgePrev);
				while(*edge != edgeAux) {	
					edge = &((*edge)->edgeNext);
				}
				*edge = (*edge)->edgeNext;
				// create a new destination node and copy edges
				FBNodePhone *node = newFBNodePhone(iNode++);
				++iNodesAdded;
				for(FBEdgePhone *edge = edgeAux->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
					if (edge == edgeAux) {
						continue;
					}
					FBEdgePhone *edgeDup = newFBEdgePhone(node,edge->nodeNext);	
					edgeDup->lexUnit = edge->lexUnit;
					edgeDup->iContextLeft = newBlankContext();
					edgeDup->iContextRight = NULL;
					edgeDup->iPhone = edge->iPhone;
					edgeDup->iPosition = edge->iPosition;
					memcpy(edgeDup->iContextLeft,iContextLeftAux,iContextSizeMax*sizeof(unsigned char));
					assert((edgeDup->edgePrev != NULL) || (edgeDup->edgeNext != NULL));
					lEdge.push_back(edgeDup);
					++iEdgesAdded;
				}
				edgeAux->nodeNext = node;
				node->edgePrev = edgeAux;
				edgeAux->edgeNext = NULL;
			} 
			// identical context: do nothing
			else {	
			}
			delete [] iContextLeftAux;
		}
	}
	if (iContextSizeMax > 1) {	
		delete [] edgeFake->iContextLeft;
	}
	delete edgeFake;
	
	//print(nodeInitial);
	
	// (2) right context
	
	assert(nodeFinal != NULL);
	assert(lEdge.empty() == true);
	
 	// create a fake edge to make the algorithm simpler
	edgeFake = new FBEdgePhone;
	edgeFake->nodePrev = nodeFinal;
	if (iContextSizeMax > 1) {
		edgeFake->iContextRight = new unsigned char[iContextSizeMax];
		for(int i=0 ; i < iContextSizeMax ; ++i) {
			edgeFake->iContextRight[i] = m_phoneSet->getPhoneIndexSilence();
		}
	}
	edgeFake->lexUnit = m_lexiconManager->getLexUnitSilence();	
	edgeFake->iPhone = m_phoneSet->getPhoneIndexSilence();
	edgeFake->iPosition = WITHIN_WORD_POSITION_MONOPHONE;	
	lEdge.push_back(edgeFake);
	
	// process the edges one by one
	while(lEdge.empty() == false) {
		FBEdgePhone *edgeAux = lEdge.front();
		lEdge.pop_front();
		
		// does the come from the inital node?
		if (edgeAux->nodePrev->edgePrev == NULL) {
			assert(edgeAux->nodePrev == nodeInitial);
			continue;	
		}
		
		// predecessor edges do not have a context yet	
		bool bContextRight = false;
		for(FBEdgePhone *edge = edgeAux->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
			if (edge == edgeAux) {
				continue;
			}
			if (edge->iContextRight != NULL) {
				bContextRight = true;
				break;
			}
		}	
		if (bContextRight == false) {		
			// propagate context and put the edges in the queue
			for(FBEdgePhone *edge = edgeAux->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
				edge->iContextRight = newBlankContext();
				shiftRightContext(edgeAux,edge->iContextRight);
				lEdge.push_back(edge);
			}
		} 
		// predecessor edges already have a context
		else {
			unsigned char *iContextRightAux = newBlankContext();
			shiftRightContext(edgeAux,iContextRightAux);
			// different context: duplicate node and edges and put the edges in the queue
			if (memcmp(iContextRightAux,edgeAux->nodePrev->edgePrev->iContextRight,iContextSizeMax*sizeof(unsigned char)) != 0) {	
				// extract the auxiliar edge from the original destination node
				FBEdgePhone **edge = &(edgeAux->nodePrev->edgeNext);
				while(*edge != edgeAux) {	
					edge = &((*edge)->edgePrev);
				}
				*edge = (*edge)->edgePrev;	
				// create a new source node and copy edges
				FBNodePhone *node = newFBNodePhone(iNode++);
				++iNodesAdded;
				for(FBEdgePhone *edge = edgeAux->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
					if (edge == edgeAux) {
						continue;
					}
					FBEdgePhone *edgeDup = newFBEdgePhone(edge->nodePrev,node);	
					edgeDup->lexUnit = edge->lexUnit;
					edgeDup->iContextRight = new unsigned char[iContextSizeMax];
					edgeDup->iContextLeft = new unsigned char[iContextSizeMax];
					edgeDup->iPhone = edge->iPhone;
					edgeDup->iPosition = edge->iPosition;
					memcpy(edgeDup->iContextRight,iContextRightAux,iContextSizeMax*sizeof(unsigned char));
					memcpy(edgeDup->iContextLeft,edge->iContextLeft,iContextSizeMax*sizeof(unsigned char));
					lEdge.push_back(edgeDup);
					++iEdgesAdded;
				}
				edgeAux->nodePrev = node;
				edgeAux->edgePrev = NULL;
				node->edgeNext = edgeAux;
			} 
			// identical context: do nothing
			else {	
			}
			
			delete [] iContextRightAux;
		}
	}
	if (iContextSizeMax > 1) {
		delete [] edgeFake->iContextRight;
	}
	delete edgeFake;
	
	//print(nodeInitial);
	//renumberNodes(nodeInitial);
	
	// sanity check
	check(nodeInitial);
	//print(nodeInitial);
	
	/*printf("# nodes added: %d\n",iNodesAdded);
	printf("# edges added: %d\n",iEdgesAdded);*/
	
	*iNodesFinal = iNode;
}

// check that for each node all the successors and predecessors have the same phone attached
void HMMGraph::check(FBNodePhone *nodeInitial) {

	bool *bPhones = new bool[m_phoneSet->size()];

	MFBNodePhone mNode;
	VFBNodePhone vNode;
	vNode.push_back(nodeInitial);
	mNode.insert(MFBNodePhone::value_type(nodeInitial,true));
	
	while(vNode.empty() == false) {
		
		FBNodePhone *nodeAux = vNode.back();
		vNode.pop_back();	
		
		// check successor edges
		unsigned char iPhone = UCHAR_MAX;
		for(FBEdgePhone *edge = nodeAux->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (iPhone == UCHAR_MAX) {
				iPhone = edge->iPhone;
			} else {
				assert((edge->iPhone == iPhone) || (nodeAux == nodeInitial));
			}
		}
		// check predecessor edges
		iPhone = UCHAR_MAX;
		for(FBEdgePhone *edge = nodeAux->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
			if (iPhone == UCHAR_MAX) {
				iPhone = edge->iPhone;
			} else {
				assert((edge->iPhone == iPhone) || (nodeAux->edgeNext == NULL));
			}	
		}
		
		for(FBEdgePhone *edge = nodeAux->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (mNode.find(edge->nodeNext) == mNode.end()) {
				vNode.push_back(edge->nodeNext);
				mNode.insert(MFBNodePhone::value_type(edge->nodeNext,true));
			}
		}
	}
	
	delete [] bPhones;

	//printf("check completed: %d nodes\n",mNode.size());
}


// enumerate the nodes
int HMMGraph::renumberNodes(FBNodePhone *nodeInitial) {

	int iNode = 0;
	MFBNodePhone mNode;
	VFBNodePhone vNode;
	vNode.push_back(nodeInitial);
	mNode.insert(MFBNodePhone::value_type(nodeInitial,true));
	
	while(vNode.empty() == false) {
		
		FBNodePhone *nodeAux = vNode.back();
		vNode.pop_back();
		nodeAux->iNode = iNode++;
		
		// check redundancy
		for(FBEdgePhone *edge = nodeAux->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		
			if (mNode.find(edge->nodeNext) == mNode.end()) {
				vNode.push_back(edge->nodeNext);
				mNode.insert(MFBNodePhone::value_type(edge->nodeNext,true));
			}
		}
	}

	return iNode;
}

// create a HMM-graph from a phone-sausage 
// - expansion might be needed because of the phonetic context
// - for global statistic accumulation it is necessary to keep the accumulator (it will
//   be necessary for the edge-merging state)
FBNodeHMM *HMMGraph::create(FBNodePhone *nodePhoneInitial, int *iNodesHMM) {

	int iNode = 0;
	MFBNodePhoneNodeHMM mNode;
	vector<pair<FBNodePhone*,FBNodeHMM*> > vNodePair;		// keeps pairs of phone-hmm nodes
	FBNodeHMM *nodeHMMInitial = newFBNodeHMM(iNode++);
	vNodePair.push_back(pair<FBNodePhone*,FBNodeHMM*>(nodePhoneInitial,nodeHMMInitial));
	mNode.insert(MFBNodePhoneNodeHMM::value_type(nodePhoneInitial,nodeHMMInitial));	
	
	while(vNodePair.empty() == false) {
		
		pair<FBNodePhone*,FBNodeHMM*> nodePair = vNodePair.back();
		vNodePair.pop_back();
		
		// convert the phone-edges to a sequence of HMM-nodes
		for(FBEdgePhone *edge = nodePair.first->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		
			FBNodeHMM *nodePrev = nodePair.second;
			
			// create an edge for each HMM-state/accumulator
			for(int iState=0 ; iState < NUMBER_HMM_STATES ;++iState) {
				VAccumulator vAccumulator;
				HMMState *hmmStateEstimation = NULL;
				HMMState *hmmStateUpdate = NULL;	
				// HMM-state (for both logical and physical n-phones)				
				// mode: hmm-parameter estimation
				if (m_hmmManagerUpdate->getPurpose() == HMM_PURPOSE_ESTIMATION) {
					// logical accumulator
					if (m_hmmManagerUpdate->areAccumulatorsLogical()) {
						Accumulator *accumulator = m_hmmManagerUpdate->getAccumulator(edge->iContextLeft,
							edge->iPhone,edge->iContextRight,edge->iPosition,iState); 
						vAccumulator.push_back(accumulator);
						assert(accumulator != NULL);
					}
					// physical accumulator 
					else {	
						hmmStateUpdate = m_hmmManagerUpdate->getHMMState(edge->iContextLeft,edge->iPhone,edge->iContextRight,
							edge->iPosition,iState);
						m_hmmManagerUpdate->getAccumulators(hmmStateUpdate->getId(),&vAccumulator);
					}			
					hmmStateEstimation = m_hmmManagerEstimation->getHMMState(edge->iContextLeft,edge->iPhone,edge->iContextRight,
						edge->iPosition,iState);
				} 
				// mode: hmm evaluation
				else {					
					assert(m_hmmManagerUpdate->getPurpose() == HMM_PURPOSE_EVALUATION);
					assert(m_hmmManagerEstimation == m_hmmManagerUpdate);
					hmmStateEstimation = (HMMState*)m_hmmManagerEstimation->getHMMStateDecoding(edge->iContextLeft,edge->iPhone,
						edge->iContextRight,edge->iPosition,iState);
					hmmStateUpdate = hmmStateEstimation;
				}
				if (iState == NUMBER_HMM_STATES-1) {
					MFBNodePhoneNodeHMM::iterator it = mNode.find(edge->nodeNext);
					if (it == mNode.end()) {
						FBEdgeHMM *edgeHMM = newFBEdgeHMM(nodePrev,newFBNodeHMM(iNode++),
							hmmStateEstimation,hmmStateUpdate,vAccumulator,edge->lexUnit,edge->iPosition);
						vNodePair.push_back(pair<FBNodePhone*,FBNodeHMM*>(edge->nodeNext,edgeHMM->nodeNext));
						mNode.insert(MFBNodePhoneNodeHMM::value_type(edge->nodeNext,edgeHMM->nodeNext));
					} else {
						newFBEdgeHMM(nodePrev,it->second,
							hmmStateEstimation,hmmStateUpdate,vAccumulator,edge->lexUnit,edge->iPosition);
					}
				} else {
					FBEdgeHMM *edgeHMM = newFBEdgeHMM(nodePrev,newFBNodeHMM(iNode++),
						hmmStateEstimation,hmmStateUpdate,vAccumulator,edge->lexUnit,edge->iPosition);
					nodePrev = edgeHMM->nodeNext;
				}
			}	
		}
	}
	
	*iNodesHMM = iNode;

	return nodeHMMInitial;
}



// print the graph of phones
void HMMGraph::print(FBNodePhone *nodeInitial) {

	MFBNodePhone mNode;

	VFBNodePhone vNode;
	vNode.push_back(nodeInitial);
	mNode.insert(MFBNodePhone::value_type(nodeInitial,true));
	
	while(vNode.empty() == false) {
		
		FBNodePhone *nodeAux = vNode.back();
		vNode.pop_back();
		
		for(FBEdgePhone *edge = nodeAux->edgeNext ; edge != NULL ; edge = edge->edgePrev) {		
			printf("(%d -> %d) [%10s] ",edge->nodePrev->iNode,edge->nodeNext->iNode,m_phoneSet->getStrPhone(edge->iPhone));	
			if (edge->iContextLeft != NULL) {
				for(int i=0 ; i < max(m_iContextSizeWW,m_iContextSizeCW) ; ++i) {
					if (edge->iContextLeft[i] == m_iPhoneContextPadding) {
						printf("l= %s ","<pad>");
					} else {
						printf("l= %s ",m_phoneSet->getStrPhone(edge->iContextLeft[i]));
					}
				}
				printf("| ");
			} 
			if (edge->iContextRight != NULL) {
				for(int i=0 ; i < max(m_iContextSizeWW,m_iContextSizeCW) ; ++i) {
					if (edge->iContextRight[i] == m_iPhoneContextPadding) {
						printf("r= %s ","<pad>");
					} else {
						printf("r= %s ",m_phoneSet->getStrPhone(edge->iContextRight[i]));
					}
				}
			}
			printf("\n");
			if (mNode.find(edge->nodeNext) == mNode.end()) {
				vNode.push_back(edge->nodeNext);
				mNode.insert(MFBNodePhone::value_type(edge->nodeNext,true));
			}
		}
	}
	
	cout << "# nodes: " << mNode.size() << endl;
}

// print the graph of HMMs
void HMMGraph::print(FBNodeHMM *nodeInitial) {

	int iEdges = 0;
	MFBNodeHMM mNode;
	VFBNodeHMM vNode;
	vNode.push_back(nodeInitial);
	mNode.insert(MFBNodeHMM::value_type(nodeInitial,true));
	
	printf("-----------------------------------------\n");
	while(vNode.empty() == false) {
		
		FBNodeHMM *nodeAux = vNode.back();
		vNode.pop_back();
		
		for(FBEdgeHMM *edge = nodeAux->edgeNext ; edge != NULL ; edge = edge->edgePrev) {		
			++iEdges;
			
			char strLexUnitPronunciation[MAX_LEXUNIT_LENGTH+1];
			m_lexiconManager->getStrLexUnitPronunciation(edge->lexUnit,strLexUnitPronunciation);
			
			if (m_hmmManagerUpdate->getPurpose() == HMM_PURPOSE_ESTIMATION) {
				printf("%5d %5d %5d %20s %d\n",edge->nodePrev->iNode,edge->hmmStateEstimation->getId(),
					edge->nodeNext->iNode,strLexUnitPronunciation,edge->iPosition);
			} else {
				HMMStateDecoding *hmmStateDecoding = (HMMStateDecoding*)edge->hmmStateEstimation;
				printf("%5d %5d %5d %20s %d\n",edge->nodePrev->iNode,hmmStateDecoding->getId(),
					edge->nodeNext->iNode,strLexUnitPronunciation,edge->iPosition);
			}
				
			if (mNode.find(edge->nodeNext) == mNode.end()) {
				vNode.push_back(edge->nodeNext);
				mNode.insert(MFBNodeHMM::value_type(edge->nodeNext,true));
			}
			assert(edge->nodePrev == nodeAux);	
		}
	}
	
	cout << "# nodes: " << mNode.size() << endl;
	cout << "# edges: " << iEdges << endl;
	cout << "-----------------------------------------" << endl;
}

// print the sausage of lexical units
void HMMGraph::print(FBNodeLexUnit *nodeInitial) {

	MFBNodeLexUnit mNode;
	mNode.insert(MFBNodeLexUnit::value_type(nodeInitial,true));
	LFBNodeLexUnit lNode;
	lNode.push_back(nodeInitial);
	while(lNode.empty() == false) {
	
		FBNodeLexUnit *node = lNode.front();
		lNode.pop_front();				
		
		for(FBEdgeLexUnit *edge = node->edgeNext ; edge != NULL ; edge = edge->edgePrev) {		
		
			char strLexUnitPronunciation[MAX_LEXUNIT_LENGTH+1];
			m_lexiconManager->getStrLexUnitPronunciation(edge->lexUnit,strLexUnitPronunciation);	
			cout << "(" << edge->nodePrev << " -> " << edge->nodeNext << ") " << setw(20) 
				<< strLexUnitPronunciation << " " << setw(2) << edge->lexUnit->vPhones.size() << " phones" << endl;
			if (mNode.find(edge->nodeNext) == mNode.end()) {
				lNode.push_back(edge->nodeNext);
				mNode.insert(MFBNodeLexUnit::value_type(edge->nodeNext,true));
			}			
		}		
	}
}


// destroy the lex-unit sausage
void HMMGraph::destroy(FBNodeLexUnit *nodeInitial) {

	FBNodeLexUnit *nodePrev = nodeInitial;		
	while(nodePrev->edgeNext != NULL) {
				
		FBNodeLexUnit *nodeAux = nodePrev->edgeNext->nodeNext;
		for(FBEdgeLexUnit *edge = nodePrev->edgeNext ; edge != NULL ; ) {
			FBEdgeLexUnit *edgeToDelete = edge;
			edge = edge->edgePrev;
			delete edgeToDelete;
		}
		delete nodePrev;
		nodePrev = nodeAux;
	}
	delete nodePrev;
}

// destroy the lex-unit graph
void HMMGraph::destroy(FBNodeLexUnit *nodeInitial, int iNodes) {

	FBNodeLexUnit **nodeLexUnit = new FBNodeLexUnit*[iNodes];
	for(int i=0 ; i < iNodes ; ++i) {
		nodeLexUnit[i] = NULL;
	}
	
	VFBNodeLexUnit vNodes;
	vNodes.push_back(nodeInitial);
	while(vNodes.empty() == false) {
		
		FBNodeLexUnit *node = vNodes.back();
		vNodes.pop_back();
		
		nodeLexUnit[node->iNode] = node;
		
		for(FBEdgeLexUnit *edge = node->edgeNext ; edge != NULL ; ) {
			if (nodeLexUnit[edge->nodeNext->iNode] == NULL) {
				nodeLexUnit[edge->nodeNext->iNode] = edge->nodeNext;
				vNodes.push_back(edge->nodeNext);
			}
			FBEdgeLexUnit *edgeToDelete = edge;
			edge = edge->edgePrev;
			delete edgeToDelete;
		}	
	}
	
	for(int i=0 ; i < iNodes ; ++i) {
		if (nodeLexUnit[i] != NULL) {
			delete nodeLexUnit[i];
		}
	}
	
	delete [] nodeLexUnit;
}

// compact the HMM-graph by applying forward/backward node merging
FBNodeHMM **HMMGraph::compact(FBNodeHMM *nodeInitial, int iNodes, int *iNodesAfterCompacting, 
	int *iEdgesAfterCompacting, FBNodeHMM **nodeHMMFinal) {

	LFBNodeHMM lNode;
	int iNodesDeleted = 0;
	int iEdgesDeleted = 0;
	FBNodeHMM *nodeFinal = NULL;
	
	// nodes processed: this array will be used to 
	FBNodeHMM **nodesProcessed = new FBNodeHMM*[iNodes];
	int iNodesProcessed = 0;	
	
	// keep the state of each node
	unsigned char *iNodeState = new unsigned char[iNodes];
	// set all the nodes as "not visited"
	for(int i=0 ; i < iNodes ; ++i) {
		iNodeState[i] = NODE_STATE_NOT_VISITED;
	}
	// nodes to be deleted at the end
	FBNodeHMM **nodesToDelete = new FBNodeHMM*[iNodes];
	int iNodesToDelete = 0;
	
	// forward compacting
	lNode.push_back(nodeInitial);
	iNodeState[nodeInitial->iNode] = NODE_STATE_QUEUED;
	
	while(lNode.empty() == false) {
		
		FBNodeHMM *nodeAux = lNode.front();
		lNode.pop_front();

		if (iNodeState[nodeAux->iNode] == NODE_STATE_REMOVED) {
			continue;
		}
		
		// check that all predecessor edges have been processed
		bool bReady = true;
		for(FBEdgeHMM *edge = nodeAux->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
			if (iNodeState[edge->nodePrev->iNode] != NODE_STATE_PROCESSED) {
				bReady = false;
				break;
			}
		}
		if (bReady == false) {
			lNode.push_back(nodeAux);
			continue;
		}
		
		// keep the final node
		if (nodeAux->edgeNext == NULL) {
			assert(nodeFinal == NULL);
			nodeFinal = nodeAux;
		}
				
		// look for edges to merge
		for(FBEdgeHMM *edge1 = nodeAux->edgeNext ; edge1 != NULL ; ) {
			bool bMerged = false;
			for(FBEdgeHMM *edge2 = nodeAux->edgeNext ; edge2 != edge1 ; edge2 = edge2->edgePrev) {
				//if (edge1->hmmState == edge2->hmmState) {	
				if (equal(edge1,edge2) == true) {	
					// check that destination nodes have the same predecessor edges
					if (samePredecessors(edge1->nodeNext,edge2->nodeNext) == true) {
						bMerged = true;
						// non-shared destination node: remove edge and node and transfer successor edges
						if (edge1->nodeNext != edge2->nodeNext) {
							//printf("(f) merging happens!! %d -> %d and %d -> %d\n",edge1->nodePrev->iNode,edge1->nodeNext->iNode,edge2->nodePrev->iNode,edge2->nodeNext->iNode);	
							iNodeState[edge1->nodeNext->iNode] = NODE_STATE_REMOVED;
							// merge the edges and nodes	
							// - move successors from the redundant node to the original
							FBEdgeHMM **edgeTail = &(edge1->nodeNext->edgeNext);
							for(FBEdgeHMM *edge3 = edge1->nodeNext->edgeNext ; edge3 != NULL ; edge3 = edge3->edgePrev) {
								edge3->nodePrev = edge2->nodeNext;
								edgeTail = &(edge3->edgePrev);
							}
							assert(*edgeTail == NULL);
							*edgeTail = edge2->nodeNext->edgeNext;
							edge2->nodeNext->edgeNext = edge1->nodeNext->edgeNext;
							// remove incoming edges that go to the redundant node (there are equivalent edges that go
							// to the original node)
							nodesToDelete[iNodesToDelete++] = edge1->nodeNext;
							for(FBEdgeHMM *edge3 = edge1->nodeNext->edgePrev ; edge3 != NULL ; ) {
								// - remove the redundant edge from the list of successor edges
								bool bFound = true;
								for(FBEdgeHMM **edge4 = &(edge3->nodePrev->edgeNext) ; *edge4 != NULL ; edge4 = &((*edge4)->edgePrev)) {
									if (*edge4 == edge3) {
										FBEdgeHMM *edgeToDelete = edge3;
										if (*edge4 == edge1) {
											edge1 = (*edge4)->edgePrev;
										}
										*edge4 = (*edge4)->edgePrev;
										edge3 = edge3->edgeNext;
										delete edgeToDelete;
										bFound = true;
										break;
									}
								}
								assert(bFound == true);
								++iEdgesDeleted;
							}	
							break;
						}
						// destination node is shared: just remove the edge
						else {
							assert(iNodeState[edge1->nodeNext->iNode] == NODE_STATE_QUEUED);
							// disconnect the edge from the source node
							bool bFound = false;
							for(FBEdgeHMM **edge3 = &(edge1->nodePrev->edgeNext) ; *edge3 != NULL ; edge3 = &((*edge3)->edgePrev)) {
								if (*edge3 == edge1) {
									*edge3 = (*edge3)->edgePrev;
									bFound = true;
									break;
								}
							}	
							assert(bFound == true);
							// disconnect the edge from the source node
							bFound = false;
							for(FBEdgeHMM **edge3 = &(edge1->nodeNext->edgePrev) ; *edge3 != NULL ; edge3 = &((*edge3)->edgeNext)) {
								if (*edge3 == edge1) {
									*edge3 = (*edge3)->edgeNext;
									bFound = true;
									break;
								}
							}	
							assert(bFound == true);		
							// go to the next outgoing edge of the node being processed
							FBEdgeHMM *edgeToDelete = edge1;	
							edge1 = edge1->edgePrev;
							delete edgeToDelete;
							++iEdgesDeleted;
							break;
						}
					}	
				} 
			}
			if (bMerged == false) {
				if (iNodeState[edge1->nodeNext->iNode] != NODE_STATE_QUEUED) {
					lNode.push_back(edge1->nodeNext);
					iNodeState[edge1->nodeNext->iNode] = NODE_STATE_QUEUED;
				}
				edge1 = edge1->edgePrev;	
			}
		}
		
		iNodeState[nodeAux->iNode] = NODE_STATE_PROCESSED;
	}
	
	assert(nodeFinal != NULL);
	assert(lNode.empty() == true);
	
	// set all the nodes as "not visited"
	for(int i=0 ; i < iNodes ; ++i) {
		iNodeState[i] = NODE_STATE_NOT_VISITED;
	}	
	
	// backward compacting
	lNode.push_back(nodeFinal);
	iNodeState[nodeFinal->iNode] = NODE_STATE_QUEUED;	
	
	while(lNode.empty() == false) {
		
		FBNodeHMM *nodeAux = lNode.front();
		lNode.pop_front();
		
		if (iNodeState[nodeAux->iNode] == NODE_STATE_REMOVED) {
			continue;
		}
		
		// check that all successor edges have been processed
		bool bReady = true;
		for(FBEdgeHMM *edge = nodeAux->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (iNodeState[edge->nodeNext->iNode] != NODE_STATE_PROCESSED) {
				bReady = false;
				break;
			}
		}
		if (bReady == false) {
			lNode.push_back(nodeAux);
			continue;
		}
		
		// look for edges to merge
		for(FBEdgeHMM *edge1 = nodeAux->edgePrev ; edge1 != NULL ; ) {
			bool bMerged = false;
			for(FBEdgeHMM *edge2 = nodeAux->edgePrev ; edge2 != edge1 ; edge2 = edge2->edgeNext) {
				//if (edge1->hmmState == edge2->hmmState) {	
				if (equal(edge1,edge2) == true) {	
					// check that destination nodes have the same successor edges
					if (sameSuccessors(edge1->nodePrev,edge2->nodePrev) == true) {
						assert(edge1->nodePrev != edge2->nodePrev);
						bMerged = true;
						//printf("(b) merging happens!!\n");
						iNodeState[edge1->nodePrev->iNode] = NODE_STATE_REMOVED;
						// merge the edges and nodes	
						// - move successors from the redundant node to the original
						FBEdgeHMM **edgeTail = &(edge1->nodePrev->edgePrev);
						for(FBEdgeHMM *edge3 = edge1->nodePrev->edgePrev ; edge3 != NULL ; edge3 = edge3->edgeNext) {
							edge3->nodeNext = edge2->nodePrev;
							edgeTail = &(edge3->edgeNext);
						}
						assert(*edgeTail == NULL);
						*edgeTail = edge2->nodePrev->edgePrev;
						edge2->nodePrev->edgePrev = edge1->nodePrev->edgePrev;
						// remove successor edges that go to the redundant node (they also go
						// to the original node)
						nodesToDelete[iNodesToDelete++] = edge1->nodePrev;
						for(FBEdgeHMM *edge3 = edge1->nodePrev->edgeNext ; edge3 != NULL ; ) {
							// - remove the redundant edge from the list of predecessor edges
							bool bFound = true;
							for(FBEdgeHMM **edge4 = &(edge3->nodeNext->edgePrev) ; *edge4 != NULL ; edge4 = &((*edge4)->edgeNext)) {
								if (*edge4 == edge3) {
									FBEdgeHMM *edgeToDelete = edge3;
									if (*edge4 == edge1) {
										edge1 = (*edge4)->edgeNext;
									}	
									*edge4 = (*edge4)->edgeNext;
									edge3 = edge3->edgePrev;
									delete edgeToDelete;
									bFound = true;
									break;
								}
							}
							assert(bFound == true);
							++iEdgesDeleted;
						}
						break;
					}	
				}
			}
			if (bMerged == false) {
				if (iNodeState[edge1->nodePrev->iNode] != NODE_STATE_QUEUED) {
					lNode.push_back(edge1->nodePrev);
					iNodeState[edge1->nodePrev->iNode] = NODE_STATE_QUEUED;
				}
				edge1 = edge1->edgeNext;	
			}
		}	
		
		iNodeState[nodeAux->iNode] = NODE_STATE_PROCESSED;
		nodesProcessed[iNodesProcessed++] = nodeAux;
	}
	
	// clean-up (note that nodes need to be removed at the end since 
	// when the node-merging is done they might be queued)
	for(int i=0 ; i < iNodesToDelete ; ++i) {
		delete nodesToDelete[i];
		++iNodesDeleted;
	}
	delete [] nodesToDelete;
	delete [] iNodeState;
	
	assert(iNodes-iNodesDeleted == iNodesProcessed);
	*iNodesAfterCompacting = iNodes-iNodesDeleted;
	
	FBNodeHMM **nodesCompacted = new FBNodeHMM*[iNodesProcessed];
	*iEdgesAfterCompacting = 0;
	int j=0;
	for(int i=0 ; i < iNodesProcessed ; ++i) {
		nodesCompacted[i] = nodesProcessed[i];
		nodesCompacted[i]->iNode = i;
		// count the edges
		for(FBEdgeHMM *edge = nodesCompacted[i]->edgeNext ; edge != NULL ; edge = edge->edgePrev, ++j) {
			++(*iEdgesAfterCompacting);
			edge->iEdge = j;
		}
	}	
	
	*nodeHMMFinal = nodeFinal;
	
	delete [] nodesProcessed;
	
	/*printf("-- HMM compacting ------------------------------------\n");
	printf("# nodes: %d\n",iNodes-iNodesDeleted);
	//printf("# edges: %d\n",iEdges-iEdgesDeleted);
	printf("# nodes deleted: %d\n",iNodesDeleted);
	printf("# edges deleted: %d\n",iEdgesDeleted);
	printf("------------------------------------------------------\n");
	//printf("nodes processed: %d\n",mNodeProcessed.size());
	print(nodeInitial);*/

	return nodesCompacted;
}

// return whether the nodes have the same set of predecessors
bool HMMGraph::samePredecessors(FBNodeHMM *node1, FBNodeHMM *node2) {

	// same # of predecessors?
	int iPredecessors = 0;
	for(FBEdgeHMM *edge1 = node1->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
		++iPredecessors;
	}
	for(FBEdgeHMM *edge2 = node2->edgePrev ; edge2 != NULL ; edge2 = edge2->edgeNext) {
		--iPredecessors;
	}
	if (iPredecessors != 0) {
		return false;
	}
		
	// same predecessors?
	for(FBEdgeHMM *edge1 = node1->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
		bool bFound = false;
		for(FBEdgeHMM *edge2 = node2->edgePrev ; edge2 != NULL ; edge2 = edge2->edgeNext) {
			//if ((edge1->hmmState == edge2->hmmState) && (edge1->nodePrev == edge2->nodePrev)) {
			if ((equal(edge1,edge2) == true) && (edge1->nodePrev == edge2->nodePrev)) {
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

// return whether the nodes have the same set of successors
bool HMMGraph::sameSuccessors(FBNodeHMM *node1, FBNodeHMM *node2) {

	// same # of successors?
	int iSuccessors = 0;
	for(FBEdgeHMM *edge1 = node1->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
		++iSuccessors;
	}
	for(FBEdgeHMM *edge2 = node2->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
		--iSuccessors;
	}
	if (iSuccessors != 0) {
		return false;
	}
		
	// same successors?
	for(FBEdgeHMM *edge1 = node1->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
		bool bFound = false;
		for(FBEdgeHMM *edge2 = node2->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
			//if ((edge1->hmmState == edge2->hmmState) && (edge1->nodeNext == edge2->nodeNext)) {
			if ((equal(edge1,edge2) == true) && (edge1->nodeNext == edge2->nodeNext)) {
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

// destroy the HMM-graph
void HMMGraph::destroy(FBNodeHMM **nodes, int iNodes) {

	for(int i=0 ; i < iNodes ; ++i) {
		for(FBEdgeHMM *edge = nodes[i]->edgeNext ; edge != NULL ; ) {
			FBEdgeHMM *edgeToDelete = edge;
			edge = edge->edgePrev;
			delete edgeToDelete;
		}
		delete nodes[i];
	}
	
	delete [] nodes;
}

// destroy the HMM-graph
void HMMGraph::destroy(FBNodeHMM *nodeInitial, int iNodes) {

	FBNodeHMM **nodeHMM = new FBNodeHMM*[iNodes];
	for(int i=0 ; i < iNodes ; ++i) {
		nodeHMM[i] = NULL;
	}
	
	VFBNodeHMM vNodes;
	vNodes.push_back(nodeInitial);
	while(vNodes.empty() == false) {
		
		FBNodeHMM *node = vNodes.back();
		vNodes.pop_back();
		
		nodeHMM[node->iNode] = node;
		
		for(FBEdgeHMM *edge = node->edgeNext ; edge != NULL ; ) {
			if (nodeHMM[edge->nodeNext->iNode] == NULL) {
				nodeHMM[edge->nodeNext->iNode] = edge->nodeNext;
				vNodes.push_back(edge->nodeNext);
			}
			FBEdgeHMM *edgeToDelete = edge;
			edge = edge->edgePrev;
			delete edgeToDelete;
		}	
	}
	
	for(int i=0 ; i < iNodes ; ++i) {
		if (nodeHMM[i] != NULL) {
			delete nodeHMM[i];
		}
	}
	
	delete [] nodeHMM;
}

// destroy the phone-graph
void HMMGraph::destroy(FBNodePhone *nodeInitial, int iNodes) {

	FBNodePhone **nodes = new FBNodePhone*[iNodes];
	for(int i=0 ; i < iNodes ; ++i) {
		nodes[i] = NULL;
	}
	
	VFBNodePhone vNodes;
	vNodes.push_back(nodeInitial);
	nodes[nodeInitial->iNode] = nodeInitial;	
	while(vNodes.empty() == false) {
			
		FBNodePhone *node = vNodes.back();
		vNodes.pop_back();
		
		for(FBEdgePhone *edge = node->edgeNext ; edge != NULL ; ) {
			assert(edge->nodeNext->iNode < iNodes);
			if (nodes[edge->nodeNext->iNode] == NULL) {
				nodes[edge->nodeNext->iNode] = edge->nodeNext;
				vNodes.push_back(edge->nodeNext);
			}
			FBEdgePhone *edgeToDelete = edge;	
			edge = edge->edgePrev;
			if (edgeToDelete->iContextLeft != NULL) {
				delete [] edgeToDelete->iContextLeft;
			}
			if (edgeToDelete->iContextRight != NULL) {
				delete [] edgeToDelete->iContextRight;
			}
			delete edgeToDelete;
		}
	}
	
	for(int i=0 ; i < iNodes ; ++i) {
		assert(nodes[i] != NULL);
		delete nodes[i];
	}
	
	delete [] nodes;
}

// compute the distance from each edge to the initial and final edges in the graph
// - it applies Dijkstra
void HMMGraph::computeDistances(FBNodeHMM **nodes, FBNodeHMM *nodeInitial, 
	FBNodeHMM *nodeFinal, int iNodes, int iEdges) {
	
	// (1) compute the distance from each node to the final node
	
	// intialize the state of all nodes to "not visited"
	unsigned char *iNodeState = new unsigned char[iNodes];
	for(int i=0 ; i <iNodes ; ++i) {
		iNodeState[i] = NODE_STATE_NOT_VISITED;
	}
	
	// initialize all distances to inf except for the final node
	for(int i=0 ; i<iNodes ; ++i) {
		nodes[i]->iDistanceEnd = INT_MAX;
	}
	nodeFinal->iDistanceEnd = 0;	
	iNodeState[nodeFinal->iNode] = NODE_STATE_QUEUED;
	
	VFBNodeHMM vNodes;
	vNodes.push_back(nodeFinal);
	while(vNodes.empty() == false) {
		
		FBNodeHMM *node = vNodes.back();
		vNodes.pop_back();
		
		iNodeState[node->iNode] = NODE_STATE_PROCESSED;
			
		for(FBEdgeHMM *edge = node->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
			// relaxation
			if (edge->nodePrev->iDistanceEnd > (node->iDistanceEnd+1)) {
				edge->nodePrev->iDistanceEnd = node->iDistanceEnd+1;
				if (iNodeState[edge->nodePrev->iNode] != NODE_STATE_QUEUED) {
					vNodes.push_back(edge->nodePrev);
					iNodeState[edge->nodePrev->iNode] = NODE_STATE_QUEUED;
				}
			}
		}
	}
		
	// (1) compute the distance from each node to the initial node
	
	// intialize the state of all nodes to "not visited"
	for(int i=0 ; i <iNodes ; ++i) {
		iNodeState[i] = NODE_STATE_NOT_VISITED;
	}
	
	// initialize all distances to inf except for the initial node
	for(int i=0 ; i<iNodes ; ++i) {
		nodes[i]->iDistanceStart = INT_MAX;
	}
	nodeInitial->iDistanceStart = 0;	
	iNodeState[nodeInitial->iNode] = NODE_STATE_QUEUED;
	
	assert(vNodes.empty() == true);
	vNodes.push_back(nodeInitial);
	while(vNodes.empty() == false) {
		
		FBNodeHMM *node = vNodes.back();
		vNodes.pop_back();
		
		iNodeState[node->iNode] = NODE_STATE_PROCESSED;
			
		for(FBEdgeHMM *edge = node->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			// relaxation
			if (edge->nodeNext->iDistanceStart > (node->iDistanceStart+1)) {
				edge->nodeNext->iDistanceStart = node->iDistanceStart+1;
				if (iNodeState[edge->nodeNext->iNode] != NODE_STATE_QUEUED) {
					vNodes.push_back(edge->nodeNext);
					iNodeState[edge->nodeNext->iNode] = NODE_STATE_QUEUED;
				}
			}
		}
	}
	
	// sanity check
	for(int i=0 ; i<iNodes ; ++i) {
		assert(nodes[i]->iDistanceStart != INT_MAX);
		assert(nodes[i]->iDistanceEnd != INT_MAX);
	}
	
	delete [] iNodeState;
}

};	// end-of-namespace
