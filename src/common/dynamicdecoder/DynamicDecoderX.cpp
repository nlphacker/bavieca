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


#include "DynamicDecoderX.h"

#include "BestPath.h"
#include "LMFSM.h"
#include "HMMManager.h"
#include "LMManager.h"
#include "LMLookAhead.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
DynamicDecoderX::DynamicDecoderX(PhoneSet *phoneSet, HMMManager *hmmManager, 
			LexiconManager *lexiconManager, LMManager *lmManager, float fLMScalingFactor, 
			DynamicNetworkX *dynamicNetwork, int iMaxActiveNodes, 
			int iMaxActiveNodesWE, int iMaxActiveTokensNode, float fBeamWidthNodes, 
			float fBeamWidthNodesWE, float fBeamWidthTokensNode, bool bWordGraphGeneration, 
			int iMaxWordSequencesState)
{
	m_phoneSet = phoneSet;
	m_hmmManager = hmmManager;
	m_lexiconManager = lexiconManager;
	m_lmManager = lmManager;
	m_lmFSM = lmManager->getFSM();
	m_fLMScalingFactor = fLMScalingFactor;
	m_iNGram = m_lmFSM->getNGramOrder();
	m_dynamicNetwork = dynamicNetwork;
	
	// network properties
	m_arcs = m_dynamicNetwork->getArcs(&m_iArcs);
	m_nodes = m_dynamicNetwork->getNodes(&m_iNodes);
	
	m_iTimeCurrent = -1;
	
	// pruning parameters
	m_iMaxActiveNodes = iMaxActiveNodes;
	m_iMaxActiveNodesWE = iMaxActiveNodesWE;
	m_iMaxActiveTokensNode = iMaxActiveTokensNode;	
	m_fBeamWidthNodes = fBeamWidthNodes;	
	m_fBeamWidthNodesWE = fBeamWidthNodesWE;	
	m_fBeamWidthTokensNode = fBeamWidthTokensNode;
		
	// active nodes
	m_nodesActiveCurrent = NULL;
	m_nodesActiveNext = NULL;
	
	// history items
	m_iHistoryItems = 0;
	m_historyItems = NULL;
	m_iHistoryItemBegSentence = -1;
	m_iHistoryItemAvailable = -1;
	m_iTimeGarbageCollectionLast = -1;
	
	// word-graph generation
	m_bLatticeGeneration = bWordGraphGeneration;
	m_iMaxWordSequencesState = iMaxWordSequencesState;
	m_iWSHashBuckets = UINT_MAX;
	m_iWSHashEntries = UINT_MAX;	
	m_wshashEntries = NULL;
	m_iWSHashEntryCollisionAvailable = -1;
	
	// wrod-graph tokens
	m_iWGTokens = 0;
	m_wgTokens = NULL;
	m_iWGTokenAvailable = -1;
	m_iWordSequenceAux = NULL;
	
	m_iLexUnitPronUnknown = m_lexiconManager->m_lexUnitUnknown->iLexUnitPron;
	
	// mark it as uninitialized
	m_bInitialized = false;
}

// destructor
DynamicDecoderX::~DynamicDecoderX()
{
}

// initialization
void DynamicDecoderX::initialize() {
	
	// max tokens per active arc
	m_iTokensNodeMax = m_iMaxActiveTokensNode*2;

	// allocate memory for the tables of active arcs
	m_iNodesActiveCurrentMax = m_iMaxActiveNodes*5;
	m_iNodesActiveNextMax = m_iMaxActiveNodes*5;	
	m_nodesActiveCurrent = new DNode*[m_iNodesActiveCurrentMax];
	m_nodesActiveNext = new DNode*[m_iNodesActiveNextMax];
	
	// allocate memory for the tokens
	m_iTokensMax = m_iNodesActiveCurrentMax*10;
	m_tokensCurrent = new Token[m_iTokensMax];
	m_tokensNext = new Token[m_iTokensMax];
	m_iTokensNext = 0;
	
	// pre-allocate memory for the active tokens 
	m_iActiveTokenMax = m_iNodesActiveCurrentMax*m_iTokensNodeMax;
	m_activeTokenCurrent = new ActiveToken[m_iActiveTokenMax];
	m_activeTokenNext = new ActiveToken[m_iActiveTokenMax];
	m_iActiveTokenTables = 0;
	
	// history item management (there exists garbage collection for history items)
	m_iHistoryItems = 10000;
	m_historyItems = new HistoryItem[m_iHistoryItems];
	m_iHistoryItemsAuxBuffer = new int[m_iTokensNodeMax];
	m_iHistoryItemsAux = NULL;
	m_iHistoryItemsAuxSize = -1;
	
	// language model look-ahead
	m_lmLookAhead = new LMLookAhead(m_lexiconManager,m_lmManager,m_dynamicNetwork,this,m_iTokensMax);
	m_lmLookAhead->initialize();
		
	// word-graph generation
	if (m_bLatticeGeneration) {
		// allocate memory for the hash table
		m_iWSHashBuckets = 100000;
		m_iWSHashEntries = 2*m_iWSHashBuckets;
		m_wshashEntries = new WSHashEntry[m_iWSHashEntries];
		// set all the entries to "old"
		for(unsigned int i=0 ; i < m_iWSHashEntries ; ++i) {
			m_wshashEntries[i].iTime = -1;	
		}	
		m_iWSHashEntryCollisionAvailable = m_iWSHashBuckets;
		
		// word-graph token management (there exists garbage collection for wg-tokens)
		m_iWGTokens = 5000*m_iMaxWordSequencesState;
		m_wgTokens = new WGToken[m_iWGTokens];
		for(unsigned int i=0 ; i < m_iWGTokens-1 ; i += m_iMaxWordSequencesState) {
			m_wgTokens[i].iActive = -1;
			m_wgTokens[i].iPrev = i+m_iMaxWordSequencesState;
		}
		m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iActive = -1;
		m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iPrev = -1;
		m_iWGTokenAvailable = 0;
		m_iWordSequenceAux = new int[m_iTokensNodeMax];
	}
	
	m_bInitialized = true;	
}

// uninitialize
void DynamicDecoderX::uninitialize() {
	
	assert(m_bInitialized);
	
	delete [] m_nodesActiveCurrent;
	delete [] m_nodesActiveNext;
	delete [] m_tokensCurrent;
	delete [] m_tokensNext;
	delete [] m_activeTokenCurrent;
	delete [] m_activeTokenNext;
	delete [] m_historyItems;
	delete [] m_iHistoryItemsAuxBuffer;
	delete m_lmLookAhead;
	// word-graph generation?
	if (m_bLatticeGeneration) {
		delete [] m_wshashEntries;
		delete [] m_wgTokens;
		delete [] m_iWordSequenceAux;
	}

	m_bInitialized = false;
}

// begin utterance
void DynamicDecoderX::beginUtterance() {

	assert(m_bInitialized);

	// mark all the nodes as inactive
	for(int i=0 ; i < m_iNodes; ++i) {
		m_nodes[i].iActiveTokensCurrent = 0;
		m_nodes[i].iActiveTokensCurrentBase = -1;
		m_nodes[i].iActiveTokensNext = 0;
		m_nodes[i].iActiveTokensNextBase = -1;
	}
	
	// reset history items
	for(unsigned int i=0 ; i < m_iHistoryItems-1 ; ++i) {
		m_historyItems[i].iActive = -1;
		m_historyItems[i].iPrev = i+1;
	}
	m_historyItems[m_iHistoryItems-1].iActive = -1;
	m_historyItems[m_iHistoryItems-1].iPrev = -1;
	m_iHistoryItemAvailable = 0;
	m_iTimeGarbageCollectionLast = -1;	
	m_iHistoryItemBegSentence = -1;
	
	// active tokens
	m_iActiveTokenTables = 0;
	m_iTokensNext = 0;
	
	m_iNodesActiveCurrent = 0;
	m_iNodesActiveNext = 0;
	
	m_fScoreBest = -FLT_MAX;
	m_fScoreBestWE = -FLT_MAX;
	
	// utterance information
	m_iFeatureVectorsUtterance = 0;
	
	m_hmmManager->resetHMMEmissionProbabilityComputation();
	
	// lattice generation
	if (m_bLatticeGeneration) {
		
		// reset word-graph tokens
		for(unsigned int i=0 ; i < m_iWGTokens-1 ; i += m_iMaxWordSequencesState) {
			m_wgTokens[i].iActive = -1;
			m_wgTokens[i].iPrev = i+m_iMaxWordSequencesState;
		}
		m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iActive = -1;
		m_wgTokens[m_iWGTokens-m_iMaxWordSequencesState].iPrev = -1;
		m_iWGTokenAvailable = 0;	
		
		// set all the entries as "old" (hash table of word-sequences)
		for(unsigned int i=0 ; i < m_iWSHashEntries ; ++i) {
			m_wshashEntries[i].iTime = -1;	
		}	
		m_iWSHashEntryCollisionAvailable = m_iWSHashBuckets;
	}	
}
		
// end utterance
void DynamicDecoderX::endUtterance() {

	if (m_bLatticeGeneration) {
		for(unsigned int i=0 ; i < m_iWSHashEntries ; ++i) {
			if (m_wshashEntries[i].iTime != -1) {
				assert(m_wshashEntries[i].iLexUnit);
				delete [] m_wshashEntries[i].iLexUnit;
			}
		}
	}

	assert(m_bInitialized);
}

// process input feature vectors
void DynamicDecoderX::process(MatrixBase<float> &mFeatures) {

	assert(m_bInitialized);
	assert(mFeatures.getRows() > 0);
	
	// Viterbi search
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	unsigned int t = 0;
	if (m_iFeatureVectorsUtterance == 0) {
		m_iTimeCurrent = 0;
		// root node expansion	
		VectorStatic<float> vFeatureVector = mFeatures.getRow(0);
		expandRoot(vFeatureVector);
		// output search status
		BVC_VERB << "t= " << setw(5) << m_iTimeCurrent << " nodes= " << setw(6) << m_iNodesActiveCurrent << 
			" bestScore= " << FLT(12,4) << m_fScoreBest << " we: " << FLT(12,4) << m_fScoreBestWE;
		++m_iTimeCurrent;
		++t;
	}
	
	// regular expansion
	for( ; t < mFeatures.getRows() ; ++t, ++m_iTimeCurrent) {	
	
		// prune active nodes/tokens
		pruning();
		
		// next frame expansion
		VectorStatic<float> vFeatureVector = mFeatures.getRow(t);
		expand(vFeatureVector,m_iTimeCurrent);	
		
		// output search status
		if (m_iTimeCurrent % 100 == 0) {
			BVC_VERB << "t= " << setw(5) << m_iTimeCurrent << " nodes= " << setw(6) << m_iNodesActiveCurrent << 
				" bestScore= " << FLT(12,4) << m_fScoreBest << " we: " << FLT(12,4) << m_fScoreBestWE;
			//printHashsStats();
			//printHashContents();
		}
	}
	
	m_iFeatureVectorsUtterance += mFeatures.getRows();	
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	
	BVC_VERB << "decoding time: " << FLT(8,2) << dTimeSeconds << " seconds (RTF: " <<  
		FLT(5,2) << dTimeSeconds/(((float)mFeatures.getRows())/100.0) << ")";	
}

// root-node expansion
void DynamicDecoderX::expandRoot(VectorBase<float> &vFeatureVector) {

	DNode *nodeRoot = m_dynamicNetwork->getRootNode();
	float fScore;
	
	int iLMStateInitial = m_lmFSM->getInitialState();
	
	// create the <s> history item
	m_iHistoryItemBegSentence = newHistoryItem();
	HistoryItem *historyItemBegSentence = m_historyItems+m_iHistoryItemBegSentence;
	historyItemBegSentence->iLexUnitPron = m_lexiconManager->m_lexUnitBegSentence->iLexUnitPron;
	historyItemBegSentence->iEndFrame = INT_MIN;
	historyItemBegSentence->fScore = 0.0;
	historyItemBegSentence->iPrev = -1;
	historyItemBegSentence->iActive = -1;
	historyItemBegSentence->iWGToken = -1;	
	
	// expand the root node
	DArc *arcEnd = m_arcs+(nodeRoot+1)->iArcNext;
	for(DArc *arc = m_arcs+nodeRoot->iArcNext ; arc != arcEnd ; ++arc) {
		assert(arc->iType == ARC_TYPE_HMM);
		
		DNode *nodeDest = m_nodes+arc->iNodeDest;
			
		// compute emission probability
		fScore = arc->state->computeEmissionProbability(vFeatureVector.getData(),0);	
		
		// apply insertion-penalty
		fScore += m_dynamicNetwork->getIP((m_nodes+arc->iNodeDest)->iIPIndex);	
		
		if (fScore > m_fScoreBest-m_fBeamWidthNodes) {
		
			if (fScore > m_fScoreBest) {
				m_fScoreBest = fScore;
			}	
		
			int iToken = newToken();
			Token *token = m_tokensNext+iToken;
			token->fScore = fScore;
			token->state = arc->state;
			token->iLMState = iLMStateInitial;
			token->iLexUnitPron = m_iLexUnitPronUnknown;
			token->iNode = arc->iNodeDest;
			token->iHistoryItem = (int)(historyItemBegSentence-m_historyItems);
			token->iLANode = -1;
			token->fLAScores = NULL;
			
			// word graph generation?
			if (m_bLatticeGeneration) {	
				// create the root word-graph token
				token->iWGToken = newWGToken();	
				(token->iWGToken+m_wgTokens)[0].iWordSequence = -2;
				(token->iWGToken+m_wgTokens)[0].iLexUnitPron = m_iLexUnitPronUnknown;
				(token->iWGToken+m_wgTokens)[0].fScore = fScore;
				(token->iWGToken+m_wgTokens)[0].iHistoryItem = (int)(historyItemBegSentence-m_historyItems);
				(token->iWGToken+m_wgTokens)[1].iWordSequence = -1;
			} else {
				token->iWGToken = -1;
			}
			
			// activate the token
			assert(nodeDest->iActiveTokensNextBase == -1);
			nodeDest->iActiveTokensNextBase = newActiveTokenTable();
			(m_activeTokenNext+nodeDest->iActiveTokensNextBase)[0].iLMState = token->iLMState;
			(m_activeTokenNext+nodeDest->iActiveTokensNextBase)[0].iToken = iToken;
			nodeDest->iActiveTokensNext = 1;	
			m_nodesActiveNext[m_iNodesActiveNext++] = nodeDest;
			assert(m_iNodesActiveNext < m_iNodesActiveNextMax);
		}
	}	
}

// regular expansion
void DynamicDecoderX::expand(VectorBase<float> &vFeatureVector, int t) {

	m_fScoreBest = -FLT_MAX;
	m_fScoreBestWE = -FLT_MAX;
	
	float fScore;
	float fScoreToken;
	
	// expand nodes in the active states
	
	// (1) self-loop (hmm is in the token) (no token recombination is needed since all LM-states in the arc are different)
	for(int l=0 ; l < m_iNodesActiveCurrent ; ++l) {
		
		DNode *node = m_nodesActiveCurrent[l];
		ActiveToken *activeTokensCurrent = m_activeTokenCurrent+node->iActiveTokensCurrentBase;
		HMMStateDecoding *state = (m_tokensCurrent+activeTokensCurrent[0].iToken)->state;
		
		// (1) self loop (the hmm-state is in the token)
		
		// compute emission probability
		fScore = state->computeEmissionProbability(vFeatureVector.getData(),t);	
	
		// regular-node
		float *fScoreBest = &m_fScoreBest;
		float fBeamWidth = m_fBeamWidthNodes;
		// we-node
		/*if (node->bWordEnd) {
			fScoreBest = &m_fScoreBestWE;
			fBeamWidth = m_fBeamWidthArcsWE;
		}*/	
		
		// propagate tokens within the arc
		for(int i=0 ; i < node->iActiveTokensCurrent ; ++i) {	
		
			Token *token = m_tokensCurrent+(activeTokensCurrent+i)->iToken;
			fScoreToken = token->fScore+fScore;
			
			if (fScoreToken > (*fScoreBest-fBeamWidth)) {
			
				// keep higher score
				if (fScoreToken > *fScoreBest) {
					*fScoreBest = fScoreToken;
				}
				
				// create expanded token
				int iToken = newToken();
				Token *tokenAux = m_tokensNext+iToken;
				tokenAux->fScore = fScoreToken;
				tokenAux->state = token->state;
				tokenAux->iLMState = token->iLMState;
				tokenAux->iLexUnitPron = token->iLexUnitPron;
				tokenAux->iNode = (int)(node-m_nodes);
				tokenAux->iHistoryItem = token->iHistoryItem;
				tokenAux->iLANode = token->iLANode;
				tokenAux->fLAScores = token->fLAScores;
				
				// word-graph generation?
				if (m_bLatticeGeneration) {
				 	assert(token->iWGToken != -1);
				 	tokenAux->iWGToken = newWGToken(token->iWGToken,fScore);	
				} else {
					tokenAux->iWGToken = -1;	
				}
				
				// activate the arc
				if (node->iActiveTokensNextBase == -1) {
					node->iActiveTokensNextBase = newActiveTokenTable();
					node->iActiveTokensNext = 0;
					m_nodesActiveNext[m_iNodesActiveNext++] = node;
					assert(m_iNodesActiveNext < m_iNodesActiveNextMax);
				}
				(m_activeTokenNext+node->iActiveTokensNextBase)[node->iActiveTokensNext].iLMState = token->iLMState;
				(m_activeTokenNext+node->iActiveTokensNextBase)[node->iActiveTokensNext].iToken = iToken;
				node->iActiveTokensNext++;
				assert(node->iActiveTokensNext < m_iTokensNodeMax);
			}
		}
	}
		
	// (2) outgoing arcs
	for(int i=0 ; i < m_iNodesActiveCurrent ; ++i) {
	
		DNode *node = m_nodesActiveCurrent[i];
		
		// word-end? if so, get ready to extend word history
		m_iHistoryItemsAux = NULL;
		if (node->bWordEnd) {
			m_iHistoryItemsAuxSize = node->iActiveTokensCurrent;
			assert(node->iActiveTokensCurrent < m_iTokensNodeMax);
			m_iHistoryItemsAux = m_iHistoryItemsAuxBuffer;
			for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
				m_iHistoryItemsAux[j] = -1;
			}
			if (m_bLatticeGeneration) {
				for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
					m_iWordSequenceAux[j] = -1;
				}
			}
		}	
			
		DArc *arcEnd = m_arcs+(node+1)->iArcNext;
		for(DArc *arcNext = m_arcs+node->iArcNext ; arcNext != arcEnd ; ++arcNext) {
			
			// hmm-arc
			if (arcNext->iType == ARC_TYPE_HMM) {
				expandToHMM(node,arcNext,vFeatureVector,t);	
			} 
			// word-arc
			else if (arcNext->iType == ARC_TYPE_WORD) {	
			
				// lm-transition
				LMTransition *lmTransition = new LMTransition[node->iActiveTokensCurrent];
				for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
					lmTransition[j].iLMState = -1;
				}	
							
				DNode *node2 = m_nodes+arcNext->iNodeDest;
				DArc *arcEnd2 = m_arcs+(node2+1)->iArcNext;
				for(DArc *arcNext2 = m_arcs+node2->iArcNext ; arcNext2 != arcEnd2 ; ++arcNext2) {
				
					// hmm-arc
					if (arcNext2->iType == ARC_TYPE_HMM) {	
						expandToHMMNewWord(node,arcNext2,arcNext->lexUnit,lmTransition,vFeatureVector,t);
					} 
					// null-arc
					else {
						assert(arcNext2->iType == ARC_TYPE_NULL);
						
						DNode *node3 = m_nodes+arcNext2->iNodeDest;	
						DArc *arcEnd3 = m_arcs+(node3+1)->iArcNext;
						for(DArc *arcNext3 = m_arcs+node3->iArcNext ; arcNext3 != arcEnd3 ; ++arcNext3) {
							// hmm-arc
							assert(arcNext3->iType == ARC_TYPE_HMM);
							expandToHMMNewWord(node,arcNext3,arcNext->lexUnit,lmTransition,vFeatureVector,t);
						}	
					}
				}
				
				// there can be multiple words (homophones) getting to a starting node
				if (node->bWordEnd) {
					for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
						m_iHistoryItemsAux[j] = -1;
					}
					if (m_bLatticeGeneration) {
						for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
							m_iWordSequenceAux[j] = -1;
						}
					}
				}	
				
				delete [] lmTransition;
			}
			// null-arc
			else {
				assert(arcNext->iType == ARC_TYPE_NULL);
			
				DNode *node2 = m_nodes+arcNext->iNodeDest;	
				DArc *arcEnd2 = m_arcs+(node2+1)->iArcNext;
				for(DArc *arcNext2 = m_arcs+node2->iArcNext ; arcNext2 != arcEnd2 ; ++arcNext2) {
			
					// hmm-arc
					if (arcNext2->iType == ARC_TYPE_HMM) {	
						expandToHMM(node,arcNext2,vFeatureVector,t);
					} 
					// word-arc
					else {
						assert(arcNext2->iType == ARC_TYPE_WORD);	
						
						// lm-transition
						LMTransition *lmTransition = new LMTransition[node->iActiveTokensCurrent];
						for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
							lmTransition[j].iLMState = -1;
						}	
						
						DNode *node3 = m_nodes+arcNext2->iNodeDest;	
						DArc *arcEnd3 = m_arcs+(node3+1)->iArcNext;
						for(DArc *arcNext3 = m_arcs+node3->iArcNext ; arcNext3 != arcEnd3 ; ++arcNext3) {
							
							// hmm-arc
							if (arcNext3->iType == ARC_TYPE_HMM) {
								expandToHMMNewWord(node,arcNext3,arcNext2->lexUnit,lmTransition,vFeatureVector,t);	
							} 
							// null-arc
							else {
								assert(arcNext3->iType == ARC_TYPE_NULL);
							
								DNode *node4 = m_nodes+arcNext3->iNodeDest;	
								DArc *arcEnd4 = m_arcs+(node4+1)->iArcNext;
								for(DArc *arcNext4 = m_arcs+node4->iArcNext ; arcNext4 != arcEnd4 ; ++arcNext4) {
									// hmm-arc
									assert(arcNext3->iType == ARC_TYPE_HMM);
									expandToHMMNewWord(node,arcNext4,arcNext2->lexUnit,lmTransition,vFeatureVector,t);	
								}
							}	
						}
						
						// there can be multiple words (homophones) getting to a starting node
						if (node->bWordEnd) {
							for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
								m_iHistoryItemsAux[j] = -1;
							}
							if (m_bLatticeGeneration) {
								for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
									m_iWordSequenceAux[j] = -1;
								}
							}
						}	
						
						delete [] lmTransition;
					}
				}	
			}
		}
		
		node->iActiveTokensCurrent = 0;
		node->iActiveTokensCurrentBase = -1;
	}
}

// expand a series of tokens to a hmm-state
void DynamicDecoderX::expandToHMM(DNode *node, DArc *arcNext, VectorBase<float> &vFeatureVector, int t) {

	float fScore;
	float fScoreToken;
	ActiveToken *activeTokensCurrent = m_activeTokenCurrent+node->iActiveTokensCurrentBase;
	DNode *nodeNext = m_nodes+arcNext->iNodeDest;
	ActiveToken *activeTokensNext = m_activeTokenNext+nodeNext->iActiveTokensNextBase;

	// compute emission probability
	fScore = arcNext->state->computeEmissionProbability(vFeatureVector.getData(),t);	

	bool bWordEnd = (nodeNext->iIPIndex != -1);

	// apply insertion-penalty?
	if (bWordEnd) {
		fScore += m_dynamicNetwork->getIP(m_nodes[arcNext->iNodeDest].iIPIndex);
	}

	// regular-node
	float *fScoreBest = &m_fScoreBest;
	float fBeamWidth = m_fBeamWidthNodes;
	// we-node
	/*if (nodeNext->bWordEnd) {
		fScoreBest = &m_fScoreBestWE;
		fBeamWidth = m_fBeamWidthNodesWE;
	}*/
	
	// propagate tokens within the arc
	for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
	
		Token *token = m_tokensCurrent+(activeTokensCurrent+j)->iToken;
		fScoreToken = token->fScore+fScore;
		
		// compute look-ahead scores?
		if ((arcNext->iLANode != -1) && (token->iLANode == -1)) {
			token->iLANode = arcNext->iLANode;
			//m_lmLookAhead->getLAScores(token->iLMState);
		}	
		
		if (fScoreToken > (*fScoreBest-fBeamWidth)) {
		
			// keep higher score
			if (fScoreToken > *fScoreBest) {
				*fScoreBest = fScoreToken;
			}
			
			// token recombination?
			bool bFound = false;
			for(int k=0 ; k < nodeNext->iActiveTokensNext ; ++k) {
				if (activeTokensNext[k].iLMState == token->iLMState) {
					Token *tokenRec = m_tokensNext+activeTokensNext[k].iToken;
					// no lattice-generation
					if (m_bLatticeGeneration == false) {	
						if (fScoreToken > tokenRec->fScore) {
							// recombine
							tokenRec->fScore = fScoreToken;
							tokenRec->iLexUnitPron = token->iLexUnitPron;
							//assert(tokenRec->state == arcNext->state);
							// (a) not an end-of-word
							//if (m_iHistoryItemsAux == NULL) {
							if (bWordEnd == false) {
								tokenRec->iHistoryItem = token->iHistoryItem;
								tokenRec->iWGToken = -1;	
							} 
							// (b) end-of-word
							else {
								assert(m_iHistoryItemsAux != NULL);
								if (m_iHistoryItemsAux[j] == -1) {
									m_iHistoryItemsAux[j] = newHistoryItem();
									HistoryItem &historyItem = m_historyItems[m_iHistoryItemsAux[j]];
									historyItem.iLexUnitPron = token->iLexUnitPron;
									historyItem.iEndFrame = t-1;
									historyItem.fScore = token->fScore;
									historyItem.iPrev = token->iHistoryItem;
									assert(m_historyItems[historyItem.iPrev].iEndFrame < historyItem.iEndFrame);
									historyItem.iActive = -1;
									historyItem.iWGToken = -1;
								}
								tokenRec->iHistoryItem = m_iHistoryItemsAux[j];
								tokenRec->iWGToken = -1;
							}
						}
					} 
					// lattice-generation
					else {
						
						int iWGToken;
					
						// (a) not an end-of-word	
						if (bWordEnd == false) {
							assert(token->iWGToken != -1);
							iWGToken = newWGToken(token->iWGToken,fScore);	
						} 
						// (b) end-of-word
						else {
							if (m_iHistoryItemsAux[j] == -1) {	
								m_iHistoryItemsAux[j] = newHistoryItem();
								HistoryItem &historyItem = m_historyItems[m_iHistoryItemsAux[j]];
								historyItem.iLexUnitPron = token->iLexUnitPron;
								historyItem.iEndFrame = t-1;
								historyItem.fScore = token->fScore;
								historyItem.iPrev = token->iHistoryItem;
								assert(m_historyItems[historyItem.iPrev].iEndFrame < historyItem.iEndFrame);
								historyItem.iActive = -1;
								historyItem.iWGToken = token->iWGToken;
								assert(m_wgTokens[historyItem.iWGToken].iLexUnitPron != m_iLexUnitPronUnknown);
								m_iWordSequenceAux[j] = hashWordSequence(&historyItem);
							}
							iWGToken = newWGToken(m_iWordSequenceAux[j],fScoreToken,m_iHistoryItemsAux[j]);
						}
						if (m_historyItems[tokenRec->iHistoryItem].iWGToken != -1) {
							assert(m_historyItems[tokenRec->iHistoryItem].iLexUnitPron == m_wgTokens[m_historyItems[tokenRec->iHistoryItem].iWGToken].iLexUnitPron);
						}
						// merge word sequences arriving to the tokens
						bool bReturn = mergeWordSequences(tokenRec->iWGToken,iWGToken);
						if (!bReturn) {
							tokenRec->iLexUnitPron = token->iLexUnitPron;
						}
						// update the score based on the result from the merging process
						tokenRec->fScore = (tokenRec->iWGToken+m_wgTokens)[0].fScore;
						tokenRec->iHistoryItem = (tokenRec->iWGToken+m_wgTokens)[0].iHistoryItem;
						deleteWGToken(iWGToken);
					}
					bFound = true;
					break;
				}
			}
			if (bFound) {
				continue;
			}
		
			// create expanded token (no recombination)
			int iToken = newToken();
			Token *tokenAux = m_tokensNext+iToken;
			tokenAux->fScore = fScoreToken;
			tokenAux->state = arcNext->state;
			tokenAux->iLMState = token->iLMState;
			tokenAux->iLexUnitPron = token->iLexUnitPron;
			tokenAux->iNode = arcNext->iNodeDest;
			tokenAux->iLANode = token->iLANode;
			tokenAux->fLAScores = token->fLAScores;
			// (a) not an end-of-word
			if (bWordEnd == false) {
				tokenAux->iHistoryItem = token->iHistoryItem;
				// word-graph generation?
				if (m_bLatticeGeneration) {
					assert(token->iWGToken != -1);
					tokenAux->iWGToken = newWGToken(token->iWGToken,fScore);	
				} else {
					tokenAux->iWGToken = -1;	
				}
			} 
			// (b) end-of-word
			else {
				assert(m_iHistoryItemsAux != NULL);	
				if (m_iHistoryItemsAux[j] == -1) {
					m_iHistoryItemsAux[j] = newHistoryItem();
					HistoryItem &historyItem = m_historyItems[m_iHistoryItemsAux[j]];
					historyItem.iLexUnitPron = token->iLexUnitPron;
					historyItem.iEndFrame = t-1;
					historyItem.fScore = token->fScore;
					historyItem.iPrev = token->iHistoryItem;
					assert(m_historyItems[historyItem.iPrev].iEndFrame < historyItem.iEndFrame);
					historyItem.iActive = -1;
					// word-graph generation?
					if (m_bLatticeGeneration) {
						historyItem.iWGToken = token->iWGToken;
						m_iWordSequenceAux[j] = hashWordSequence(&historyItem);
					} else {
						historyItem.iWGToken = -1;
					}
				}
				tokenAux->iHistoryItem = m_iHistoryItemsAux[j];
				// word-graph generation?
				if (m_bLatticeGeneration) {
					tokenAux->iWGToken = newWGToken(m_iWordSequenceAux[j],fScoreToken,m_iHistoryItemsAux[j]);
				} else {
					tokenAux->iWGToken = -1;
				}
			}
			
			// activate the node
			if (nodeNext->iActiveTokensNextBase == -1) {
				nodeNext->iActiveTokensNextBase = newActiveTokenTable();
				nodeNext->iActiveTokensNext = 0;
				activeTokensNext = m_activeTokenNext+nodeNext->iActiveTokensNextBase;
				m_nodesActiveNext[m_iNodesActiveNext++] = nodeNext;	
				assert(m_iNodesActiveNext < m_iNodesActiveNextMax);
			}
			activeTokensNext[nodeNext->iActiveTokensNext].iLMState = token->iLMState;
			activeTokensNext[nodeNext->iActiveTokensNext].iToken = iToken;
			nodeNext->iActiveTokensNext++;
			if (nodeNext->iActiveTokensNext >= m_iTokensNodeMax) {
				pruneExtraTokens(nodeNext);
			}
			assert(nodeNext->iActiveTokensNext < m_iTokensNodeMax);
		}
	}
}

// expand a series of tokens to a hmm-state after obsering a new word
void DynamicDecoderX::expandToHMMNewWord(DNode *node, DArc *arcNext, LexUnit *lexUnit, LMTransition *lmTransition, 
	VectorBase<float> &vFeatureVector, int t) {

	float fScore;
	float fScoreLM;
	float fScoreToken;
	ActiveToken *activeTokensCurrent = m_activeTokenCurrent+node->iActiveTokensCurrentBase;
	DNode *nodeNext = m_nodes+arcNext->iNodeDest;
	ActiveToken *activeTokensNext = m_activeTokenNext+nodeNext->iActiveTokensNextBase;
	
	// compute emission probability
	fScore = arcNext->state->computeEmissionProbability(vFeatureVector.getData(),t);	

	bool bWordEnd = (nodeNext->iIPIndex != -1);

	// apply insertion-penalty?
	if (bWordEnd) {
		fScore += m_dynamicNetwork->getIP(m_nodes[arcNext->iNodeDest].iIPIndex);
	}

	// regular-arc
	float *fScoreBest = &m_fScoreBest;
	float fBeamWidth = m_fBeamWidthNodes;
	// we-arc
	if (nodeNext->bWordEnd) {
		//fScoreBest = &m_fScoreBestWE;
		//fBeamWidth = m_fBeamWidthNodesWE;
	}

	// propagate tokens within the node
	for(int j=0 ; j < node->iActiveTokensCurrent ; ++j) {
	
		Token *token = m_tokensCurrent+activeTokensCurrent[j].iToken;
		int iLMState = -1;

		// standard lexical unit: apply LM-score and insertion penalty
		if (m_lexiconManager->isStandard(lexUnit)) {
			// compute LM-score if needed
			if (lmTransition[j].iLMState == -1) {
				lmTransition[j].iLMState = m_lmFSM->updateLMState(token->iLMState,
					lexUnit->iLexUnit,&lmTransition[j].fScoreLM);
				lmTransition[j].fScoreLM *= m_fLMScalingFactor;
			}
			iLMState = lmTransition[j].iLMState;
			fScoreLM = lmTransition[j].fScoreLM;
		} 
		// filler lexical unit: keep original lm state 
		else {
			iLMState = token->iLMState;
			fScoreLM = 0.0;
		}
		
		fScoreToken = token->fScore+fScore+fScoreLM;
		if (fScoreToken > (*fScoreBest-fBeamWidth)) {
		
			// keep higher score
			if (fScoreToken > *fScoreBest) {
				*fScoreBest = fScoreToken;
			}
			
			// token recombination?
			bool bFound = false;
			for(int k=0 ; k < nodeNext->iActiveTokensNext ; ++k) {
				if (activeTokensNext[k].iLMState == iLMState) {
					Token *tokenRec = m_tokensNext+activeTokensNext[k].iToken;
					assert(tokenRec->state == arcNext->state);

					// no lattice-generation
					if (m_bLatticeGeneration == false) {
						
						if (fScoreToken > tokenRec->fScore) {
							// recombine
							tokenRec->fScore = fScoreToken;
							tokenRec->iLexUnitPron = lexUnit->iLexUnitPron;
							// (a) not an end-of-word
							//if (m_iHistoryItemsAux == NULL) {
							if (bWordEnd == false) {
								tokenRec->iHistoryItem = token->iHistoryItem;
								tokenRec->iWGToken = -1;	
							} 
							// (b) end-of-word
							else {
								assert(m_iHistoryItemsAux != NULL);
								if (m_iHistoryItemsAux[j] == -1) {
									m_iHistoryItemsAux[j] = newHistoryItem();
									HistoryItem &historyItem = m_historyItems[m_iHistoryItemsAux[j]];	
									historyItem.iLexUnitPron = lexUnit->iLexUnitPron;
									historyItem.iEndFrame = t-1;
									historyItem.fScore = token->fScore+fScoreLM;
									historyItem.iPrev = token->iHistoryItem;
									assert(m_historyItems[historyItem.iPrev].iEndFrame < historyItem.iEndFrame);
									historyItem.iActive = -1;
									historyItem.iWGToken = -1;
								}
								tokenRec->iHistoryItem = m_iHistoryItemsAux[j];
								tokenRec->iWGToken = -1;
							}
						}
					} 
					// lattice-generation
					else {
						
						int iWGToken;
					
						// (a) not an end-of-word	
						if (bWordEnd == false) {
							assert(token->iWGToken != -1);
							iWGToken = newWGToken(token->iWGToken,fScore+fScoreLM);
							attachLexUnit(iWGToken,lexUnit);	
						} 
						// (b) end-of-word
						else {
							if (m_iHistoryItemsAux[j] == -1) {
								// (IMP: both newHistoryItem and newWGToken can activate garbage collection and invalidate pointers, but not
								// indices)
								// (IMP: newWGToken must be called before newHistoryItem, since otherwise it could invalidate the historyItem)	
								int iWGTokenAux = newWGToken(token->iWGToken);
								m_iHistoryItemsAux[j] = newHistoryItem();
								// create a copy of the wgToken (token can be expanded to different word ends, e.g. homophonic, prefix words)
								HistoryItem &historyItem = m_historyItems[m_iHistoryItemsAux[j]];
								historyItem.iLexUnitPron = lexUnit->iLexUnitPron;
								historyItem.iEndFrame = t-1;
								historyItem.fScore = token->fScore+fScoreLM;
								historyItem.iPrev = token->iHistoryItem;
								assert(m_historyItems[historyItem.iPrev].iEndFrame < historyItem.iEndFrame);
								historyItem.iActive = -1;	
								historyItem.iWGToken = iWGTokenAux;	
								attachLexUnit(historyItem.iWGToken,lexUnit);
								m_iWordSequenceAux[j] = hashWordSequence(&historyItem);
							}
							iWGToken = newWGToken(m_iWordSequenceAux[j],fScoreToken,m_iHistoryItemsAux[j]);
						}
						if (m_historyItems[tokenRec->iHistoryItem].iWGToken != -1) {
							assert(m_historyItems[tokenRec->iHistoryItem].iLexUnitPron == m_wgTokens[m_historyItems[tokenRec->iHistoryItem].iWGToken].iLexUnitPron);
						}
						// merge word sequences arriving to the tokens
						bool bReturn = mergeWordSequences(tokenRec->iWGToken,iWGToken);
						if (!bReturn) {
							tokenRec->iLexUnitPron = lexUnit->iLexUnitPron;	
						}
						// update the score based on the result from the merging process
						tokenRec->fScore = (tokenRec->iWGToken+m_wgTokens)[0].fScore;
						tokenRec->iHistoryItem = (tokenRec->iWGToken+m_wgTokens)[0].iHistoryItem;
						deleteWGToken(iWGToken);
					}
					bFound = true;
					break;
				}
			}
			if (bFound) {
				continue;
			}
		
			// create expanded token (no recombination)
			int iToken = newToken();
			Token *tokenAux = m_tokensNext+iToken;
			tokenAux->fScore = fScoreToken;
			tokenAux->state = arcNext->state;
			tokenAux->iLMState = iLMState;
			tokenAux->iLexUnitPron = lexUnit->iLexUnitPron;
			tokenAux->iNode = (int)(nodeNext-m_nodes);
			tokenAux->iLANode = -1;
			tokenAux->fLAScores = NULL;
			// (a) not an end-of-word
			if (bWordEnd == false) {
				tokenAux->iHistoryItem = token->iHistoryItem;
				// word-graph generation?
				if (m_bLatticeGeneration) {
					assert(token->iWGToken != -1);
					tokenAux->iWGToken = newWGToken(token->iWGToken,fScore+fScoreLM);	
					attachLexUnit(tokenAux->iWGToken,lexUnit);
				} else {
					tokenAux->iWGToken = -1;	
				}
			} 
			// (b) end-of-word
			else {
				assert(m_iHistoryItemsAux != NULL);
				if (m_iHistoryItemsAux[j] == -1) {
					m_iHistoryItemsAux[j] = newHistoryItem();
					HistoryItem &historyItem = m_historyItems[m_iHistoryItemsAux[j]];
					historyItem.iLexUnitPron = lexUnit->iLexUnitPron;
					historyItem.iEndFrame = t-1;
					historyItem.fScore = token->fScore+fScoreLM;
					historyItem.iPrev = token->iHistoryItem;
					assert(m_historyItems[historyItem.iPrev].iEndFrame < historyItem.iEndFrame);
					historyItem.iActive = -1;
					// word-graph generation?
					if (m_bLatticeGeneration) {
						// create a copy of the wgToken (token can be expanded to different word ends, e.g. homophonic, prefix words)
						// IMP: this assignment to be done in two steps since newWGToken(...) can reallocate m_historyItems
						int iWGTokenAux = newWGToken(token->iWGToken);
						m_historyItems[m_iHistoryItemsAux[j]].iWGToken = iWGTokenAux;
						attachLexUnit(m_historyItems[m_iHistoryItemsAux[j]].iWGToken,lexUnit);
						m_iWordSequenceAux[j] = hashWordSequence(&m_historyItems[m_iHistoryItemsAux[j]]);
					} else {
						historyItem.iWGToken = -1;
					}
				}
				tokenAux->iHistoryItem = m_iHistoryItemsAux[j];
				// word-graph generation?
				if (m_bLatticeGeneration) {
					tokenAux->iWGToken = newWGToken(m_iWordSequenceAux[j],fScoreToken,m_iHistoryItemsAux[j]);
				} else {
					tokenAux->iWGToken = -1;
				}
				assert((tokenAux->iHistoryItem+m_historyItems)->iWGToken >= -1);
			}
				
			// activate the node
			if (nodeNext->iActiveTokensNextBase == -1) {
				assert(nodeNext->iActiveTokensNext == 0);
				nodeNext->iActiveTokensNextBase = newActiveTokenTable();
				nodeNext->iActiveTokensNext = 0;
				activeTokensNext = m_activeTokenNext+nodeNext->iActiveTokensNextBase;
				m_nodesActiveNext[m_iNodesActiveNext++] = nodeNext;
				assert(m_iNodesActiveNext < m_iNodesActiveNextMax);
			}
			activeTokensNext[nodeNext->iActiveTokensNext].iLMState = iLMState;
			activeTokensNext[nodeNext->iActiveTokensNext].iToken = iToken;
			nodeNext->iActiveTokensNext++;
			if (nodeNext->iActiveTokensNext >= m_iTokensNodeMax) {
				pruneExtraTokens(nodeNext);
			}	
			assert(nodeNext->iActiveTokensNext < m_iTokensNodeMax);
		}
	}
}

// pruning (token based)
void DynamicDecoderX::pruningOriginal() {

	// (1) reset active nodes from current time frame
	
	assert(m_fScoreBestWE <= m_fScoreBest);
	
	m_fScoreBestWE = m_fScoreBest; //	HACK
	
	// compute threshold for beam based pruning
	float fThresholdRegular = m_fScoreBest-m_fBeamWidthNodes;
	float fThresholdWE = m_fScoreBestWE-m_fBeamWidthNodesWE;
	int iNumberBins = NUMBER_BINS_HISTOGRAM;
	
	// (1) get the best scoring token within each active arc
	float *fScoreBestNode = new float[m_iNodesActiveNext];
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
		fScoreBestNode[i] = -FLT_MAX;
	}
		
	// get the best token-score for each active arc
	int iWillSurviveRegular = 0;
	int iWillSurviveWE = 0;
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
		// get the best score
		assert(m_nodesActiveNext[i]->iActiveTokensNext > 0);
		for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
			Token *token = m_tokensNext+(m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase)[j].iToken;
			if (token->fScore > fScoreBestNode[i]) {
				fScoreBestNode[i] = token->fScore;
			}
		}
		// keep track of surviving arcs
		if (m_nodesActiveNext[i]->bWordEnd) {
			if (fScoreBestNode[i] > fThresholdWE) {
				++iWillSurviveWE;
			}
		} else {
			if (fScoreBestNode[i] > fThresholdRegular) {
				++iWillSurviveRegular;
			}
		}	
	}
	
	// (2) create the histogram
	//if (iWillSurviveRegular > m_iMaxActiveNodes) {
	
		// (2.1) compute the size of each bin and initialize them
		float fLengthRegular = m_fBeamWidthNodes+1;
		float fLengthWE = m_fBeamWidthNodesWE+1;
		float fBinSizeRegular = ((float)fLengthRegular)/((float)iNumberBins);
		float fBinSizeWE = ((float)fLengthWE)/((float)iNumberBins);
		assert(fBinSizeRegular > 0);
		assert(fBinSizeWE > 0);
		int *iBinsRegular = new int[iNumberBins];
		int *iBinsWE = new int[iNumberBins];
		int iRegular = 0;
		int iWE = 0;
		for(int i = 0 ; i < iNumberBins ; ++i) {
			iBinsRegular[i] = 0;
			iBinsWE[i] = 0;
		}
		// (2.2) fill the bins, the first bin keeps the best tokens
		int iBin;
		float fAuxRegular = ((float)iNumberBins)/fLengthRegular;
		float fAuxWE = ((float)iNumberBins)/fLengthWE;
		for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
			// WE-arc
			if (m_nodesActiveNext[i]->bWordEnd) {	
				for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
					Token *token = m_tokensNext+(m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase)[j].iToken;
					if (token->fScore <= fThresholdWE) {
						continue;
					}	
					iBin = (int)(fabs(token->fScore-m_fScoreBestWE)*fAuxWE);
					assert((iBin >= 0) && (iBin < iNumberBins));
					iBinsWE[iBin]++;
					++iWE;
				}
			}
			// regular arc
			else {
				for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
					Token *token = m_tokensNext+(m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase)[j].iToken;
					if (token->fScore <= fThresholdRegular) {
						continue;
					}	
					iBin = (int)(fabs(token->fScore-m_fScoreBest)*fAuxRegular);
					assert((iBin >= 0) && (iBin < iNumberBins));
					iBinsRegular[iBin]++;
					++iRegular;
				}
			}
		}	
		// (2.3) get the threshold
		int iSurvivorsRegular = 0;
		int iSurvivorsWE = 0;
		float fThresholdHistogramRegular = fThresholdRegular;
		float fThresholdHistogramWE = fThresholdWE;
		for(int i = 0 ; i < iNumberBins-1 ; ++i) {
			iSurvivorsRegular += iBinsRegular[i];
			// this is the cut-off
			if (iSurvivorsRegular >= m_iMaxActiveNodes) {
				fThresholdHistogramRegular = m_fScoreBest-(((float)(i+1))*(fLengthRegular/((float)iNumberBins)));
				break;
			}
		}
		for(int i = 0 ; i < iNumberBins-1 ; ++i) {
			iSurvivorsWE += iBinsWE[i];
			// this is the cut-off
			if (iSurvivorsWE >= m_iMaxActiveNodesWE) {
				fThresholdHistogramWE = m_fScoreBestWE-(((float)(i+1))*(fLengthWE/((float)iNumberBins)));
				break;
			}
		}
		
		delete [] iBinsRegular;
		delete [] iBinsWE;

		fThresholdRegular = max(fThresholdRegular,fThresholdHistogramRegular);
		fThresholdWE = max(fThresholdWE,fThresholdHistogramWE);
		
		fThresholdWE = fThresholdRegular;		// HACK
		
		//printf("survivorsWE:      %12d (%12d) %12.4f\n",iSurvivorsWE,iWE,fThresholdHistogramWE);
		//printf("survivorsRegular: %12d (%12d) %12.4f\n",iSurvivorsRegular,iRegular,fThresholdHistogramRegular);
	//}
	
	// prune arcs
	
	
	// for each active arc
	/*for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
		// get the best score
		float fScoreBestNode = -FLT_MAX;
		for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
			if (m_nodesActiveNext[i]->activeTokensNext[j].token->fScore > fScoreBestNode) {
				fScoreBestNode = m_nodesActiveNext[i]->activeTokensNext[j].token->fScore;
			}
		}
		// should the arc be pruned?
		if (fScoreBestNode < (m_fScoreBest-m_fBeamWidthArcs)) {
			
		}	
		// histogram pruning 
		iBin = (int)(fabs(fScoreBestNode-)*fAux);
		assert((iBin >= 0) && (iBin < iNumberBins));
		iBins[iBin]++;
		
		
		// pruning within the arc
		for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
			if (m_nodesActiveNext[i]->activeTokensNext[j].token->fScore > fScoreBestNode) {
				fScoreBestNode = m_nodesActiveNext[i]->activeTokensNext[j].token->fScore;
			}
		}
	}	*/	
	
	// apply pruning to each arc
	m_iNodesActiveCurrent = 0;
	float fThresholdArc;
	int iPrunedTotal = 0;	
	/*int iSurvivorsH[1000];
	for(int i=0 ; i < 1000 ; ++i) {	
		iSurvivorsH[i] = 0;
	}*/
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {	
	
		DNode *node = m_nodesActiveNext[i];
		assert(node->iActiveTokensNext > 0);
		
		/*if (arc->bWordEnd == false) {
			fThresholdArc = max(fThresholdRegular,fScoreBestNode[i]-m_fBeamWidthTokensArc);		
		} else {
			fThresholdArc = max(fThresholdWE,fScoreBestNode[i]-m_fBeamWidthTokensArc);	
		}
		
		if (fScoreBestNode[i] < fThresholdArc) {
			iPrunedTotal += arc->iActiveTokensNext;	
			assert(arc->activeTokensCurrent == NULL);
			assert(arc->iActiveTokensCurrent == 0);	
			arc->activeTokensNext = NULL;
			arc->iActiveTokensNext = 0;
			continue;
		}*/			
		
		//fThresholdArc = fScoreBestNode[i]-m_fBeamWidthTokensArc;
		
		// compute bin-size for histogram pruning
		/*float fLength = m_fBeamWidthTokensNode+1;
		float fBinSize = ((float)fLength)/((float)iNumberBins);
		assert(fBinSize > 0);
		int iBins[iNumberBins];
		int iBin;
		float fAux = ((float)iNumberBins)/fLength;
		
		// within-arc histogram building
		if (arc->iActiveTokensNext > m_iMaxActiveTokensArc) {
		
			// (2.1) compute the size of each bin and initialize them
			for(int j = 0 ; j < iNumberBins ; ++j) {
				iBins[j] = 0;
			}
			// (2.2) fill the bins, the first bin keeps the best tokens
			for(int j=0 ; j < arc->iActiveTokensNext ; ++j) {
				Token *token = m_tokensNext+arc->activeTokensNext[j].iToken;
				if (token->fScore > fThresholdArc) {
					iBin = (int)(fabs(token->fScore-fScoreBestNode[i])*fAux);
					assert((iBin >= 0) && (iBin < iNumberBins));
					iBins[iBin]++;
				}
			}
			// (2.3) get the threshold
			int iSurvivors = 0;
			float fThresholdHistogram = fThresholdArc;
			for(int j = 0 ; j < iNumberBins-1 ; ++j) {
				iSurvivors += iBins[j];
				// this is the cut-off
				if (iSurvivors >= m_iMaxActiveTokensArc) {
					fThresholdHistogram = fScoreBestNode[i]-(((float)(j+1))*(fLength/((float)iNumberBins)));
					break;
				}
			}
			fThresholdArc = max(fThresholdArc,fThresholdHistogram);
		}*/
		
		if (node->bWordEnd) {
			fThresholdArc = fThresholdWE;
		} else {
			fThresholdArc = fThresholdRegular;
		}		
		
		// actual within-arc token pruning
		int iPruned = 0;
		int iAvailable = -1;
		ActiveToken *activeTokensNext = m_activeTokenNext+node->iActiveTokensNextBase;
		for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
			Token *token = m_tokensNext+activeTokensNext[j].iToken;
			if (token->fScore < fThresholdArc) {	
				if (iAvailable == -1) {
					iAvailable = j;	
				}
				++iPruned;
			} else if (iAvailable != -1) {
				activeTokensNext[iAvailable] = activeTokensNext[j];
				++iAvailable;
			}
		}
		iPrunedTotal += iPruned;
		assert(node->iActiveTokensCurrentBase == -1);
		assert(node->iActiveTokensCurrent == 0);
		int iSurvivors = node->iActiveTokensNext-iPruned;
		//iSurvivorsH[iSurvivors]++;
		if (iSurvivors > 0) {
			assert(node->iActiveTokensNextBase != -1);
			node->iActiveTokensCurrentBase = node->iActiveTokensNextBase;
			node->iActiveTokensCurrent = iSurvivors;
			m_nodesActiveCurrent[m_iNodesActiveCurrent++] = node;
		}
		node->iActiveTokensNextBase = -1;
		node->iActiveTokensNext = 0;
	}
	
	/*int iAcc = 0;
	for(int i=1 ; i < 1000 ; ++i) {	
		iAcc += i*iSurvivorsH[i];
		if (iSurvivorsH[i] != 0) {
			printf("%d -> %d (%d)\n",i,iSurvivorsH[i],iAcc);
		}
	}*/
	
	//printf("# arcs active: (%d -> %d) # tokens active: (%d -> %d)\n",
	//	m_iNodesActiveNext,m_iNodesActiveCurrent,m_iTokensNext,m_iTokensNext-iPrunedTotal);
	
	// clean the table of next active states
	m_iNodesActiveNext = 0;	
	
	// swap the token tables
	swapTokenTables();
	
	delete [] fScoreBestNode;
}


// pruning (token based)
void DynamicDecoderX::pruning() {

	assert(m_fScoreBestWE <= m_fScoreBest);
	
	m_fScoreBestWE = m_fScoreBest; //	HACK
	
	// compute threshold for beam based pruning
	int iNumberBins = NUMBER_BINS_HISTOGRAM;
	int iSurvivorsRegular = 0;
	float fThresholdRegular = m_fScoreBest-m_fBeamWidthNodes;	
	
	// (1) apply within-node pruning
		
	// apply pruning to each arc
	m_iNodesActiveCurrent = 0;
	int iPrunedTotal = 0;	
	int *iBins = new int[iNumberBins];
	int iBin;
	float fScoreBestNode = -FLT_MAX;
	float *fThresholdNode = new float[m_iNodesActiveNext];
	int iSurvivorsAll = 0;
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
	
		DNode *node = m_nodesActiveNext[i];
		assert(node->iActiveTokensNext > 0);
		fThresholdNode[i] = fThresholdRegular;
		
		// (a) histogram pruning within the node
		if (node->iActiveTokensNext > m_iMaxActiveTokensNode) {
		
			// get the best score within the node
			fScoreBestNode = -FLT_MAX;
			for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
				Token *token = m_tokensNext+(m_activeTokenNext+node->iActiveTokensNextBase)[j].iToken;
				if (token->fScore > fScoreBestNode) {
					fScoreBestNode = token->fScore;
				}
			}
		
			// compute bin-size for histogram pruning
			float fLength = m_fBeamWidthTokensNode+1;
			float fBinSize = fLength/((float)iNumberBins);
			assert(fBinSize > 0);
					
			fThresholdNode[i] = max(fScoreBestNode-m_fBeamWidthTokensNode,fThresholdRegular);
		
			// (2.1) compute the size of each bin and initialize them
			for(int j = 0 ; j < iNumberBins ; ++j) {
				iBins[j] = 0;
			}
			// (2.2) fill the bins, the first bin keeps the best tokens
			ActiveToken *activeTokensNext = m_activeTokenNext+node->iActiveTokensNextBase;
			for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
				Token *token = m_tokensNext+activeTokensNext[j].iToken;
				if (token->fScore > fThresholdNode[i]) {
					iBin = (int)(fabs(token->fScore-fScoreBestNode)/fBinSize);
					assert((iBin >= 0) && (iBin < iNumberBins));
					iBins[iBin]++;
				}
			}
			// (2.3) get the threshold
			int iSurvivors = 0;
			float fThresholdHistogram = fThresholdRegular;
			for(int j = 0 ; j < iNumberBins-1 ; ++j) {
				iSurvivors += iBins[j];
				// this is the cut-off
				if (iSurvivors >= m_iMaxActiveTokensNode) {
					fThresholdHistogram = fScoreBestNode-((j+1)*fBinSize);
					iSurvivorsAll += iSurvivors;
					break;
				}
			}
			fThresholdNode[i] = max(fThresholdNode[i],fThresholdHistogram);
		} else {
			iSurvivorsAll += node->iActiveTokensNext;
		}
	}
	
	// (2) global-pruning
	
	// histogram pruning?
	if (iSurvivorsAll > m_iMaxActiveNodes) {
	
		// (2.1) compute the size of each bin and initialize them
		float fLengthRegular = m_fBeamWidthNodes+1;
		float fBinSizeRegular = ((float)fLengthRegular)/((float)iNumberBins);
		assert(fBinSizeRegular > 0);
		int *iBinsRegular = new int[iNumberBins];
		int iRegular = 0;
		for(int i = 0 ; i < iNumberBins ; ++i) {
			iBinsRegular[i] = 0;
		}
		// (2.2) fill the bins, the first bin keeps the best tokens
		int iBin;
		
		for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
			ActiveToken *activeTokens = m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase;
			for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
				Token *token = m_tokensNext+activeTokens[j].iToken;
				if (token->fScore <= fThresholdNode[i]) {
					continue;
				}	
				iBin = (int)(fabs(token->fScore-m_fScoreBest)/fBinSizeRegular);
				assert((iBin >= 0) && (iBin < iNumberBins));
				iBinsRegular[iBin]++;
				++iRegular;
			}
		}	
		// (2.3) get the threshold
		float fThresholdHistogramRegular = fThresholdRegular;
		for(int i = 0 ; i < iNumberBins-1 ; ++i) {
			iSurvivorsRegular += iBinsRegular[i];
			// this is the cut-off
			if (iSurvivorsRegular >= m_iMaxActiveNodes) {
				fThresholdHistogramRegular = m_fScoreBest-((i+1)*fBinSizeRegular);
				break;
			}
		}
		fThresholdRegular = max(fThresholdRegular,fThresholdHistogramRegular);

		delete [] iBinsRegular;
	}	
		
	// (3) actual pruning
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
	
		DNode *node = m_nodesActiveNext[i];
		fThresholdNode[i] = max(fThresholdNode[i],fThresholdRegular);	
		
		int iPruned = 0;
		int iAvailable = -1;
		ActiveToken *activeTokensNext = m_activeTokenNext+node->iActiveTokensNextBase;
		for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
			Token *token = m_tokensNext+activeTokensNext[j].iToken;
			if (token->fScore < fThresholdNode[i]) {
				if (iAvailable == -1) {
					iAvailable = j;	
				}
				++iPruned;
			} else if (iAvailable != -1) {
				activeTokensNext[iAvailable] = activeTokensNext[j];
				++iAvailable;
			}
		}
		iPrunedTotal += iPruned;
		assert(node->iActiveTokensCurrentBase == -1);
		assert(node->iActiveTokensCurrent == 0);
		int iSurvivors = node->iActiveTokensNext-iPruned;
		//iSurvivorsH[iSurvivors]++;
		if (iSurvivors > 0) {
			assert(node->iActiveTokensNextBase != -1);
			node->iActiveTokensCurrentBase = node->iActiveTokensNextBase;
			node->iActiveTokensCurrent = iSurvivors;
			m_nodesActiveCurrent[m_iNodesActiveCurrent++] = node;
		}
		node->iActiveTokensNextBase = -1;
		node->iActiveTokensNext = 0;
	}
	
	delete [] iBins;
	delete [] fThresholdNode;
	
	//printf("# tokens active: (%d -> %d)\n",m_iTokensNext,m_iTokensNext-iPrunedTotal);
	
	// clean the table of next active states
	m_iNodesActiveNext = 0;	
	
	// swap the token tables
	swapTokenTables();
}
		
// return the BestPath
BestPath *DynamicDecoderX::getBestPath() {

	assert(m_bInitialized);

	// (1) get the best scoring token
	float fScoreBestWE = -FLT_MAX;
	Token *tokenBest = NULL;
	LexUnit *lexUnitLast = NULL;
	float fScoreToken;
	float fScoreLM1;
	int iLMState;
	
	// for each active token
	for(int i=0 ; i < m_iTokensNext ; ++i) {
		Token *token = &m_tokensNext[i];
		// check if pruned by "pruneExtraTokens"
		if (token->iNode == -1) {	
			continue;
		}
		DNode *node = m_nodes+token->iNode;
		// only tokens in word-end nodes
		if (node->bWordEnd == false) {
			continue;
		}
		// two possibilities: 
		// (a) monophones: the node goes to a word-arc before any hmm-node is seen 
		// (b) multi-phones: the node goes directly to an hmm-state
		VLexUnit vLexUnitDest;
		getDestinationMonophoneLexUnits(node,vLexUnitDest);
		// (a) monophones
		if (vLexUnitDest.empty() == false) {
			for(VLexUnit::iterator it = vLexUnitDest.begin() ; it != vLexUnitDest.end() ; ++it) {	
				fScoreLM1 = 0.0;
				iLMState = -1;
				if (m_lexiconManager->isStandard(*it)) {
					iLMState = m_lmFSM->updateLMState(token->iLMState,(*it)->iLexUnit,&fScoreLM1);
					fScoreLM1 *= m_fLMScalingFactor;
				} else {
					iLMState = token->iLMState;
				}
				float fScoreLM2 = m_lmFSM->toFinalState(iLMState)*m_fLMScalingFactor;
				fScoreToken = token->fScore + fScoreLM1 + fScoreLM2;
				if (fScoreToken > fScoreBestWE) {
					fScoreBestWE = fScoreToken;
					tokenBest = token;
					lexUnitLast = *it;
				}
			}
		}
		// (b) multi-phones 
		else {
			fScoreLM1 = m_lmFSM->toFinalState(token->iLMState)*m_fLMScalingFactor;
			fScoreToken = token->fScore + fScoreLM1;
			if (fScoreToken > fScoreBestWE) {
				fScoreBestWE = fScoreToken;
				tokenBest = token;
				lexUnitLast = m_lexiconManager->getLexUnitPron(token->iLexUnitPron);
			}
		}
	}
	
	// no active tokens at terminal nodes:
	if (tokenBest == NULL) {
		BVC_WARNING << "unable to retrieve the best decoding path, no active tokens found at terminal nodes";
		// TODO in this scenario, which is usually very rare, it would be possible to generate a best path by 
		// modifying the path in the best token by doing a) or b):
		// a) changing the alignment of the last word/silence so it ends a the last HMM-state (use forced alignment)
		// b) remove last word and run force alignment on the rest in order to compute best path score
		return NULL;
	}	
	
	// (2) get the best sequence of lexical units from the best scoring token
	BestPath *bestPath = new BestPath(m_lexiconManager,fScoreBestWE);
	int iHistoryItem = tokenBest->iHistoryItem;
	int iFrameStart;
	int iFrameEnd;
	float fScore;
	while (iHistoryItem != -1) {
		HistoryItem *historyItem = iHistoryItem+m_historyItems;
		if (historyItem->iPrev != -1) {
			iFrameStart = max(0,m_historyItems[historyItem->iPrev].iEndFrame+1);			// (INT_MIN is used for the initial item)
			iFrameEnd = historyItem->iEndFrame;
			assert(iFrameStart < iFrameEnd);
			fScore = historyItem->fScore;
		} else {
			iFrameStart = -1;
			iFrameEnd = -1;
			fScore = 0.0;
		}
		// get the observed lexical unit
		LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(historyItem->iLexUnitPron);
		// add a new element
		bestPath->newElementFront(iFrameStart,iFrameEnd,fScore,0.0,0.0,0.0,lexUnit,0.0);
		iHistoryItem = m_historyItems[iHistoryItem].iPrev;
	}

	// add the final lexical unit if any
	assert(lexUnitLast != NULL);
	if (m_historyItems[tokenBest->iHistoryItem].iLexUnitPron != m_lexiconManager->m_lexUnitBegSentence->iLexUnitPron) {
		iFrameStart = m_historyItems[tokenBest->iHistoryItem].iEndFrame+1;
	} else {
		iFrameStart = 0;
	}
	bestPath->newElementBack(iFrameStart,m_iFeatureVectorsUtterance-1,tokenBest->fScore,0.0,0.0,0.0,lexUnitLast,0.0);
	
	// add the end of sentence
	bestPath->newElementBack(-1,-1,0.0,0.0,0.0,0.0,m_lexiconManager->m_lexUnitEndSentence,0.0);
	
	// TODO attach lm-scores and compute real am-scores by substracting lm-score and insertion penalty
	
	BVC_VERB	<< "best score:    " << FLT(12,4) << m_fScoreBest;
	BVC_VERB	<< "best WE score: " << FLT(12,4) << fScoreBestWE << " " << FLT(12,4) << m_fScoreBestWE;
	
	return bestPath;
}

// get monophone lexical units accessible right aftet the given hmm-node 
void DynamicDecoderX::getDestinationMonophoneLexUnits(DNode *node, VLexUnit &vLexUnitDest) {

	map<int,bool> mLexUnitSeen;

	DArc *arcEnd = m_arcs+(node+1)->iArcNext;
	for(DArc *arcNext = m_arcs+node->iArcNext ; arcNext != arcEnd ; ++arcNext) {
		
		// word-arc
		if (arcNext->iType == ARC_TYPE_WORD) {
			vLexUnitDest.push_back(arcNext->lexUnit);
		}
		// null-arc
		else if (arcNext->iType == ARC_TYPE_NULL) {
		
			DNode *node2 = m_nodes+arcNext->iNodeDest;	
			DArc *arcEnd2 = m_arcs+(node2+1)->iArcNext;
			for(DArc *arcNext2 = m_arcs+node2->iArcNext ; arcNext2 != arcEnd2 ; ++arcNext2) {	
				if (arcNext2->iType == ARC_TYPE_WORD) {
					if (arcNext2->lexUnit->vPhones.size() == 1) {
						if (mLexUnitSeen.find(arcNext2->lexUnit->iLexUnit) == mLexUnitSeen.end()) {
							//m_lexiconManager->print(arcNext2->lexUnit);	
							mLexUnitSeen.insert(map<int,bool>::value_type(arcNext2->lexUnit->iLexUnit,true));
							vLexUnitDest.push_back(arcNext2->lexUnit);
						}
					}
				}
			}	
		}
	}
}

// garbage collection of history items
// (1) it starts by marking the active items by traversing back items from the active states
// (2) it adds inactive items to the queue of available items
void DynamicDecoderX::historyItemGarbageCollection() {

	int iItemsActive = 0;
	
	// (1) check if garbage collection was already run within the current time frame
	// note: this is an undesirable situation because it requires an extra pass over the complete array of items
	// it should be avoided by allocating a larger number of entries from the beginning
	if (m_iTimeGarbageCollectionLast == m_iTimeCurrent) {
		// mark all the history items as inactive
		for(unsigned int i=0 ; i < m_iHistoryItems ; ++i) {
			m_historyItems[i].iActive = -1;
		}	
	}	
	
	// (2) mark items coming from active arcs as active
	// (2.1) active arcs for current time frame
	for(int i=0 ; i < m_iNodesActiveCurrent ; ++i) {
		ActiveToken *activeTokens = m_activeTokenCurrent+m_nodesActiveCurrent[i]->iActiveTokensCurrentBase;
		for(int j=0 ; j < m_nodesActiveCurrent[i]->iActiveTokensCurrent ; ++j) {
			Token *token = m_tokensCurrent+activeTokens[j].iToken;
			int iHistoryItem = token->iHistoryItem;	
			while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
				(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;	
				iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
				++iItemsActive;
			}
		}
	}	
	// (2.2) active arcs for next time frame
	// TODO: it would be faster to traverse the array of next tokens, but this should not be too bad
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
		ActiveToken *activeTokens = m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase;
		for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
			Token *token = m_tokensNext+activeTokens[j].iToken;
			int iHistoryItem = token->iHistoryItem;	
			while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
				(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;	
				iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
				++iItemsActive;
			}
		}
	}
	// check also auxiliar arrays for history items in use
	if (m_iHistoryItemsAux != NULL) {
		for(int i=0 ; i < m_iHistoryItemsAuxSize ; ++i) {
			if (m_iHistoryItemsAux[i] != -1) {
				if ((m_historyItems+m_iHistoryItemsAux[i])->iActive != m_iTimeCurrent) {	
					(m_historyItems+m_iHistoryItemsAux[i])->iActive = m_iTimeCurrent;
					++iItemsActive;
				}
			}
		}	
	}
	
	// (3) if a certain percentage* of the items are active then we need to allocate 
	// a bigger data structure to keep the new history items
	// (* if we wait until all the items are active there will be many calls to the garbage collector
	// when the array of reach a high occupation, that would introduce substantial overhead)
	assert(iItemsActive <= (int)m_iHistoryItems);
	if (iItemsActive >= (0.20*m_iHistoryItems)) {

		//printf("history item garbage collection...\n");	
		//printf("allocating space for new items (item sused: %d existing: %d)\n",iItemsActive,m_iHistoryItems);	
		
		// allocate a new data structure with double capacity
		HistoryItem *historyItems = NULL;
		try {
			historyItems = new HistoryItem[m_iHistoryItems*2];
		} 
		catch (const std::bad_alloc&) {
			int iBytes = m_iHistoryItems*2*sizeof(WGToken);
			BVC_ERROR << "unable to allocate memory for history items, " << iBytes << " Bytes needed";
		}
		
		// copy the active items from the old data structure
		for(unsigned int i=0 ; i < m_iHistoryItems ; ++i) {
			historyItems[i].iLexUnitPron = m_historyItems[i].iLexUnitPron;
			historyItems[i].iEndFrame = m_historyItems[i].iEndFrame;
			historyItems[i].fScore = m_historyItems[i].fScore;
			historyItems[i].iActive = m_iTimeCurrent;
			historyItems[i].iWGToken = m_historyItems[i].iWGToken;
			historyItems[i].iPrev = m_historyItems[i].iPrev;
		}
		
		// create the linked list of available items
		for(unsigned int i=m_iHistoryItems ; i < (2*m_iHistoryItems)-1 ; ++i) {
			historyItems[i].iPrev = i+1;
			historyItems[i].iActive = -1;	
		}
		historyItems[(2*m_iHistoryItems)-1].iPrev = -1;
		historyItems[(2*m_iHistoryItems)-1].iActive = -1;
		
		delete [] m_historyItems;
		m_historyItems = historyItems;
		m_iHistoryItemAvailable = m_iHistoryItems;
		m_iHistoryItems *= 2;
	}
	// (3') there are inactive items: create a linked list with them
	else {
		int *iHistoryItemAux = &m_iHistoryItemAvailable;
		for(unsigned int i = 0 ; i < m_iHistoryItems ; ++i) {
			if (m_historyItems[i].iActive != m_iTimeCurrent) {
				m_historyItems[i].iActive = -1;
				m_historyItems[i].iEndFrame = -1;
				*iHistoryItemAux = i;
				iHistoryItemAux = &m_historyItems[i].iPrev;	
			}
		}
		*iHistoryItemAux = -1;
	}
	
	m_iTimeGarbageCollectionLast = m_iTimeCurrent;
	
	//printf("%d used %d total\n",iItemsActive,m_iHistoryItems);	
}
		
// marks unused history items as available (lattice generation)
void DynamicDecoderX::historyItemGarbageCollectionLattice(bool bRecycleHistoryItems, bool bRecycleWGTokens) {

	unsigned int iItemsActive = 0;
	unsigned int iTokensActive = 0;
	assert(bRecycleHistoryItems != bRecycleWGTokens);
	
	//printf("doing garbage collection\n");
	
	// (0) check if garbage collection was already run within the current time frame
	// note: this is an undesirable situation because it requires an extra pass over the complete array of items
	// it should be avoided by allocating a larger number of entries from the beginning
	if (m_iTimeGarbageCollectionLast == m_iTimeCurrent) {
		// mark all the history items as inactive
		for(unsigned int i=0 ; i < m_iHistoryItems ; ++i) {
			m_historyItems[i].iActive = -1;
		}	
		// mark all the word-graph tokens as inactive
		for(unsigned int i=0 ; i < m_iWGTokens ; i += m_iMaxWordSequencesState) {
			m_wgTokens[i].iActive = -1;
		}			
	}
		
	bool *bTokensActive = new bool[m_iWGTokens];
	for(unsigned int i=0 ; i < m_iWGTokens ; ++i) {
		bTokensActive[i] = false;
	}
	
	// array to keep the active history items
	int iHistoryItemActiveSize = 0;
	int *iHistoryItemActive = new int[m_iHistoryItems];
	
	// (1) mark items coming from active arcs as active
	// (1.1) active arcs for current time frame
	for(int i=0 ; i < m_iNodesActiveCurrent ; ++i) {
		ActiveToken *activeTokens = m_activeTokenCurrent+m_nodesActiveCurrent[i]->iActiveTokensCurrentBase;
		for(int j=0 ; j < m_nodesActiveCurrent[i]->iActiveTokensCurrent ; ++j) {
			Token *token = m_tokensCurrent+activeTokens[j].iToken;
			int iHistoryItem = token->iHistoryItem;	
			while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
				(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;
				assert((m_historyItems+iHistoryItem)->iWGToken >= -1);
				iHistoryItemActive[iHistoryItemActiveSize++] = iHistoryItem;
				iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
				++iItemsActive;
			}
				
			// mark the wg-token as active
			assert(token->iWGToken != -1);
			(token->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
			assert((bTokensActive[token->iWGToken] == false));
			bTokensActive[token->iWGToken] = true;
			++iTokensActive;
			
			// keep history items at the current wg-token
			for(int i=0 ; (i < m_iMaxWordSequencesState) && ((token->iWGToken+m_wgTokens)[i].iWordSequence != -1) ; ++i) {
				HistoryItem *historyItem = m_historyItems+(token->iWGToken+m_wgTokens)[i].iHistoryItem;
				if (historyItem->iActive != m_iTimeCurrent) {
					historyItem->iActive = m_iTimeCurrent;
					assert(historyItem->iWGToken >= -1);
					iHistoryItemActive[iHistoryItemActiveSize++] = (int)(historyItem-m_historyItems);
					++iItemsActive;
				}
			}
		}
	}	
	
	//printf("(1) # active wg-tokens: %d\n",iTokensActive);	
	
	// (1.1) active arcs for next time frame
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
		ActiveToken *activeTokens = m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase;	
		for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
			Token *token = m_tokensNext+activeTokens[j].iToken;
			int iHistoryItem = token->iHistoryItem;	
			while((iHistoryItem != -1) && ((m_historyItems+iHistoryItem)->iActive != m_iTimeCurrent)) {
				(m_historyItems+iHistoryItem)->iActive = m_iTimeCurrent;	
				assert((m_historyItems+iHistoryItem)->iWGToken >= -1);
				iHistoryItemActive[iHistoryItemActiveSize++] = iHistoryItem;
				iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
				++iItemsActive;
			}
			
			assert(token->iWGToken != -1);
			(token->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
			assert(bTokensActive[token->iWGToken] == false);
			bTokensActive[token->iWGToken] = true;
			++iTokensActive;
			
			// keep history items at the current wg-token	
			for(int i=0 ; (i < m_iMaxWordSequencesState) && ((token->iWGToken+m_wgTokens)[i].iWordSequence != -1) ; ++i) {
				HistoryItem *historyItem = m_historyItems+(token->iWGToken+m_wgTokens)[i].iHistoryItem;
				if (historyItem->iActive != m_iTimeCurrent) {
					historyItem->iActive = m_iTimeCurrent;
					assert(historyItem->iWGToken >= -1);
					iHistoryItemActive[iHistoryItemActiveSize++] = (int)(historyItem-m_historyItems);	
					++iItemsActive;
				}
			}
		}
	}
	
	// check auxiliar arrays for history items in use
	if (m_iHistoryItemsAux != NULL) {
		for(int i=0 ; i < m_iHistoryItemsAuxSize ; ++i) {
			if (m_iHistoryItemsAux[i] != -1) {
				if ((m_historyItems+m_iHistoryItemsAux[i])->iActive != m_iTimeCurrent) {	
					(m_historyItems+m_iHistoryItemsAux[i])->iActive = m_iTimeCurrent;
					iHistoryItemActive[iHistoryItemActiveSize++] = m_iHistoryItemsAux[i];	
					++iItemsActive;
				}
			}
		}	
	}	
	
	// check the list used to generate the lattice, it may have items in use in the rare case that the lattice
	// is being generated
	for(list<int>::iterator it = m_lHistoryItem.begin() ; it != m_lHistoryItem.end() ; ++it) {
		HistoryItem *historyItem = (*it)+m_historyItems;
		historyItem->iActive = m_iTimeCurrent;
		iHistoryItemActive[iHistoryItemActiveSize++] = *it;	
		++iItemsActive;	
		// corresponding wg-token
		int iWGToken = historyItem->iWGToken;
		assert(iWGToken != -1);
		(iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
		assert(bTokensActive[iWGToken] == false);
		bTokensActive[iWGToken] = true;
		++iTokensActive;
		/////////////
		for(int i=0 ; (i < m_iMaxWordSequencesState) && ((historyItem->iWGToken+m_wgTokens)[i].iWordSequence != -1) ; ++i) {
			HistoryItem *historyItem2 = m_historyItems+(historyItem->iWGToken+m_wgTokens)[i].iHistoryItem;
			historyItem2->iActive = m_iTimeCurrent;
			iHistoryItemActive[iHistoryItemActiveSize++] = (historyItem->iWGToken+m_wgTokens)[i].iHistoryItem;	
			++iItemsActive;	
		}	
	}
		
	while(iHistoryItemActiveSize > 0) {
	
		HistoryItem *historyItem = m_historyItems+iHistoryItemActive[iHistoryItemActiveSize-1];
		--iHistoryItemActiveSize;
		
		if (historyItem->iWGToken != -1) {
			
			(historyItem->iWGToken+m_wgTokens)->iActive = m_iTimeCurrent;
			
			if (bTokensActive[historyItem->iWGToken] == false) {
				bTokensActive[historyItem->iWGToken] = true;
				++iTokensActive;	
			
				for(int i=0 ; ((i < m_iMaxWordSequencesState) && 
					((historyItem->iWGToken+m_wgTokens)[i].iWordSequence != -1)) ; ++i) {
				
					HistoryItem *historyItem2 = m_historyItems+(historyItem->iWGToken+m_wgTokens)[i].iHistoryItem;	
					if (historyItem2->iActive != m_iTimeCurrent) {
						historyItem2->iActive = m_iTimeCurrent;
						iHistoryItemActive[iHistoryItemActiveSize++] = (int)(historyItem2-m_historyItems);	
						++iItemsActive;
					}
				}
			}
		}	
	}	
	delete [] iHistoryItemActive;
	delete [] bTokensActive;		
	
	//printf("# active wg-tokens: %d  %d #items: %d\n",iTokensActive,m_iTokensNext,iItemsActive);
	
	if (bRecycleHistoryItems) {
		// (2) if a certain percentage* of the items are active then we need to allocate 
		// a bigger data structure to keep the new history items
		// (* otherwise if we wait until all the items are active there will be many calls to the garbage collector
		// when the array reaches a high occupation, that would introduce substantial overhead)
		assert(iItemsActive <= m_iHistoryItems);
		if (iItemsActive >= (0.20*m_iHistoryItems)) {
	
			//printf("history item garbage collection...\n");	
			//printf("allocating space for new items (item sused: %d existing: %d)\n",iItemsActive,m_iHistoryItems);	
			
			// allocate a new data structure with double capacity	
			HistoryItem *historyItems = NULL;
			try {
				historyItems = new HistoryItem[m_iHistoryItems*2];
			} 
			catch (const std::bad_alloc&) {
				int iBytes = m_iHistoryItems*2*sizeof(WGToken);
				BVC_ERROR << "unable to allocate memory for history items, " << iBytes << " Bytes needed";
			}		
			
			// copy the active items from the old data structure
			for(unsigned int i=0 ; i < m_iHistoryItems ; ++i) {
				historyItems[i].iLexUnitPron = m_historyItems[i].iLexUnitPron;
				historyItems[i].iEndFrame = m_historyItems[i].iEndFrame;
				historyItems[i].fScore = m_historyItems[i].fScore;
				historyItems[i].iActive = m_iTimeCurrent;
				historyItems[i].iWGToken = m_historyItems[i].iWGToken;
				historyItems[i].iPrev = m_historyItems[i].iPrev;
			}
			
			// create the linked list of available items
			int iLastIndex = (2*m_iHistoryItems)-1;
			for(int i=m_iHistoryItems ; i < iLastIndex ; ++i) {
				historyItems[i].iPrev = i+1;
				historyItems[i].iActive = -1;	
			}
			historyItems[iLastIndex].iPrev = -1;
			historyItems[iLastIndex].iActive = -1;
			
			delete [] m_historyItems;
			m_historyItems = historyItems;
			m_iHistoryItemAvailable = m_iHistoryItems;
			m_iHistoryItems *= 2;
		}
		// (2') there are inactive items: create a linked list with them
		else {
			int *iHistoryItemAux = &m_iHistoryItemAvailable;
			for(unsigned int i = 0 ; i < m_iHistoryItems ; ++i) {
				if (m_historyItems[i].iActive != m_iTimeCurrent) {
					m_historyItems[i].iActive = -1;
					m_historyItems[i].iEndFrame = -1;
					*iHistoryItemAux = i;
					iHistoryItemAux = &m_historyItems[i].iPrev;	
				}
			}
			*iHistoryItemAux = -1;
		}
	}
	
	// word-graph token garbage collection
	if (bRecycleWGTokens) {
		assert(iTokensActive <= (m_iWGTokens/m_iMaxWordSequencesState));
		if (iTokensActive >= 0.20*(m_iWGTokens/m_iMaxWordSequencesState)) {
			
			// allocate a new data structure with double capacity	
			WGToken *wgTokens = NULL;
			try {
				wgTokens = new WGToken[m_iWGTokens*2];
			}
			catch (const std::bad_alloc&) {
				int iBytes = m_iWGTokens*2*sizeof(WGToken);
				BVC_ERROR << "unable to allocate memory for the lattice tokens (" << iBytes << " Bytes needed)";
			}		
			//printf("WGTokens: from %d to %d\n",m_iWGTokens*sizeof(WGToken),m_iWGTokens*2*sizeof(WGToken));
			
			// copy the active items from the old data structure
			for(unsigned int i=0 ; i < m_iWGTokens ; ++i) {
				wgTokens[i].iWordSequence = m_wgTokens[i].iWordSequence;
				wgTokens[i].iLexUnitPron = m_wgTokens[i].iLexUnitPron;
				wgTokens[i].fScore = m_wgTokens[i].fScore;
				wgTokens[i].iHistoryItem = m_wgTokens[i].iHistoryItem;
				wgTokens[i].iPrev = -1;
				wgTokens[i].iActive = m_iTimeCurrent;
			}
			// create the linked list of available items
			int iLastIndex = ((2*m_iWGTokens)-m_iMaxWordSequencesState);
			for(int i=m_iWGTokens ; i < iLastIndex ; i += m_iMaxWordSequencesState) {
				wgTokens[i].iPrev = i+m_iMaxWordSequencesState;
				wgTokens[i].iActive = -1;	
			}
			wgTokens[iLastIndex].iPrev = -1;
			wgTokens[iLastIndex].iActive = -1;
			
			// swap structures and delete the old one
			delete [] m_wgTokens;
			m_wgTokens = wgTokens;
			m_iWGTokenAvailable = m_iWGTokens;
			m_iWGTokens *= 2;
		}
		// there are inactive tokens: create a linked list with them
		else {
			//printf("WGTokens: relinking\n");
			// the linked list goes from lower to higher memory addresses
			int *iAux = &m_iWGTokenAvailable;
			for(unsigned int i = 0 ; i < m_iWGTokens ; i += m_iMaxWordSequencesState) {
				if (m_wgTokens[i].iActive != m_iTimeCurrent) {	
					*iAux = i;
					iAux = &m_wgTokens[i].iPrev;	
				}	
			}
			*iAux = -1;
		}
	}	
	
	m_iTimeGarbageCollectionLast = m_iTimeCurrent;
	
	//printf("%d used %d total\n",iItemsActive,m_iHistoryItems);	
}

// build a hypothesis lattice for the utterance
HypothesisLattice *DynamicDecoderX::getHypothesisLattice() {

	assert(m_bInitialized);
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
		
	float fScoreBest = -FLT_MAX;
	float fScoreToken;
	float fScoreLM1;
	int iLMState;
	assert(m_lHistoryItem.empty());
	
	// (1) get all the history items at terminal arcs while keeping the best history item
	// for each unique word-sequence	
	
	// for each active token
	for(int i=0 ; i < m_iTokensNext ; ++i) {
		Token *token = &m_tokensNext[i];
		// check if pruned by "pruneExtraTokens"
		if (token->iNode == -1) {	
			continue;
		}		
		DNode *node = m_nodes+token->iNode;
		// only tokens at word-end nodes
		if (node->bWordEnd == false) {
			continue;
		}
		// two possibilities: 
		// (a) monophones: the node goes to a word-arc before any hmm-node is seen 
		// (b) multi-phones: the node goes directly to an hmm-state
		VLexUnit vLexUnitDest;
		getDestinationMonophoneLexUnits(node,vLexUnitDest);
		// (a) monophones
		if (vLexUnitDest.empty() == false) {
			for(VLexUnit::iterator it = vLexUnitDest.begin() ; it != vLexUnitDest.end() ; ++it) {	
				fScoreLM1 = 0.0;
				iLMState = -1;
				if (m_lexiconManager->isStandard(*it)) {
					iLMState = m_lmFSM->updateLMState(token->iLMState,(*it)->iLexUnit,&fScoreLM1);
					fScoreLM1 *= m_fLMScalingFactor;
				} else {
					iLMState = token->iLMState;
				}
				float fScoreLM2 = m_lmFSM->toFinalState(iLMState)*m_fLMScalingFactor;
				fScoreToken = token->fScore + fScoreLM1 + fScoreLM2;
				// keep best final score or discard if necessary
				if (fScoreToken > fScoreBest) {
					fScoreBest = fScoreToken;	
				} else if (fScoreToken < (fScoreBest-m_fBeamWidthNodesWE)) {
					continue;
				}					
				// create a history item for the token
				int iWGTokenAux = newWGToken(token->iWGToken,fScoreLM1+fScoreLM2);	// add the lm-score
				int iHistoryItem = newHistoryItem();
				HistoryItem *historyItem = m_historyItems+iHistoryItem;
				historyItem->iLexUnitPron = (*it)->iLexUnitPron;
				historyItem->iEndFrame = m_iFeatureVectorsUtterance-1;
				historyItem->fScore = fScoreToken;
				historyItem->iPrev = token->iHistoryItem;
				historyItem->iActive = m_iTimeCurrent;
				historyItem->iWGToken = iWGTokenAux;
				m_lHistoryItem.push_back(historyItem-m_historyItems);
				// keep the best history-item for each word-sequence
				keepBestHistoryItem(historyItem-m_historyItems);
			}
		}
		// (b) multi-phones 
		else {
			fScoreLM1 = m_lmFSM->toFinalState(token->iLMState)*m_fLMScalingFactor;
			fScoreToken = token->fScore + fScoreLM1;
			// keep best final score or discard if necessary
			if (fScoreToken > fScoreBest) {
				fScoreBest = fScoreToken;	
			} else if (fScoreToken < (fScoreBest-m_fBeamWidthNodesWE)) {
				continue;
			}	
			// create a history item for the token	
			int iWGTokenAux = newWGToken(token->iWGToken,fScoreLM1);		// add the lm-score
			int iHistoryItem = newHistoryItem();
			HistoryItem *historyItem = m_historyItems+iHistoryItem;
			historyItem->iLexUnitPron = token->iLexUnitPron;
			historyItem->iEndFrame = m_iFeatureVectorsUtterance-1;
			historyItem->fScore = fScoreToken;
			historyItem->iPrev = token->iHistoryItem;
			historyItem->iActive = m_iTimeCurrent;
			historyItem->iWGToken = iWGTokenAux;
			m_lHistoryItem.push_back(historyItem-m_historyItems);
			// keep the best history-item for each word-sequence
			keepBestHistoryItem(historyItem-m_historyItems);
		}
	}
	
	// if the map is empty the lattice cannot be built
	if (m_mWSHistoryItem.empty()) {
		BVC_WARNING << "unable to create the lattice, no history items were collected at terminal nodes";
		return NULL;
	}
	
	// create the initial node 
	LNode *lnodeInitial = HypothesisLattice::newNode(-1);
	
	// create the final node 
	LNode *lnodeFinal = HypothesisLattice::newNode(m_iFeatureVectorsUtterance-1);
	
	MHistoryItemLNode mHistoryItemLNode;	
	int iNodes = 0;
	int iEdges = 0;
	
	// disable wg-token entries that do not keep the best score for their word sequence
	for(list<int>::iterator it = m_lHistoryItem.begin() ; it != m_lHistoryItem.end() ; ) {
		
		int iHistoryItem = *it;	
		HistoryItem *historyItem = (*it)+m_historyItems;
		bool bSurvive = false;
		
		// above threshold: invalidate redundant word-sequences with lower likelihood
		if (historyItem->fScore > (fScoreBest-m_fBeamWidthNodesWE)) {
			
			// multiple tokens at the last time frame may hold the same lexical unit and share the same history-items
			// in their wg-token, thus it is necessary to disable those with lower scores
			WGToken *wgToken = m_wgTokens+historyItem->iWGToken;
			// keep original backpointer for the history-item
			int iHistoryItemBackpointer = historyItem->iPrev;
			for(int i=0 ; (i < m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1) ; ++i) {
				// get the word-sequence by moving the pointer
				historyItem->iPrev = wgToken[i].iHistoryItem;
				int iWordSequence = hashWordSequence(historyItem);
				map<int,pair<float,int> >::iterator it = m_mWSHistoryItem.find(iWordSequence);
				assert(it != m_mWSHistoryItem.end());
				// invalidate it if necessary
				if (it->second.second != iHistoryItem) {
					wgToken[i].iWordSequence = INT_MIN; 
				} else {
					bSurvive = true;
				}
			}	
			// recover original backpointer
			historyItem->iPrev = iHistoryItemBackpointer;	
		}
	
		if (bSurvive) {
			assert(historyItem != NULL);
			assert(historyItem->iWGToken != -1);	
			//printf("%x connected to final node with lex unit: \n");
			mHistoryItemLNode.insert(MHistoryItemLNode::value_type(historyItem,lnodeFinal));
			++it;
		} else {
			it = m_lHistoryItem.erase(it);
		}
	}
	m_mWSHistoryItem.clear();
	
	// process all the history items
	while(m_lHistoryItem.empty() == false) {
		
		HistoryItem *historyItemAux = m_lHistoryItem.back()+m_historyItems;
		m_lHistoryItem.pop_back();
		
		// get the graph node for this history item
		MHistoryItemLNode::iterator it = mHistoryItemLNode.find(historyItemAux);
		assert(it != mHistoryItemLNode.end());
		LNode *lnodeAux = it->second;	
		
		// items that go to this item
		map<int,bool> mItemSeen;
		for(int i=0 ; ((i < m_iMaxWordSequencesState) && ((historyItemAux->iWGToken+m_wgTokens)[i].iWordSequence != -1)) ; ++i) {
			WGToken &wgToken = (historyItemAux->iWGToken+m_wgTokens)[i];
			
			if (wgToken.iWordSequence == INT_MIN) {
				continue;
			}
			
			// get lexical unit
			LexUnit *lexUnit = m_lexiconManager->getLexUnitPron(wgToken.iLexUnitPron);
			if (lexUnit->iLexUnitPron == m_iLexUnitPronUnknown) {
				assert(lnodeAux == lnodeFinal);
				lexUnit = m_lexiconManager->getLexUnitPron(historyItemAux->iLexUnitPron);
			}
			
			// get the history item
			HistoryItem *historyItemPrev = m_historyItems+wgToken.iHistoryItem;	
			
			// connect to initial node?
			if (wgToken.iHistoryItem == m_iHistoryItemBegSentence) {
				LEdge *ledge = HypothesisLattice::newEdge(0,historyItemAux->iEndFrame,lexUnit,0.0,0.0,0.0);
				HypothesisLattice::connectEdge(lnodeInitial,ledge,lnodeAux);
				break;
			}
			
			assert(mItemSeen.find(wgToken.iHistoryItem) == mItemSeen.end());
			mItemSeen.insert(map<int,bool>::value_type(wgToken.iHistoryItem,true));
			LNode *lnodePrev = NULL;
			LEdge *ledgePrev = NULL;
			MHistoryItemLNode::iterator jt = mHistoryItemLNode.find(historyItemPrev);
			// the history item is not in the graph: create a graph node for it
			if (jt == mHistoryItemLNode.end()) {
				lnodePrev = HypothesisLattice::newNode(historyItemPrev->iEndFrame);
				ledgePrev = HypothesisLattice::newEdge(historyItemPrev->iEndFrame+1,historyItemAux->iEndFrame,lexUnit,0.0,0.0,0.0);
				mHistoryItemLNode.insert(MHistoryItemLNode::value_type(historyItemPrev,lnodePrev));
				m_lHistoryItem.push_back(historyItemPrev-m_historyItems);
				++iNodes;
			}
			// the history item is in the graph: create a link
			else {
				ledgePrev = HypothesisLattice::newEdge(historyItemPrev->iEndFrame+1,historyItemAux->iEndFrame,lexUnit,0.0,0.0,0.0);	
				lnodePrev = jt->second;
			}
			++iEdges;
			// make the connection
			HypothesisLattice::connectEdge(lnodePrev,ledgePrev,lnodeAux);
		}
 	}
 	
 	assert(lnodeInitial->edgeNext != NULL);
 	
 	HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager);
 	hypothesisLattice->buildContainer(lnodeInitial,lnodeFinal);
 	//hypothesisLattice->attachLMProbabilities(m_lmManager);
 	
 	double dTimeEnd = TimeUtils::getTimeMilliseconds();
 	double dTime = (dTimeEnd-dTimeBegin)/1000.0;
 	
 	BVC_VERB << "Lattice building time: " << FLT(8,4) << dTime << "s";

	return hypothesisLattice;
}

// keeps the best history item for each unique word-sequence (auxiliar method)
void DynamicDecoderX::keepBestHistoryItem(int iHistoryItem) {
		
	// multiple tokens at the last time frame may hold the same lexical unit and share the same history-items
	// in their wg-token, thus it is necessary to disable those with lower scores
	WGToken *wgToken = m_wgTokens+(m_historyItems+iHistoryItem)->iWGToken;
	// keep original backpointer for the history-item
	int iHistoryItemBackpointer = (m_historyItems+iHistoryItem)->iPrev;
	for(int i=0 ; (i < m_iMaxWordSequencesState) && (wgToken[i].iWordSequence != -1) ; ++i) {
		// get the word-sequence by moving the pointer
		(m_historyItems+iHistoryItem)->iPrev = wgToken[i].iHistoryItem;
		int iWordSequence = hashWordSequence(m_historyItems+iHistoryItem);
		map<int,pair<float,int> >::iterator it = m_mWSHistoryItem.find(iWordSequence);
		if (it == m_mWSHistoryItem.end()) {
			m_mWSHistoryItem.insert(map<int,pair<float,int> >::value_type(iWordSequence,
				pair<float,int>(wgToken[i].fScore,iHistoryItem)));	
		} else if (it->second.first < wgToken[i].fScore) {
			it->second.first = wgToken[i].fScore;
			it->second.second = iHistoryItem;
		}
	}	
	// recover original backpointer
	(m_historyItems+iHistoryItem)->iPrev = iHistoryItemBackpointer;
}

// merge two sets of word sequences by keeping the N best unique word sequences in wgToken1 (not commutative)
// 1) sorting: both sets are sorted so it is very efficient: O(n)
// 2) unique: this is linear too
bool DynamicDecoderX::mergeWordSequences(int iWGToken1, int iWGToken2) {

	WGToken *wgTokenTable = NULL;
	WGToken *wgToken1 = iWGToken1+m_wgTokens;
	WGToken *wgToken2 = iWGToken2+m_wgTokens;
	
	assert(iWGToken1 != iWGToken2);

	//assert(wgToken1[0].iWordSequence < 1000000);
	
	int iLength1 = 0;
	for( ; ((iLength1 < m_iMaxWordSequencesState) && (wgToken1[iLength1].iWordSequence != -1)) ; ++iLength1);
	int iLength2 = 0;
	for( ; ((iLength2 < m_iMaxWordSequencesState) && (wgToken2[iLength2].iWordSequence != -1)) ; ++iLength2);
	
	/*WGToken wgTokenAux1[m_iMaxWordSequencesState];
	WGToken wgTokenAux2[m_iMaxWordSequencesState];
	memcpy(&wgTokenAux1,wgToken1,sizeof(WGToken)*m_iMaxWordSequencesState);
	memcpy(&wgTokenAux2,wgToken2,sizeof(WGToken)*m_iMaxWordSequencesState);*/
	
	assert((iLength1 > 0) && (iLength2 > 0));
		
	int k = iLength1-1;	
	int j = iLength2-1;
	int iTotal = iLength1+iLength2;
	wgTokenTable = new WGToken[iTotal];
	
	// sort
	bool bReturn = false;
	for(int i=iTotal-1 ; i >= 0 ; --i) {
		if (wgToken1[k].fScore < wgToken2[j].fScore) {
			memcpy(&wgTokenTable[i],&wgToken1[k],sizeof(WGToken));
			--k;
			if (k < 0) {
				for(int h=i-1 ; h >= 0 ; --h) {
					assert(j>=0);
					memcpy(&wgTokenTable[h],&wgToken2[j--],sizeof(WGToken));
				}	
				bReturn = false;
				break;
			}
		} else {
			memcpy(&wgTokenTable[i],&wgToken2[j],sizeof(WGToken));
			--j;
			if (j < 0) {
				for(int h=i-1 ; h >= 0 ; --h) {
					assert(k>=0);
					memcpy(&wgTokenTable[h],&wgToken1[k--],sizeof(WGToken));
				}
				bReturn = true;
				break;
			}
		}	
	}
	
	// unique 
	int iUniqueElements = 0;
	memcpy(&wgToken1[0],&wgTokenTable[0],sizeof(WGToken));
	for(int i=1 ; ((i < iTotal) && (iUniqueElements+1 < m_iMaxWordSequencesState)) ; ++i) {
		bool bDuplicated = false;
		for(int j=0 ; j < iUniqueElements+1 ; ++j) {
			if (wgTokenTable[i].iWordSequence == wgToken1[j].iWordSequence) {
				bDuplicated = true;
				break;
			}
		}
		if (bDuplicated == false) {	
			memcpy(&wgToken1[++iUniqueElements],&wgTokenTable[i],sizeof(WGToken));
		}
	}	
	
	if (iUniqueElements+1 < m_iMaxWordSequencesState) {
		wgToken1[iUniqueElements+1].iWordSequence = -1;
	}
	
	//printWGToken(wgTokenTable);
	delete [] wgTokenTable;
	
	// sanity check
	/*int iElements = 0;
	for(int i=0 ; ((i < m_iMaxWordSequencesState) && (wgToken1[i].iWordSequence != -1)) ; ++i, ++iElements) {
		//printf("%3d -> %12.4f\n",i,wgToken1[i].fScore);
		if (i > 0) {
			//if (wgToken1[i-1].fScore <= wgToken1[i].fScore) {
			//	printWGToken(wgTokenAux1);
			//	printWGToken(wgTokenAux2);
			//	printWGToken(wgToken1);
			//}
			assert(wgToken1[i-1].fScore >= wgToken1[i].fScore);
		}
	}*/
	//assert(iElements == std::min(iLength1+iLength2,m_iMaxWordSequencesState));

	return bReturn;
}

// compute the load factor of the hash table containing unque word sequences
float DynamicDecoderX::computeLoadFactorHashWordSequences(int *iBucketsUsed, int *iCollisions) {

	// count 
	*iBucketsUsed = 0;
	*iCollisions = 0;
	for(unsigned int i=0 ; i < m_iWSHashEntries ; ++i) {
		if (m_wshashEntries[i].iTime > -1) {
			if (i < m_iWSHashBuckets) {
				++(*iBucketsUsed);
			} else {
				++(*iCollisions);
			}
		}
	}

	// load factor
	return ((float)(*iBucketsUsed))/((float)m_iWSHashBuckets);
}

// show hash-occupation information (debugging)
void DynamicDecoderX::printHashsStats() {

	int iBucketsUsed = 0;
	int iCollisions = 0;
	float fLoadFactor = computeLoadFactorHashWordSequences(&iBucketsUsed,&iCollisions);
	
	BVC_VERB << "- hash table containing word-sequences ---";
	BVC_VERB << " # buckets:      " << setw(8) << m_iWSHashBuckets;
	BVC_VERB << " # entries:      " << setw(8) << m_iWSHashEntries;
	BVC_VERB << " # buckets used: " << setw(8) << iBucketsUsed << " (" << FLT(5,2) << 
		100.0*((float)iBucketsUsed)/((float)m_iWSHashBuckets) << "%)";
	BVC_VERB << " # collisions:   " << setw(8) << iCollisions << " (" << FLT(5,2) <<
		100.0*((float)iCollisions)/((float)(iBucketsUsed+iCollisions)) << "%)";
	BVC_VERB << " # elements:     " << setw(8) << iBucketsUsed+iCollisions;
	BVC_VERB << " load factor:    " << FLT(8,4) << fLoadFactor;
	BVC_VERB << "------------------------------------------";
}

// print the hash-contents (debugging)
void DynamicDecoderX::printHashContents() {

	for(unsigned int i=0 ; i < m_iWSHashEntries ; ++i) {
		WSHashEntry *entry = &m_wshashEntries[i];
		if (entry->iTime != -1) {
			BVC_VERB << "[BUCKET]";
			do {
				BVC_VERB << "entry:";
				for(int j=0 ; j < entry->iLexUnits ; ++j) {
					m_lexiconManager->print(m_lexiconManager->getLexUnit(entry->iLexUnit[j]));
				}
				if (entry->iNext == -1) {
					break;
				}
				entry = m_wshashEntries+entry->iNext;
			} while(1);
		}
	}
}

// prune active tokens that are not in the top-N within a node 
void DynamicDecoderX::pruneExtraTokens(DNode *node) {

	assert(node->iActiveTokensNext >= m_iTokensNodeMax);

	// get the best score within the node
	float fScoreBestNode = -FLT_MAX;
	for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
		Token *token = m_tokensNext+(m_activeTokenNext+node->iActiveTokensNextBase)[j].iToken;
		if (token->fScore > fScoreBestNode) {
			fScoreBestNode = token->fScore;
		}
	}
		
	// compute bin-size for histogram pruning
	int iNumberBins = NUMBER_BINS_HISTOGRAM_WITHIN_NODE;
	int *iBins = new int[iNumberBins];
	int iBin;
	float fLength = m_fBeamWidthTokensNode+1;
	float fBinSize = fLength/((float)iNumberBins);
	assert(fBinSize > 0);
			
	float fThresholdLikelihood = fScoreBestNode-m_fBeamWidthTokensNode;

	// compute the size of each bin and initialize them
	for(int j = 0 ; j < iNumberBins ; ++j) {
		iBins[j] = 0;
	}
	// fill the bins, the first bin keeps the best tokens
	ActiveToken *activeTokensNext = m_activeTokenNext+node->iActiveTokensNextBase;
	for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
		Token *token = m_tokensNext+activeTokensNext[j].iToken;
		if (token->fScore >= fThresholdLikelihood) {
			iBin = (int)(fabs(token->fScore-fScoreBestNode)/fBinSize);
			assert((iBin >= 0) && (iBin < iNumberBins));
			iBins[iBin]++;
		}
	}
	// get the threshold
	int iSurvivors = 0;
	float fThresholdHistogram = -FLT_MAX;
	for(int j = 0 ; j < iNumberBins ; ++j) {
		int iSum = iSurvivors+iBins[j];
		// this is the cut-off
		if (iSum > m_iMaxActiveTokensNode) {	
			fThresholdHistogram = fScoreBestNode-(j*fBinSize);
			break;
		}
		iSurvivors = iSum;
	}
	float fThresholdNode = max(fThresholdLikelihood,fThresholdHistogram);
	
	// actual pruning
	int iPruned = 0;
	int iAvailable = -1;
	for(int j=0 ; j < node->iActiveTokensNext ; ++j) {
		Token *token = m_tokensNext+activeTokensNext[j].iToken;
		if (token->fScore < fThresholdNode) {
			if (iAvailable == -1) {
				iAvailable = j;	
			}
			token->iNode = -1;
			++iPruned;	
		} else if (iAvailable != -1) {
			activeTokensNext[iAvailable] = activeTokensNext[j];
			++iAvailable;
		}
	}
	node->iActiveTokensNext = node->iActiveTokensNext-iPruned;
	
	// this should never happen
	if ((node->iActiveTokensNext >= m_iTokensNodeMax) || (node->iActiveTokensNext == 0)) {
		for(int i=0 ; i < iNumberBins ; ++i) {
			BVC_VERB << setw(3) << i << " -> " << setw(4) << iBins[i];
		}
	}
	
	assert(node->iActiveTokensNext < m_iTokensNodeMax);

	delete [] iBins;
}

// return the active lm-states at the current time (lm-state in active tokens)
void DynamicDecoderX::getActiveLMStates(map<int,bool> &mLMState) {

	assert(m_bInitialized);

	// (1) active arcs for current time frame
	for(int i=0 ; i < m_iNodesActiveCurrent ; ++i) {
		ActiveToken *activeTokens = m_activeTokenCurrent+m_nodesActiveCurrent[i]->iActiveTokensCurrentBase;
		for(int j=0 ; j < m_nodesActiveCurrent[i]->iActiveTokensCurrent ; ++j) {
			Token *token = m_tokensCurrent+activeTokens[j].iToken;
			if (token->iLANode != -1) {
				assert(token->iLMState != -1);
				mLMState[token->iLMState] = true;
			}
		}
	}	
	// (2) active arcs for next time frame
	for(int i=0 ; i < m_iNodesActiveNext ; ++i) {
		ActiveToken *activeTokens = m_activeTokenNext+m_nodesActiveNext[i]->iActiveTokensNextBase;
		for(int j=0 ; j < m_nodesActiveNext[i]->iActiveTokensNext ; ++j) {
			Token *token = m_tokensNext+activeTokens[j].iToken;
			if (token->iLANode != -1) {
				assert(token->iLMState != -1);
				mLMState[token->iLMState] = true;
			}
		}
	}
}

};	// end-of-namespace

