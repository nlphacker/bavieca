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


#include "NetworkBuilderX.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

#include <iomanip>

namespace Bavieca {

// constructor
NetworkBuilderX::NetworkBuilderX(PhoneSet *phoneSet, HMMManager *hmmManager, LexiconManager *lexiconManager)
{
	m_lexiconManager = lexiconManager;
	m_hmmManager = hmmManager;
	m_phoneSet = phoneSet;
	
	m_iNodesTemp = 0;
	m_iNodesTempRoot = 0;
	m_iNodesTempFI = 0;	
	m_iNodesTempMI = 0;
	m_iNodesTempFO = 0;	
	m_iNodesTempWord = 0;
	m_iNodesTempHMM = 0;
	m_iNodesTempNull = 0;
	m_iNodesDepth = 0;	
	
	m_iNodesTempID = 0;
		
	m_iArcsTemp = 0;
	m_iArcsNull = 0;
	m_iArcsHMM = 0;
	m_iArcsWord = 0;
   	
	m_nodeTempRoot = NULL;
	
	// language model look-ahead
	m_iLANodes = -1;
	m_iLATree = NULL;	
}

// destructor
NetworkBuilderX::~NetworkBuilderX()
{	
}

// build the network
// - the goal is to minimize the number of hmm-nodes in the network, decoding speed is directly proportional
//   to the number of active states, so we do not want redundant active states
DynamicNetworkX *NetworkBuilderX::build() {
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
   
   // get the context size
   int iContextSizeWW = m_hmmManager->getContextSizeHMM();			// within-word
   int iContextSizeCW = m_hmmManager->getContextSizeHMMCW();		// cross-word
   int iBasePhones = m_phoneSet->size();								// used for context padding
   
   assert(iContextSizeWW == iContextSizeCW);
   
   // this ensures that the network is built right for context-independent models
   if (iContextSizeWW == 0) {
	   iContextSizeWW = iContextSizeCW = 1;
   }   

   // get the lexicon 
   VLexUnit *lexicon = m_lexiconManager->getLexiconReference(); 
   
   // create the root node
   m_nodeTempRoot = newNodeTempX(NODE_TYPE_ROOT,0);
   
   // allocate memory for mapping phonemes to IP values
   float *fPhoneIP = new float[iBasePhones];
   for(int i=0 ; i < iBasePhones ; ++i) {
   	fPhoneIP[i] = FLT_MAX; 
   }
   
   // get the length of the longest lexical unit
   unsigned int iPhonesLongestLexicalUnit = 0;
   for (VLexUnit::iterator it = lexicon->begin() ; it != lexicon->end() ; ++it) {
		// keep the max number of phones
		if ((*it)->vPhones.size() > iPhonesLongestLexicalUnit) {
			iPhonesLongestLexicalUnit = (*it)->vPhones.size();
		}
   }

   // compute node-depths at different layers
   m_iDepthFI = 1;
   m_iDepthMI = 5;
   m_iDepthFO = m_iDepthMI + (iPhonesLongestLexicalUnit-2)*NUMBER_HMM_STATES + 1;
 
   // get all the possible left/right contexts (word ends/beginnings)
   MContextBool mContextLeft;
   MContextBool mContextRight;
   MContextNode mContextBeg;				// word beginnings
   MContextNode mContextEnd;				// word ends
   MContextV mContextBegV;					// word beginnings
   MContextV mContextEndV;					// word ends
   
   unsigned char *iContextLeft = new unsigned char[iContextSizeCW+1]; 
   unsigned char *iContextRight = new unsigned char[iContextSizeCW+1]; 
   unsigned char *iContextBeg = new unsigned char[iContextSizeWW+2]; 
   unsigned char *iContextBegPartial = new unsigned char[iContextSizeWW+1]; 
   unsigned char *iContextEnd = new unsigned char[iContextSizeWW+2];
   unsigned char *iContextEndPartial = new unsigned char[iContextSizeWW+1];
   for(VLexUnit::iterator it = lexicon->begin() ; it != lexicon->end() ; ++it) {
   	
   	// only standard and filler lexical units (no sentence markers)
   	if (((*it)->iType != LEX_UNIT_TYPE_STANDARD) && ((*it)->iType != LEX_UNIT_TYPE_FILLER)) {
   		continue;
   	}
   	
   	// keep the IP value
   	if (fPhoneIP[(*it)->vPhones.front()] == FLT_MAX) {
   		fPhoneIP[(*it)->vPhones.front()] = (*it)->fInsertionPenalty;
   	} else if (fPhoneIP[(*it)->vPhones.front()] != (*it)->fInsertionPenalty) {
   		BVC_ERROR << "two lexical units with different insertion penalties share their initial phone";
   	}
   	
   	int iPhones = (*it)->vPhones.size();
   	int iContextLength = min(iContextSizeCW,(int)(*it)->vPhones.size());
   	int iPhonesWW = min(iContextSizeWW+1,(int)(*it)->vPhones.size());
   	
   	// (1) get all the context
   	
   	// (1.1) left context
   	// context padding
   	for(int i=0 ; i < iContextSizeCW-iContextLength ; ++i) {
			iContextLeft[i] = iBasePhones; 
   	}
   	// actual context
   	for(int i=0 ; i < iContextLength ; ++i) {
			iContextLeft[iContextSizeCW-iContextLength+i] = (*it)->vPhones[(iPhones-iContextLength)+i]; 
   	}
   	iContextLeft[iContextSizeCW] = UCHAR_MAX;
   	if (mContextLeft.find(iContextLeft) == mContextLeft.end()) {
   		mContextLeft.insert(MContextBool::value_type(copyContext(iContextLeft,iContextSizeCW),true));
   	}
   	
   	// (1.2) right context
   	// context padding
   	for(int i=iContextLength ; i < iContextSizeCW ; ++i) {
			iContextRight[i] = iBasePhones; 
   	}
   	// actual context
   	for(int i=0 ; i < iContextLength ; ++i) {
			iContextRight[i] = (*it)->vPhones[i]; 
   	}
   	iContextRight[iContextSizeCW] = UCHAR_MAX;
   	if (mContextRight.find(iContextRight) == mContextRight.end()) {
   		mContextRight.insert(MContextBool::value_type(copyContext(iContextRight,iContextSizeCW),true));	
   	}
   	
   	// (1.3) word-ends (needed to create the FO-nodes)
		// context padding
		for(int i=0 ; i < iContextSizeWW+1-iPhonesWW ; ++i) {
			iContextEnd[i] = iBasePhones; 
		}
		// actual context
		for(int i=0 ; i < iPhonesWW ; ++i) {
			iContextEnd[iContextSizeWW+1-iPhonesWW+i] = (*it)->vPhones[(iPhones-iPhonesWW)+i]; 
		}
		iContextEnd[iContextSizeWW+1] = UCHAR_MAX;
		if (mContextEnd.find(iContextEnd) == mContextEnd.end()) {	
			memcpy(iContextEndPartial,iContextEnd+1,sizeof(unsigned char)*iContextSizeWW);
			iContextEndPartial[iContextSizeWW] = UCHAR_MAX;
			//print(iContextEndPartial);
			//print(iContextLeft);
			// no FO-node is needed for monophone lexical units
			NodeTempX *nodeFO = newNodeTempX(NODE_TYPE_FO,m_iDepthFO);
			MContextV::iterator jt = mContextEndV.find(iContextEndPartial);
			if (jt == mContextEndV.end()) {
				mContextEndV[copyContext(iContextEndPartial,iContextSizeWW)].push_back(pair<unsigned char*,NodeTempX*>
				(copyContext(iContextEnd,iContextSizeWW+1),nodeFO));
			} else {
				jt->second.push_back(pair<unsigned char*,NodeTempX*>(copyContext(iContextEnd,iContextSizeWW+1),nodeFO));
			}
			//mContextEndV[copyContext(iContextEndPartial,iContextSizeWW)].push_back(pair<unsigned char*,NodeTempX*>
			//	(copyContext(iContextEnd,iContextSizeWW+1),nodeFO));
			mContextEnd.insert(MContextNode::value_type(copyContext(iContextEnd,iContextSizeWW+1),nodeFO));
		}
		
		// (1.4) word-beginning (needed to create the MI-nodes)
		// context padding
		for(int i=iPhonesWW ; i < iContextSizeWW+1 ; ++i) {
			iContextBeg[i] = iBasePhones; 
		}
		// actual context
		for(int i=0 ; i < iPhonesWW ; ++i) {
			iContextBeg[i] = (*it)->vPhones[i]; 
		}
		iContextBeg[iContextSizeWW+1] = UCHAR_MAX;	
		if (mContextBeg.find(iContextBeg) == mContextBeg.end()) {
			memcpy(iContextBegPartial,iContextBeg,sizeof(unsigned char)*iContextSizeWW);
			iContextBegPartial[iContextSizeWW] = UCHAR_MAX;
			// no MI-node is needed for monophone lexical units
			NodeTempX *nodeMI = newNodeTempX(NODE_TYPE_MI,m_iDepthMI);
			MContextV::iterator jt = mContextBegV.find(iContextBegPartial);
			if (jt == mContextBegV.end()) {
				mContextBegV[copyContext(iContextBegPartial,iContextSizeWW)].push_back(pair<unsigned char*,NodeTempX*>
				(copyContext(iContextBeg,iContextSizeWW+1),nodeMI));
			} else {
				jt->second.push_back(pair<unsigned char*,NodeTempX*>(copyContext(iContextBeg,iContextSizeWW+1),nodeMI));
			}	
			mContextBeg.insert(MContextNode::value_type(copyContext(iContextBeg,iContextSizeWW+1),nodeMI));
		}
	}
	
	// clean-up
   delete [] iContextLeft; 
   delete [] iContextRight; 
   delete [] iContextBeg; 
   delete [] iContextBegPartial;
   delete [] iContextEnd;
   delete [] iContextEndPartial;	
	
	cout << "left contexts:     " << setw(9) << mContextLeft.size() << endl;
	cout << "right contexts:    " << setw(9) << mContextRight.size() << endl;
	cout << "word begs:         " << setw(9) << mContextBeg.size() << endl;
	cout << "word ends:         " << setw(9) << mContextEnd.size() << endl;
	cout << "word begs partial: " << setw(9) << mContextBegV.size() << endl;
	cout << "word ends partial: " << setw(9) << mContextEndV.size() << endl;
	   
   map<int,bool> mHMMUsed;
   
	double dTimeBeginFI = TimeUtils::getTimeMilliseconds(); 
	
	// (1) HEAD-NETWORK
	// - use auxiliar FI-nodes and MI-nodes, they are useful for minimization. They will be removed
	//   when connecting the head-network to the other subnetworks 
	
   MContextNode mNodesFI;
   VNodeTempX vNodesFIUnique;
		
	// allocate data structures
	NodeMergeInfo *fiNodeMergeInfo = new NodeMergeInfo[mContextLeft.size()];
	unsigned char *iPhoneLeft = new unsigned char[iContextSizeCW];
	unsigned char *iPhoneRight = new unsigned char[iContextSizeWW];	
	
	// create the head network
	for(MContextV::iterator it = mContextBegV.begin() ; it != mContextBegV.end() ; ++it) {
	
		// keep the MI-nodes
		VNodeTempX vNodeTempMI;
		for(vector<pair<unsigned char*,NodeTempX*> >::iterator jt = it->second.begin() ; jt != it->second.end() ; ++jt) {
			vNodeTempMI.push_back(jt->second);
		}
		
		int iFI = 0;
		for(MContextBool::iterator jt = mContextLeft.begin() ; jt != mContextLeft.end() ; ++jt, ++iFI) {

			// create the FI-node
			unsigned char *iKey = new unsigned char[2*iContextSizeCW+1];
			for(int i=0 ; i<iContextSizeCW ; ++i) {
				iKey[i] = jt->first[i];
				iKey[iContextSizeCW+i] = it->first[i];
			}
			iKey[2*iContextSizeCW] = UCHAR_MAX;
			NodeTempX *nodeFI = newNodeTempX(NODE_TYPE_FI,m_iDepthFI);
			nodeFI->iFI = iFI;
				
			fiNodeMergeInfo[iFI].iMerged = -1;
			fiNodeMergeInfo[iFI].node = nodeFI;	
			fiNodeMergeInfo[iFI].iKey = iKey;

			for(int i=0 ; i < iContextSizeCW ; ++i) {
				iPhoneLeft[i] = jt->first[i];
			}	
	
			// MI-nodes
			for(vector<pair<unsigned char*,NodeTempX*> >::iterator kt = it->second.begin() ; kt != it->second.end() ; ++kt) {
				
				NodeTempX *nodeMI = kt->second;
				
				// ignore MI-node for monophone lexical units
				if (kt->first[1] == iBasePhones) {	
					continue;
				}	
				
				unsigned char iPhone = kt->first[0];
				for(int i=0 ; i < iContextSizeWW ; ++i) {
					iPhoneRight[i] = kt->first[1+i];
				}
					
				// get the sequence of states for the n-phone
				NodeTempX *nodePrev = nodeFI;
				for(int iState=0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				
					// get the HMM-state
					HMMStateDecoding *state = m_hmmManager->getHMMStateDecoding(iPhoneLeft,iPhone,iPhoneRight,
						WITHIN_WORD_POSITION_START,iState);
					mHMMUsed[state->getId()] = true;
				
					// insert a node per HMM-state if they do not exist yet
					NodeTempX *nodeAux = NULL;
					for(LArcTempX::iterator kt = nodePrev->lArcNext.begin() ; kt != nodePrev->lArcNext.end() ; ++kt) {
						if ((*kt)->nodeDest->state == state) {
							nodeAux = (*kt)->nodeDest;
							break;
						}
					}	
					
					char iIPIndex = -1;
					if (iState == 0) {
						assert(fPhoneIP[iPhone] != FLT_MAX);
						iIPIndex = iPhone;
					}
					
					if (nodeAux == NULL) {
						nodeAux = newNodeTempX(NODE_TYPE_HMM,m_iDepthFI+iState+1,state,iIPIndex);
						newArcTempX(nodePrev,nodeAux,ARC_TYPE_NULL,NULL,NULL);
					}
					
					// connect to the MI-node
					if (iState == NUMBER_HMM_STATES-1) {
						newArcTempX(nodeAux,nodeMI,ARC_TYPE_NULL,NULL,NULL);	
					}
					nodePrev = nodeAux;
				}
			}	
		}
		// minimize backwards from the middle node
		mergeBackwardMI(vNodeTempMI,fiNodeMergeInfo);
		// keep surviving FI-nodes
		for(unsigned int i=0 ; i < mContextLeft.size() ; ++i) {
			int iMerged = fiNodeMergeInfo[i].iMerged;
			// (a) FI-node was not merged: keep it
			if (iMerged == -1) {
				mNodesFI.insert(MContextNode::value_type(fiNodeMergeInfo[i].iKey,fiNodeMergeInfo[i].node));
				vNodesFIUnique.push_back(fiNodeMergeInfo[i].node);
			} 
			// (b) FI-node was merged: point to the one that we will keep
			else {
				iMerged = i;
				while(fiNodeMergeInfo[iMerged].iMerged != -1) {
					iMerged = fiNodeMergeInfo[iMerged].iMerged;
				}
				mNodesFI.insert(MContextNode::value_type(fiNodeMergeInfo[i].iKey,fiNodeMergeInfo[iMerged].node));
			}
		}
   }
   
   delete [] iPhoneRight;
   delete [] iPhoneLeft; 
   delete [] fiNodeMergeInfo;
	 
	double dTimeEndFI = TimeUtils::getTimeMilliseconds();
	
	cout << "nodes FI: " << mNodesFI.size() << endl;
   
   printf("FI-done: %f12.4\n",(dTimeEndFI-dTimeBeginFI)/1000.0);
   
	print();
	
	//exit(-1);
	
	// (2) TAIL-NETWORK
	// - use auxiliar FI-nodes and MI-nodes, they are useful for minimization. They will be removed
	//   when connecting the head-network to the other subnetworks 
   
   int iLeft = mContextEndV.size();
   
   NodeMergeInfo *foNodeMergeInfo = new NodeMergeInfo[m_iNodesTempFO];
   int iIndex=0;
   for(MContextNode::iterator it = mContextEnd.begin() ; it != mContextEnd.end() ; ++it, ++iIndex) {
		foNodeMergeInfo[iIndex].node = it->second;
		it->second->iFI = iIndex;
		foNodeMergeInfo[iIndex].iMerged = -1;
		foNodeMergeInfo[iIndex].iKey = it->first;
   }
   assert(iIndex == m_iNodesTempFO);
	
	// create the head network
	for(MContextV::iterator jt = mContextEndV.begin() ; jt != mContextEndV.end() ; ++jt, --iLeft) {	
	
		for(MContextBool::iterator it = mContextRight.begin() ; it != mContextRight.end() ; ++it) {
				
			// get the FI-node
			unsigned char *iKey = new unsigned char[2*iContextSizeCW+1];
			for(int i=0 ; i<iContextSizeCW ; ++i) {
				iKey[i] = jt->first[i];
				iKey[iContextSizeCW+i] = it->first[i];
			}
			iKey[2*iContextSizeCW] = UCHAR_MAX;
			NodeTempX *nodeFI = mNodesFI[iKey];
			if (nodeFI == NULL) {
				print(iKey);
			}
			delete [] iKey;
			
			assert(nodeFI != NULL);
			
			unsigned char *iPhoneRight = new unsigned char[iContextSizeCW];
			for(int i=0 ; i < iContextSizeCW ; ++i) {
				iPhoneRight[i] = it->first[i];
			}	
				
			// FO-nodes
			for(vector<pair<unsigned char*,NodeTempX*> >::iterator kt = jt->second.begin() ; kt != jt->second.end() ; ++kt) {
				
				NodeTempX *nodeFO = kt->second;
				
				// ignore FO-node for monophone lexical units
				if (kt->first[iContextSizeCW-1] == iBasePhones) {
					continue;
				}
					
				unsigned char iPhone = kt->first[iContextSizeWW];
				assert(iPhone != iBasePhones);
				unsigned char *iPhoneLeft = new unsigned char[iContextSizeWW];
				for(int i=0 ; i < iContextSizeWW ; ++i) {
					iPhoneLeft[i] = kt->first[i];
				}
					
				// get the sequence of states for the n-phone
				NodeTempX *nodePrev = nodeFO;
				for(int iState=0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				
					// get the HMM-state
					HMMStateDecoding *state = m_hmmManager->getHMMStateDecoding(iPhoneLeft,iPhone,iPhoneRight,
						WITHIN_WORD_POSITION_END,iState);
					mHMMUsed[state->getId()] = true;
				
					// insert a node per HMM-state if they do not exist yet
					NodeTempX *nodeAux = NULL;
					for(LArcTempX::iterator kt = nodePrev->lArcNext.begin() ; kt != nodePrev->lArcNext.end() ; ++kt) {
						if ((*kt)->nodeDest->state == state) {
							nodeAux = (*kt)->nodeDest;
							break;
						}
					}	
					
					if (nodeAux == NULL) {
						nodeAux = newNodeTempX(NODE_TYPE_HMM,m_iDepthFO+iState+1,state);
						if (iState == NUMBER_HMM_STATES-1) {
							nodeAux->bWordEnd = true;
						}
						newArcTempX(nodePrev,nodeAux,ARC_TYPE_NULL,NULL,NULL);
					}
					// connect to the FI-node
					if (iState == NUMBER_HMM_STATES-1) {						
						newArcTempX(nodeAux,nodeFI,ARC_TYPE_NULL,NULL,NULL);
					}
					nodePrev = nodeAux;
				}
				delete [] iPhoneLeft;
			}
			delete [] iPhoneRight;
		}
   }
	
	// minimize backwards from the FI nodes
	printf("merging backward\n");
	mergeBackwardFI(vNodesFIUnique,foNodeMergeInfo);
	
	// keep surviving FO-nodes
	for(unsigned int i=0 ; i < mContextEnd.size() ; ++i) {
		// (b) FO-node was merged: point to the one that will be kept
		if (foNodeMergeInfo[i].iMerged != -1) {
			int iMerged = foNodeMergeInfo[i].iMerged;
			while(foNodeMergeInfo[iMerged].iMerged != -1) {
				iMerged = foNodeMergeInfo[iMerged].iMerged;
			}
			mContextEnd[foNodeMergeInfo[i].iKey] = foNodeMergeInfo[iMerged].node;
		}
	}	
	delete [] foNodeMergeInfo;
	
	print();
	
	//exit(-1);
	
	for(VNodeTempX::iterator it = vNodesFIUnique.begin() ; it != vNodesFIUnique.end() ; ++it) {
		newArcTempX(m_nodeTempRoot,*it,ARC_TYPE_NULL,NULL,NULL);
	}
	
	double dTimeEndFO = TimeUtils::getTimeMilliseconds();
   
   printf("FO-done: %f12.4\n",(dTimeEndFO-dTimeEndFI)/1000.0);
   
	print();
	
	//exit(-1);
	
	VLexUnit *vLexUnitMonophone = new VLexUnit[iBasePhones];
   
  	// (3) INTERNAL-NETWORK (insert one lexical unit at a time)
   for (VLexUnit::iterator it = lexicon->begin() ; it != lexicon->end() ; ++it) {
   	
   	// only standard and filler lexical units (no sentence markers)
   	if (((*it)->iType != LEX_UNIT_TYPE_STANDARD) && ((*it)->iType != LEX_UNIT_TYPE_FILLER)) {
   		continue;
   	}
   	
   	//m_lexiconManager->print(*it);
  	
   	// get the number of phones in the word
   	int iPhones = (*it)->vPhones.size();
   	
   	// (7.1) monophone lexical units
   	if (iPhones == 1) {
   	
   		// keep the monophone lexical unit to be processed later
   		vLexUnitMonophone[(*it)->vPhones.front()].push_back(*it); 	
   	}	
   	// (7.2) multiple-phone lexical units
   	else if (iPhones > 1) {
   	
			int iPhonesWW = min(iContextSizeWW+1,(int)(*it)->vPhones.size());
	
			// (1) get the MI-node
			unsigned char *iKeyMI = new unsigned char[iContextSizeWW+2];
			
			// padding
			for(int i=iPhonesWW ; i < iContextSizeWW+1 ; ++i) {
				iKeyMI[i] = iBasePhones; 
			}
			// actual context
			for(int i=0 ; i < iPhonesWW ; ++i) {
				iKeyMI[i] = (*it)->vPhones[i]; 
			}
			iKeyMI[iContextSizeWW+1] = UCHAR_MAX;
			
			// (2) get the FO-node
			unsigned char *iKeyFO = new unsigned char[iContextSizeWW+2];
		
			// padding
			for(int i=0 ; i < iContextSizeWW+1-iPhonesWW ; ++i) {
				iKeyFO[i] = iBasePhones; 
			}
			// actual context
			for(int i=0 ; i < iPhonesWW ; ++i) {
				iKeyFO[iContextSizeWW+1-iPhonesWW+i] = (*it)->vPhones[iPhones-(iPhonesWW-i)]; 
			}
			iKeyFO[iContextSizeWW+1] = UCHAR_MAX;
			
			// get the MI-node
			NodeTempX *nodeMI = mContextBeg[iKeyMI];	
			// get the FO-node
			NodeTempX *nodeFO = mContextEnd[iKeyFO];	
			
			//print(iKeyMI);
			//print(iKeyFO);
			
			delete [] iKeyMI;
			delete [] iKeyFO;
			
			assert(nodeMI != NULL);
			assert(nodeFO != NULL);
			
			NodeTempX *nodePrev = nodeMI;
			for(int i=1 ; i < iPhones-1 ; ++i) {
			
				unsigned char iPhone = (*it)->vPhones[i];
				unsigned char *iPhoneLeft = new unsigned char[iContextSizeWW];
				unsigned char *iPhoneRight = new unsigned char[iContextSizeWW];
	
				// left context
				unsigned char iPhonesAvailableLeft = min(iContextSizeWW,i);
				// padding
				for(int j=0 ; j < iContextSizeWW-i ; ++j) {
					iPhoneLeft[j] = iBasePhones;
				}
				// context
				for(int j=0 ; j < iPhonesAvailableLeft ; ++j) {
					iPhoneLeft[iContextSizeWW-j-1] = (*it)->vPhones[i-j-1];
				}	
				
				// right context
				unsigned char iPhonesAvailableRight = min(iContextSizeWW,iPhones-i-1);
				// context 
				for(int j=0 ; j < iPhonesAvailableRight ; ++j) {
					iPhoneRight[j] = (*it)->vPhones[i+j+1];
				}
				// padding
				for(int j=iPhonesAvailableRight ; j < iContextSizeWW ; ++j) {
					iPhoneRight[j] = iBasePhones;
				}
				
				//print(iPhoneLeft,2);
				//print(iPhoneRight,2);
			
				//printf("%s-%s+%s\n",m_phoneSet->getStrPhone(iPhoneLeft[0]),m_phoneSet->getStrPhone(iPhone),m_phoneSet->getStrPhone(iPhoneRight[0]));
				
				// get the sequence of states for the n-phone
				for(int iState=0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				
					// get the HMM-state
					HMMStateDecoding *state = m_hmmManager->getHMMStateDecoding(iPhoneLeft,iPhone,iPhoneRight,
						WITHIN_WORD_POSITION_INTERNAL,iState);
					mHMMUsed[state->getId()] = true;
					
					// insert a node per HMM-state if they do not exist yet
					NodeTempX *nodeAux = NULL;
					for(LArcTempX::iterator kt = nodePrev->lArcNext.begin() ; kt != nodePrev->lArcNext.end() ; ++kt) {
						if ((*kt)->nodeDest->state == state) {
							nodeAux = (*kt)->nodeDest;
							break;
						}
					}	
					
					if (nodeAux == NULL) {
						nodeAux = newNodeTempX(NODE_TYPE_HMM,m_iDepthMI+1+((i-1)*NUMBER_HMM_STATES)+iState,state);
						/*if (nodePrev == nodeMI) {
							for(LArcTempX::iterator kt = nodeMI->lArcPrev.begin() ; kt != nodeMI->lArcPrev.end() ; ++kt) {
								newArcTempX((*kt)->nodeSource,nodeAux,ARC_TYPE_NULL,NULL,NULL,false,-1);
							}
						} else {*/
							newArcTempX(nodePrev,nodeAux,ARC_TYPE_NULL,NULL,NULL);
						//}
					} 
					nodePrev = nodeAux;
				}
				delete [] iPhoneLeft;
				delete [] iPhoneRight;
			}
			
			/*if (nodePrev == nodeMI) {
				for(LArcTempX::iterator kt = nodeMI->lArcPrev.begin() ; kt != nodeMI->lArcPrev.end() ; ++kt) {
					newArcTempX((*kt)->nodeSource,nodeFO,ARC_TYPE_WORD,NULL,*it,false,-1);
				}
			} else {*/
				// create the WI arc to the FO-node
				newArcTempX(nodePrev,nodeFO,ARC_TYPE_WORD,NULL,*it);
			//}
			
			/*for(LArcTempX::iterator kt = nodeFO->lArcNext.begin() ; kt != nodeFO->lArcNext.end() ; ++kt) {
				newArcTempX(nodePrev,(*kt)->nodeDest,ARC_TYPE_WORD,NULL,*it,false,-1);
			}*/
		}
	}
	
	// deal with monophone lexical units
	unsigned char *iKeyFI = new unsigned char[iContextSizeCW+iContextSizeWW+1];
	unsigned char *iKeyFIDest = new unsigned char[iContextSizeCW+iContextSizeWW+1];
	//unsigned char *iKeyMI = new unsigned char[iContextSizeWW+1];
	
	
	double dTimeBeginMonophones = TimeUtils::getTimeMilliseconds();
   
	// for each possible phone	
	for(int iPhone=0 ; iPhone < iBasePhones ; ++iPhone) {
		
		if (vLexUnitMonophone[iPhone].empty()) {
			continue;
		}
		
		printf("phone: %s\n",m_phoneSet->getStrPhone(iPhone));
		int iHMMNodes = 0;
		
		hash_map<int,NodeTempX*> mFIWordNode;
		VNodeTempX vNodeWord;
		
		// get all the possible destination FI-nodes
		vector< pair<unsigned char*,NodeTempX*> > vNodeFIDest;
		for(MContextBool::iterator kt = mContextRight.begin() ; kt != mContextRight.end() ; ++kt) {
		
			// get the destination FI-node
			for(int i=0 ; i < iContextSizeCW-1 ; ++i) {
				iKeyFIDest[i] = iBasePhones;
			}
			iKeyFIDest[iContextSizeCW-1] = iPhone;	
			for(int i=iContextSizeCW ; i < iContextSizeCW+iContextSizeWW ; ++i) {
				iKeyFIDest[i] = kt->first[i-iContextSizeCW];
			}	
			iKeyFIDest[iContextSizeCW+iContextSizeWW] = UCHAR_MAX;
			NodeTempX *nodeFIDest = mNodesFI[iKeyFIDest];	
			assert(nodeFIDest != NULL);
			vNodeFIDest.push_back(pair<unsigned char*,NodeTempX*>(copyContext(iKeyFIDest,2*iContextSizeCW),nodeFIDest));	
		}	
		
		// for each possible left cross-word context	
		for(MContextBool::iterator jt = mContextLeft.begin() ; jt != mContextLeft.end() ; ++jt) {	
			
			// get the FI-node
			for(int i=0 ; i < iContextSizeCW ; ++i) {
				iKeyFI[i] = jt->first[i];
			}
			iKeyFI[iContextSizeCW] = iPhone;	
			for(int i=iContextSizeCW+1 ; i < iContextSizeCW+iContextSizeWW ; ++i) {
				iKeyFI[i] = iBasePhones;
			}	
			iKeyFI[iContextSizeCW+iContextSizeWW] = UCHAR_MAX;
			NodeTempX *nodeFI = mNodesFI[iKeyFI];	
			assert(nodeFI != NULL);	
			
			unsigned char *iPhoneLeft = new unsigned char[iContextSizeWW];	
			unsigned char *iPhoneRight = new unsigned char[iContextSizeWW];	
			
			// left context	
			for(int i=0 ; i < iContextSizeCW ; ++i) {
				iPhoneLeft[i] = jt->first[i];
			}		
			
			// for each possible destination FI-node
			for(vector< pair<unsigned char*,NodeTempX*> >::iterator kt = vNodeFIDest.begin() ; 
				kt != vNodeFIDest.end() ; ++kt) {

				NodeTempX *nodeFIDest = kt->second;	
				assert(nodeFIDest != NULL);
				
				NodeTempX *nodeWord = NULL;
				hash_map<int,NodeTempX*>::iterator lt = mFIWordNode.find(nodeFIDest->iNode);
				if (lt == mFIWordNode.end()) {
					nodeWord = newNodeTempX(NODE_TYPE_WORD,m_iDepthFI+1+NUMBER_HMM_STATES,NULL);
					vNodeWord.push_back(nodeWord);
					NodeTempX *nodeWord2 = newNodeTempX(NODE_TYPE_NULL,m_iDepthFI+1+NUMBER_HMM_STATES,NULL);
					mFIWordNode[nodeFIDest->iNode] = nodeWord;
					//printf("fi-id: %d\n",nodeFIDest->iNode);
					for(VLexUnit::iterator it = vLexUnitMonophone[iPhone].begin(); 
						it != vLexUnitMonophone[iPhone].end() ; ++it) {
						//for(LArcTempX::iterator mt = nodeFIDest->lArcNext.begin() ; mt != nodeFIDest->lArcNext.end() ; ++mt) {
						//	newArcTempX(nodeWord,(*mt)->nodeDest,ARC_TYPE_WORD,NULL,*it,false,-1);
						//}
						//newArcTempX(nodeWord,nodeFIDest,ARC_TYPE_WORD,NULL,*it,false,-1);	
						newArcTempX(nodeWord,nodeWord2,ARC_TYPE_WORD,NULL,*it);	
					}
					newArcTempX(nodeWord2,nodeFIDest,ARC_TYPE_NULL,NULL,NULL);	
				} else {	
					nodeWord = lt->second;
				}
				
				// right context	
				for(int i=0 ; i < iContextSizeCW ; ++i) {
					iPhoneRight[i] = kt->first[iContextSizeCW+i];
				}
		
				// get the sequence of states for the n-phone
				NodeTempX *nodePrev = nodeFI;
				for(int iState=0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				
					// get the HMM-state
					HMMStateDecoding *state = m_hmmManager->getHMMStateDecoding(iPhoneLeft,iPhone,iPhoneRight,
						WITHIN_WORD_POSITION_MONOPHONE,iState);
					mHMMUsed[state->getId()] = true;	
						
					// insert a node per HMM-state if they do not exist yet
					NodeTempX *nodeAux = NULL;
					for(LArcTempX::iterator kt = nodePrev->lArcNext.begin() ; kt != nodePrev->lArcNext.end() ; ++kt) {
						if ((*kt)->nodeDest->state == state) {
							nodeAux = (*kt)->nodeDest;
							break;
						}
					}	
					
					// insertion penalty?
					char iIPIndex = -1;
					if (iState == 0) {
						assert(fPhoneIP[iPhone] != FLT_MAX);
						iIPIndex = iPhone;
					}
					
					if (nodeAux == NULL) {
						// create the node and the arc
						nodeAux = newNodeTempX(NODE_TYPE_HMM,m_iDepthFI+1+iState,state,iIPIndex);
						newArcTempX(nodePrev,nodeAux,ARC_TYPE_NULL,NULL,NULL);
						++iHMMNodes;
					}	
					if (iState == NUMBER_HMM_STATES-1) {
						nodeAux->bWordEnd = true;
					}
					// connect to the destination FI-node through a word-arc
					if (iState == NUMBER_HMM_STATES-1) {
					
						bool bFound = false;
						for(LArcTempX::iterator kt = nodeAux->lArcNext.begin() ; kt != nodeAux->lArcNext.end() ; ++kt) {
							if ((*kt)->nodeDest == nodeWord) {
								bFound = true;
								break;
							}
						}
						if (bFound == false) {
							newArcTempX(nodeAux,nodeWord,ARC_TYPE_NULL,NULL,NULL);
						}
					}
					nodePrev = nodeAux;
				}
			}
			delete [] iPhoneLeft;
			delete [] iPhoneRight;	
		}
		for(vector< pair<unsigned char*,NodeTempX*> >::iterator kt = vNodeFIDest.begin() ; 
				kt != vNodeFIDest.end() ; ++kt) {
			delete [] kt->first;
		}
		
		cout << "hmm-nodes: " << iHMMNodes << endl;
		
		cout << "elements: " << mFIWordNode.size();
		
		// backward-merge from word nodes
		mergeBackwardWordNodes(vNodeWord);
	}
	delete [] vLexUnitMonophone;
	delete [] iKeyFIDest;
	delete [] iKeyFI;
	
	double dTimeEndMonophones = TimeUtils::getTimeMilliseconds();
   
   printf("monophones: %f12.4\n",(dTimeEndMonophones-dTimeBeginMonophones)/1000.0);
   
   // checks
   int iIgnoredMI = 0;
   int iIgnoredFO = 0;

	for(MContextNode::iterator it = mContextEnd.begin() ; it != mContextEnd.end() ; ++it) {
		
		NodeTempX *nodeFO = it->second;
		if (nodeFO->lArcPrev.empty() && nodeFO->lArcNext.empty()) {
			++iIgnoredFO;
			// remobe the FO-node (it should have never been created!!)
			deleteNodeTempX(nodeFO);
		} else {
			assert(nodeFO->lArcPrev.empty() == false);
			assert(nodeFO->lArcNext.empty() == false);
		}
	}
	for(MContextNode::iterator it = mContextBeg.begin() ; it != mContextBeg.end() ; ++it) {
		
		NodeTempX *nodeMI = it->second;
		if (nodeMI->lArcPrev.empty() && nodeMI->lArcNext.empty()) {
			++iIgnoredMI;
			// remobe the MI-node (it should have never been created!!)
			deleteNodeTempX(nodeMI);
		} else {
			assert(nodeMI->lArcPrev.empty() == false);
			assert(nodeMI->lArcNext.empty() == false);
		}
	}
	printf("ignored MI: %d FO: %d\n",iIgnoredMI,iIgnoredFO);
   
   //exit(-1);
   
   //removeMINodes(mContextBeg);
   
   cout << "HMMs: " << m_hmmManager->getNumberHMMStatesPhysical() << " used: " << mHMMUsed.size() << endl;
	
 
   cout << "FI-nodes: " << mNodesFI.size() << endl;
   cout << "MI-nodes: " << mContextBeg.size() << endl;
   cout << "FO-nodes: " << mContextEnd.size() << endl;
   //exit(-1);
   
   // Forward-merge
   print();
   mergeForward();
   //print();
   mergeForward();
   //print();
	
	// Word-label pushing
   pushWordLabels();
   //print();
   
   // Backward-merge
   // note: two passes are typically necessary and enough, some FI-nodes merged in the first pass 
   //       enable merging of NULL-nodes in the second pass
   mergeBackward(); 
   //print();
   mergeBackward(); 
   //print();
   
   // FI removal (FI-nodes, which were used to build the network are no longer necessary)
  	removeFINodes();
   //print(); 
   mergeForward();
   mergeBackward();
   print();
      
   // build the language model look-ahead tree (needed to compute look-ahead scores)
   m_iLANodes = -1;
   m_iLATree = buildLMLookAheadTree(&m_iLANodes);
    
   // clean-up
   vector<unsigned char*> vTrash;
   for(MContextBool::iterator it = mContextLeft.begin() ; it != mContextLeft.end() ; ++it) {
   	vTrash.push_back(it->first);
   }
   for(MContextBool::iterator it = mContextRight.begin() ; it != mContextRight.end() ; ++it) {
   	vTrash.push_back(it->first);
   }
   for(MContextNode::iterator it = mContextBeg.begin() ; it != mContextBeg.end() ; ++it) {
   	vTrash.push_back(it->first);
   }
   for(MContextNode::iterator it = mContextEnd.begin() ; it != mContextEnd.end() ; ++it) {
   	vTrash.push_back(it->first);
   }
   for(MContextNode::iterator it = mNodesFI.begin() ; it != mNodesFI.end() ; ++it) {
   	vTrash.push_back(it->first);
   }
   for(MContextV::iterator it = mContextBegV.begin() ; it != mContextBegV.end() ; ++it) {
   	vTrash.push_back(it->first);
   	for(vector< pair<unsigned char*,NodeTempX*> >::iterator jt = it->second.begin() ; jt != it->second.end(); ++jt) {
   		vTrash.push_back(jt->first);
   	}
   }
   for(MContextV::iterator it = mContextEndV.begin() ; it != mContextEndV.end() ; ++it) {
   	vTrash.push_back(it->first);
   	for(vector< pair<unsigned char*,NodeTempX*> >::iterator jt = it->second.begin() ; jt != it->second.end(); ++jt) {
   		vTrash.push_back(jt->first);
   	}
   }
   for(vector<unsigned char*>::iterator it = vTrash.begin() ; it != vTrash.end() ; ++it) {
   	delete [] *it;
   }
	
	DynamicNetworkX *dynamicNetwork = compact();
   dynamicNetwork->setIP(fPhoneIP,iBasePhones);
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	printf("building time: %.2f seconds\n",dTimeSeconds);	

	return dynamicNetwork;
}

// push the word-labels to the beginning of the network
// only word-labels of polyphonic words are pushed, labels for monophonic words stay the same
// word labels of polyphonc wors originally go from an HMM-node to a FO-node
void NetworkBuilderX::pushWordLabels() {

	// (1) collect word-nodes
	MNodeTempX mNodeSeen;
	VNodeTempX vNodeWord;

	VNodeTempX vNodes;
	vNodes.push_back(m_nodeTempRoot);
	while(vNodes.empty() == false) {
		NodeTempX *nodeAux = vNodes.back();
		vNodes.pop_back();
		bool bWordArcs = false;
		bool bOtherArcs = false;		
		for(LArcTempX::iterator it = nodeAux->lArcNext.begin() ; it != nodeAux->lArcNext.end() ; ++it) {
			// add to the queue if not seen yet
			if (mNodeSeen.find((*it)->nodeDest) == mNodeSeen.end()) {
				vNodes.push_back((*it)->nodeDest);
				mNodeSeen.insert(MNodeTempX::value_type((*it)->nodeDest,true));
			}
			if ((*it)->iType == ARC_TYPE_WORD) {
				bWordArcs = true;
			} else {
				bOtherArcs = true;
			}
		}
		// we are not interested in word-nodes, since they are used for monophones and cannot be pushed
		// further to the head of the network
		if (nodeAux->iType == NODE_TYPE_WORD) {
			continue;
		}
		// only word arcs
		if ((bWordArcs == false) || (bOtherArcs)) {
			continue;
		}
		// all word arcs should point to the same destination node
		// if there are multiple word-arcs, shifting will only be possible if all of them have identical destination
		// node (homophonic words), otherwise shifting would introduce more arcs and nodes in the network	
		bool bShift = true;
		NodeTempX *nodeDestCommon = NULL;
		for(LArcTempX::iterator it = nodeAux->lArcNext.begin() ; it != nodeAux->lArcNext.end() ; ++it) {
			assert((*it)->iType == ARC_TYPE_WORD);
			if (nodeDestCommon == NULL) {
				nodeDestCommon	= (*it)->nodeDest;
			} else if ((*it)->nodeDest != nodeDestCommon) {
				bShift = false;
				break;
			}
		}
		if (bShift == false) {
			continue;
		}
		
		vNodeWord.push_back(nodeAux);	
	}
	
	cout << "-> word nodes: " << vNodeWord.size() << endl;
	
	// (2) shift each of the word arcs
	int iShiftsTotal = 0;
	for(VNodeTempX::iterator it = vNodeWord.begin() ; it != vNodeWord.end() ; ++it) {
	
		int iShifts = 0;
		if ((*it)->lArcPrev.size() > 1) {
			continue;
		}
		ArcTempX *arc = (*it)->lArcPrev.front();
		assert(arc->iType != ARC_TYPE_WORD);
		while((arc->nodeSource->lArcNext.size() == 1) && (arc->nodeSource->lArcPrev.size() == 1)) {
			if (arc->nodeSource->iDepth > m_iDepthMI) {
				arc = arc->nodeSource->lArcPrev.front();
				++iShifts;
			} else {
				break;
			}
		}
		// can't shift?
		if (arc->nodeDest == *it) {
			continue;
		}
		ArcTempX *arcReplace = arc;
		assert(arc->iType == ARC_TYPE_NULL);
		//assert(arcReplace->nodeSource->lArcNext.size() == 1);
		assert(arcReplace->nodeSource->iDepth >= m_iDepthMI);
		iShiftsTotal += iShifts*(*it)->lArcNext.size();
		
		assert((*it)->iType == NODE_TYPE_HMM);
		
		// a NULL node needs to be created so a word-arc does not go to a hmm-arc 
		// since we want to keep the hmm's in the arcs instead of in the nodes.
		// the reason is to prevent cache-misses during token propagation
		// a node is only checked if it actually gets activated
		NodeTempX *nodeNull = newNodeTempX(NODE_TYPE_NULL,arcReplace->nodeDest->iDepth);

		// (a) extract word-arcs from original place	
		assert((*it)->lArcPrev.size() == 1);
		VArcTempX vArcWords;
		NodeTempX *nodeDest = (*it)->lArcNext.front()->nodeDest;
		for(LArcTempX::iterator jt = (*it)->lArcNext.begin() ; jt != (*it)->lArcNext.end() ; ++jt) {
			assert((*jt)->nodeDest->iType == NODE_TYPE_FO);
			bool bFound = false;
			for(LArcTempX::iterator kt = (*jt)->nodeDest->lArcPrev.begin() ; kt != (*jt)->nodeDest->lArcPrev.end() ; ++kt) {
				if (*kt == *jt) {
					(*jt)->nodeDest->lArcPrev.erase(kt);
					bFound = true;
					break;
				}
			}	
			assert(bFound);
			// keep word arcs
			vArcWords.push_back(*jt);
			// change source and dest node
			(*jt)->nodeSource = arcReplace->nodeSource;
			(*jt)->nodeDest = nodeNull;
			nodeNull->lArcPrev.push_back(*jt);	
		}		
		(*it)->lArcNext.clear();
		// (b) extract null-arc from original place
		bool bFound = false;
		for(LArcTempX::iterator jt = arcReplace->nodeSource->lArcNext.begin() ; jt != arcReplace->nodeSource->lArcNext.end() ; ++jt) {
			if (*jt == arcReplace) {
				arcReplace->nodeSource->lArcNext.erase(jt);
				bFound = true;
				break;
			}
		}
		assert(bFound);
		for(VArcTempX::iterator jt = vArcWords.begin() ; jt != vArcWords.end() ; ++jt) {
			arcReplace->nodeSource->lArcNext.push_back(*jt);
		}
		assert(arcReplace->nodeDest->lArcPrev.size() == 1);	
		arcReplace->nodeDest->lArcPrev.clear();
		newArcTempX(nodeNull,arcReplace->nodeDest,ARC_TYPE_NULL,NULL,NULL);
		arcReplace->nodeSource = *it;
		arcReplace->nodeDest = nodeDest;
		(*it)->lArcNext.push_back(arcReplace);
		nodeDest->lArcPrev.push_back(arcReplace);	
	}
	
	printf("# shifts: %d\n",iShiftsTotal);
}

// compact the network
DynamicNetworkX *NetworkBuilderX::compact() {

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	// collect network nodes
	bool *bNodesSeen = new bool[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		bNodesSeen[i] = false;
	}
	int iArcs = 0;
	VNodeTempX vNodesAll;
	VNodeTempX vNodes;
	vNodesAll.push_back(m_nodeTempRoot);
	vNodes.push_back(m_nodeTempRoot);
	while(vNodes.empty() == false) {
		NodeTempX *nodeAux = vNodes.back();
		vNodes.pop_back();
		iArcs += nodeAux->lArcNext.size();	
		for(LArcTempX::iterator it = nodeAux->lArcNext.begin() ; it != nodeAux->lArcNext.end() ; ++it) {
			// add to the queue if not seen yet
			if (bNodesSeen[(*it)->nodeDest->iNode] == false) {	
				vNodesAll.push_back((*it)->nodeDest);
				vNodes.push_back((*it)->nodeDest);
				bNodesSeen[(*it)->nodeDest->iNode] = true;
			}
		}
	}
	delete [] bNodesSeen;

	// count the number of outgoing arcs in the network
	int iNodes = vNodesAll.size();
	int iArcsNull = 0;
	int iArcsHMM = 0;
	int iArcsWord = 0;
	
	// allocate memory
	DNode *nodes = new DNode[iNodes+1];		// extra node is needed to mark the ending arc
	DArc *arcs = new DArc[iArcs];
	
	// keep already created nodes (indexed by the corresponding temporal node id)
	int *iNodeCreated = new int[m_iNodesTempID];
	int *iFirstArc = new int[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		iNodeCreated[i] = -1;
		iFirstArc[i] = -1;
	}

	// fill the structures
	int iNode = 0;
	int iArc = 0;
	int iFOSeen = 0;
	for(VNodeTempX::iterator it = vNodesAll.begin() ; it != vNodesAll.end() ; ++it) {
		DNode *nodeAux;
		int iArcBase = -1;
		if (iNodeCreated[(*it)->iNode] == -1) {
			nodeAux = &nodes[iNode];
			iNodeCreated[(*it)->iNode] = iNode;
			iNode++;
			iArcBase = iArc;
			iArc += (*it)->lArcNext.size();
		} else {		
			nodeAux = &nodes[iNodeCreated[(*it)->iNode]]; 
			iArcBase = iFirstArc[iNodeCreated[(*it)->iNode]];
		}
		// map the node-type appropiately
		if ((*it)->iType == NODE_TYPE_FO){
			++iFOSeen;
		}
		switch((*it)->iType) {
			case NODE_TYPE_FI: 
			case NODE_TYPE_MI: 
			case NODE_TYPE_FO: 
			case NODE_TYPE_HMM:	
			{
				nodeAux->iType = NODE_TYPE_NULL;
				break;
			}
			default: {
				nodeAux->iType = (*it)->iType;	
			}
		}	
		
		nodeAux->iDepth = (*it)->iDepth;
		nodeAux->iIPIndex = (*it)->iIPIndex;
		nodeAux->bWordEnd = (*it)->bWordEnd;
		nodeAux->iArcNext = iArcBase;
		// store the outgoing arcs
		int i=0;
		for(LArcTempX::iterator jt = (*it)->lArcNext.begin() ; jt != (*it)->lArcNext.end() ; ++jt,++i) {
			if (iNodeCreated[(*jt)->nodeDest->iNode] == -1) {
				iNodeCreated[(*jt)->nodeDest->iNode] = iNode;
				iFirstArc[iNode] = iArc;
				iArc += (*jt)->nodeDest->lArcNext.size();
				iNode++;
			}
			// null-arc, if it goes to an HMM-node, then convert it to an HMM-arc and keep the hmm
			if ((*jt)->iType == ARC_TYPE_NULL) {
				if ((*jt)->nodeDest->iType == NODE_TYPE_HMM) {
					arcs[iArcBase+i].iType = ARC_TYPE_HMM;
					arcs[iArcBase+i].state = (*jt)->nodeDest->state;
					++iArcsHMM;
				} else {
					arcs[iArcBase+i].iType = ARC_TYPE_NULL;
					arcs[iArcBase+i].state = NULL;
					++iArcsNull;	
				}
			}
			// lex-unit
			else {
				assert((*jt)->iType == ARC_TYPE_WORD);
				arcs[iArcBase+i].iType = ARC_TYPE_WORD;
				arcs[iArcBase+i].lexUnit = (*jt)->lexUnit;
				++iArcsWord;
			} 
			arcs[iArcBase+i].iNodeDest = iNodeCreated[(*jt)->nodeDest->iNode];
			arcs[iArcBase+i].iLANode = (*jt)->iLANode;
		}
	}
	assert(iArc == iArcs);
	assert(iNode == iNodes);
	nodes[iNodes].iArcNext = iArcs;
	
	// clean-up
	for(VNodeTempX::iterator it = vNodesAll.begin() ; it != vNodesAll.end() ; ++it) {
		for(LArcTempX::iterator jt = (*it)->lArcNext.begin() ; jt != (*it)->lArcNext.end() ; ++jt) {
			deleteArcTempX(*jt);
		}
		deleteNodeTempX(*it);
	}	
	
	delete [] iNodeCreated;
	delete [] iFirstArc;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	cout << "- compacting summary ---------------------" << endl;
	cout << "# nodes:           " << setw(12) << iNodes << endl;
	cout << "# arcs:            " << setw(12) << iArcs << endl;
	cout << "  - null arcs:     " << setw(12) << iArcsNull << endl;
	cout << "  - hmm arcs:      " << setw(12) << iArcsHMM << endl;
	cout << "  - word arcs:     " << setw(12) << iArcsWord << endl;	
	cout << "network size:      " << setw(12) << (iNodes*sizeof(DNode)+iArcs*sizeof(DArc)) << " bytes" << endl;
	cout << "processing time:   " << setw(12) << dTimeSeconds << " seconds" << endl;
	cout << "------------------------------------------" << endl;

	return new DynamicNetworkX(nodes,iNodes,arcs,iArcs,m_iLATree,m_iLANodes);
}

// print network stats
void NetworkBuilderX::print() {

	cout << "-- temporal network stats ---------------" << endl;
	cout << "# nodes:       " << setw(6) << m_iNodesTemp << endl;
	cout << " -root-nodes   " << setw(6) << m_iNodesTempRoot << endl;
	cout << " -FI-nodes:    " << setw(6) << m_iNodesTempFI << endl;
	cout << " -MI-nodes:    " << setw(6) << m_iNodesTempMI << endl;
	cout << " -FO-nodes:    " << setw(6) << m_iNodesTempFO << endl;
	cout << " -hmm-nodes:   " << setw(6) << m_iNodesTempHMM << endl;
	cout << " -word-nodes:  " << setw(6) << m_iNodesTempWord << endl;
	cout << " -null-nodes:  " << setw(6) << m_iNodesTempNull << endl;
	cout << "-----------------------------------------" << endl;
	cout << "# arcs:        " << setw(6) << m_iArcsTemp << endl;
	cout << " -null-arcs:   " << setw(6) << m_iArcsNull << endl;
	cout << " -word-arcs:   " << setw(6) << m_iArcsWord << endl;
	cout << "-----------------------------------------" << endl;
}

// merge nodes forward
void NetworkBuilderX::mergeForward() {

	int iNodesProcessed = 0;
	int iNodesMerged = 0;
	int iArcsMerged = 0;
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	// allocate a structure to keep the node states	
	unsigned char *iNodeState = new unsigned char[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		iNodeState[i] = TEMPNODE_STATE_UNSEEN;
	}	
	
	// start from the FI-nodes
	LNodeTempX lNodes;
	for(LArcTempX::iterator it = m_nodeTempRoot->lArcNext.begin() ; it != m_nodeTempRoot->lArcNext.end() ; ++it) {
		lNodes.push_back((*it)->nodeDest);
		iNodeState[(*it)->nodeDest->iNode] = TEMPNODE_STATE_QUEUED;
	}	
	//cout << "total: " << m_nodeTempRoot->lArcNext.size() << endl;
	map<int,bool> mSeen;
	
	//lNodes.push_back(m_nodeTempRoot);
	//iNodeState[m_nodeTempRoot->iNode] = TEMPNODE_STATE_QUEUED;
	iNodeState[m_nodeTempRoot->iNode] = TEMPNODE_STATE_PROCESSED;
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		NodeTempX *node = lNodes.front();
		lNodes.pop_front();
		iNodeState[node->iNode] = TEMPNODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LArcTempX::iterator it = node->lArcNext.begin() ; it != node->lArcNext.end() ; ++it) {
			if ((*it)->nodeDest->iType == NODE_TYPE_FI) {
				mSeen[(*it)->nodeDest->iNode] = true;
			}	
			LArcTempX::iterator jt = it;
			jt++;
			for( ; jt != node->lArcNext.end() ; ) {
				if (equivalentNodes((*it)->nodeDest,(*jt)->nodeDest) && equivalentArcs(*it,*jt) 
					&& samePredecessors((*it)->nodeDest,(*jt)->nodeDest)) {
					// move successors from one arc to the other (if not already there)
					for(LArcTempX::iterator kt = (*jt)->nodeDest->lArcNext.begin() ; kt != (*jt)->nodeDest->lArcNext.end() ; ++kt) {
						bool bFound = false;
						for(LArcTempX::iterator lt = (*it)->nodeDest->lArcNext.begin() ; lt != (*it)->nodeDest->lArcNext.end() ; ++lt) {
							if (((*lt)->nodeDest == (*kt)->nodeDest) && (equivalentArcs(*lt,*kt))) {
								bFound = true;
								break;
							}
						}	
						// move the arc to the other that will stay
						if (bFound == false) {
							(*kt)->nodeSource = (*it)->nodeDest;
							(*it)->nodeDest->lArcNext.push_back(*kt);
						} else {
							for(LArcTempX::iterator ht = (*kt)->nodeDest->lArcPrev.begin() ; ht != (*kt)->nodeDest->lArcPrev.end() ; ++ht) {
								if (*ht == *kt) {
									(*kt)->nodeDest->lArcPrev.erase(ht);
									deleteArcTempX(*kt);
									++iArcsMerged;
									break;
								}
							}	
						}
					}
					// remove links to predecessors from one node to the other
					for(LArcTempX::iterator kt = (*jt)->nodeDest->lArcPrev.begin() ; kt != (*jt)->nodeDest->lArcPrev.end() ; ++kt) {
						if (*kt == *jt) {
							continue;
						}
						for(LArcTempX::iterator lt = (*kt)->nodeSource->lArcNext.begin() ; lt != (*kt)->nodeSource->lArcNext.end() ; ++lt) {
							if (*lt == *kt) {
								ArcTempX *arcDelete = *lt;
								(*kt)->nodeSource->lArcNext.erase(lt);
								deleteArcTempX(arcDelete);
								++iArcsMerged;
								break;
							}
						}	
					}
					// remove the redundant arc and node
					if (iNodeState[(*it)->nodeDest->iNode] == TEMPNODE_STATE_PROCESSED) {
						iNodeState[(*it)->nodeDest->iNode] = TEMPNODE_STATE_UNSEEN;
					}	
					//assert((*jt)->nodeDest->bWordEnd == (*it)->nodeDest->bWordEnd);
					if ((*jt)->nodeDest->bWordEnd) {
						(*it)->nodeDest->bWordEnd = true;
					}
					deleteNodeTempX((*jt)->nodeDest);
					deleteArcTempX(*jt);
					jt = node->lArcNext.erase(jt);
					// stats
					++iNodesMerged;
					++iArcsMerged;
				} else {
					++jt;
				}
			}
		}
		
		// put in the queue those successor nodes that have all their predecessors processed
		for(LArcTempX::iterator it = node->lArcNext.begin() ; it != node->lArcNext.end() ; ++it) {
			if (iNodeState[(*it)->nodeDest->iNode] != TEMPNODE_STATE_UNSEEN) {
				continue;
			}
			assert((*it)->nodeDest->iNode < m_iNodesTempID);
			assert((*it)->nodeDest->iNode >= 0);	
			assert(iNodeState[(*it)->nodeDest->iNode] == TEMPNODE_STATE_UNSEEN);
			bool bReady = true;
			for(LArcTempX::iterator jt = (*it)->nodeDest->lArcPrev.begin() ; jt != (*it)->nodeDest->lArcPrev.end() ; ++jt) {	
				if (iNodeState[(*jt)->nodeSource->iNode] != TEMPNODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}
			if (bReady /*|| ((*it)->nodeDest->iType == NODE_TYPE_FI)*/) {
				lNodes.push_back((*it)->nodeDest);
				iNodeState[(*it)->nodeDest->iNode] = TEMPNODE_STATE_QUEUED;
			}
		}
	}	
	
	delete [] iNodeState;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	cout << "- forward merging summary ----------------" << endl;
	cout << "# nodes processed: " << setw(12) << iNodesProcessed << endl;
	cout << "# nodes merged:    " << setw(12) << iNodesMerged << endl;
	cout << "# arcs merged:     " << setw(12) << iArcsMerged << endl;
	cout << "processing time:   " << setw(12) << dTimeSeconds << " seconds" << endl;	
	//cout << "fi seen:           " << setw(12) << mSeen.size() << endl;
}


// remove MI nodes used to build the network
void NetworkBuilderX::removeMINodes(MContextNode &mContextBeg) {

	for(MContextNode::iterator it = mContextBeg.begin() ; it != mContextBeg.end() ; ++it) {
		
		NodeTempX *nodeMI = it->second;
		
		// remove links to the MI-node
		assert(nodeMI->lArcNext.empty());
		
		for(LArcTempX::iterator jt = nodeMI->lArcPrev.begin() ; jt != nodeMI->lArcPrev.end() ; ++jt) {
			bool bFound = false;
			for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcNext.begin() ; kt != (*jt)->nodeSource->lArcNext.end() ; ) {
				if (((*kt)->nodeDest == nodeMI) && ((*kt)->lexUnit == (*jt)->lexUnit)) {
					deleteArcTempX(*kt);
					kt = (*jt)->nodeSource->lArcNext.erase(kt);
					assert(bFound == false);
					bFound = true;
				} else {
					++kt;
				}	
			}
			assert(bFound);	
		}
				
		deleteNodeTempX(nodeMI);
	}
}

// remove FI nodes used to build the network
void NetworkBuilderX::removeFINodes() {

	int iNodesProcessed = 0;
	int iNodesMerged = 0;
	int iArcsMerged = 0;
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	// get the FI-nodes
	VNodeTempX vNodeFI;
	for(LArcTempX::iterator it = m_nodeTempRoot->lArcNext.begin() ; it != m_nodeTempRoot->lArcNext.end() ; ++it) {
		vNodeFI.push_back((*it)->nodeDest);
	}
	
	// (1) remove FI-nodes
	for(VNodeTempX::iterator it = vNodeFI.begin() ; it != vNodeFI.end() ; ++it) {
		
		NodeTempX *nodeFI = *it;
		
		// remove links coming from the FI-node
		VNodeTempX vNodeDest;
		for(LArcTempX::iterator jt = nodeFI->lArcNext.begin() ; jt != nodeFI->lArcNext.end() ; ) {
			assert((*jt)->iType == ARC_TYPE_NULL);
			bool bFound = false;
			for(LArcTempX::iterator kt = (*jt)->nodeDest->lArcPrev.begin() ; kt != (*jt)->nodeDest->lArcPrev.end() ; ) {
				if ((*kt)->nodeSource == nodeFI) {
					kt = (*jt)->nodeDest->lArcPrev.erase(kt);
					assert(bFound == false);
					bFound = true;
				} else {
					++kt;
				}
			}
			assert(bFound);	
			vNodeDest.push_back((*jt)->nodeDest);
			ArcTempX *arcToDelete = (*jt);
			jt = nodeFI->lArcNext.erase(jt);
			deleteArcTempX(arcToDelete);
		}
		// remove links going to the FI-node
		VNodeTempX vNodeSource;
		for(LArcTempX::iterator jt = nodeFI->lArcPrev.begin() ; jt != nodeFI->lArcPrev.end() ; ) {
			assert((*jt)->iType == ARC_TYPE_NULL);
			bool bFound = false;
			for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcNext.begin() ; kt != (*jt)->nodeSource->lArcNext.end() ; ) {
				if (((*kt)->nodeDest == nodeFI) && ((*kt)->lexUnit == (*jt)->lexUnit)) {
					kt = (*jt)->nodeSource->lArcNext.erase(kt);
					assert(bFound == false);
					bFound = true;
				} else {
					++kt;
				}	
			}
			assert(bFound);	
			vNodeSource.push_back((*jt)->nodeSource);
			ArcTempX *arcToDelete = (*jt);
			jt = nodeFI->lArcPrev.erase(jt);
			deleteArcTempX(arcToDelete);
		}
			
		// make direct connections
		for(VNodeTempX::iterator it = vNodeSource.begin() ; it != vNodeSource.end() ; ++it) {
			for(VNodeTempX::iterator jt = vNodeDest.begin() ; jt != vNodeDest.end() ; ++jt) {
				if (areConnected(*it,*jt) == false) {
					newArcTempX(*it,*jt,ARC_TYPE_NULL,NULL,NULL);
				}
			}
		}
		
		deleteNodeTempX(nodeFI);
	}
		
	// get the auxiliar nodes
	
	// allocate a structure to keep the node states	
	map<int,unsigned char*> mNodeState;
	
	// start from the root-node
	/*LNodeTempX lNodes;
	lNodes.push_back(m_nodeTempRoot);
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		NodeTempX *node = lNodes.front();
		lNodes.pop_front();
		mNodeState[node->iNode] = TEMPNODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LArcTempX::iterator it = node->lArcNext.begin() ; it != node->lArcNext.end() ; ++it) {
			
			
		
		
		}
		
		// put in the queue those successor nodes that have all their predecessors processed
		for(LArcTempX::iterator it = node->lArcNext.begin() ; it != node->lArcNext.end() ; ++it) {
			if ((iNodeState[(*it)->nodeDest->iNode] == TEMPNODE_STATE_QUEUED) || 
				(iNodeState[(*it)->nodeDest->iNode] == TEMPNODE_STATE_PROCESSED)) {
				continue;
			}
			assert((*it)->nodeDest->iNode < m_iNodesTempID);
			assert((*it)->nodeDest->iNode >= 0);	
			assert(iNodeState[(*it)->nodeDest->iNode] == TEMPNODE_STATE_UNSEEN);
			bool bReady = true;
			for(LArcTempX::iterator jt = (*it)->nodeDest->lArcPrev.begin() ; jt != (*it)->nodeDest->lArcPrev.end() ; ++jt) {	
				if (iNodeState[(*jt)->nodeSource->iNode] != TEMPNODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}
			if ((bReady == true) || ((*it)->nodeDest->iType == NODE_TYPE_FI)) {
				lNodes.push_back((*it)->nodeDest);
				iNodeState[(*it)->nodeDest->iNode] = TEMPNODE_STATE_QUEUED;
			}
		}

	}	*/
	
	//delete [] iNodeState;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	printf("- removal summary ----------------\n");
	printf("# nodes processed: %12d\n",iNodesProcessed);
	printf("# nodes merged:    %12d\n",iNodesMerged);
	printf("# arcs merged:     %12d\n",iArcsMerged);
	printf("processing time:   %12.2f seconds\n",dTimeSeconds);
}


// merge nodes backward
void NetworkBuilderX::mergeBackward() {

	int iNodesProcessed = 0;
	int iNodesMerged = 0;
	int iArcsMerged = 0;
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// allocate a structure to keep the node states	
	unsigned char *iNodeState = new unsigned char[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		iNodeState[i] = TEMPNODE_STATE_UNSEEN;
	}	
	
	// start from the FI-nodes
	LNodeTempX lNodes;
	for(LArcTempX::iterator it = m_nodeTempRoot->lArcNext.begin() ; it != m_nodeTempRoot->lArcNext.end() ; ++it) {
		//assert((*it)->nodeDest->iType == NODE_TYPE_FI);
		lNodes.push_back((*it)->nodeDest);
		iNodeState[(*it)->nodeDest->iNode] = TEMPNODE_STATE_QUEUED;	
	}
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		NodeTempX *node = lNodes.front();
		lNodes.pop_front();
		iNodeState[node->iNode] = TEMPNODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			LArcTempX::iterator jt = it;
			jt++;
			for( ; jt != node->lArcPrev.end() ; ) {
				if (equivalentNodes((*it)->nodeSource,(*jt)->nodeSource) && equivalentArcs(*it,*jt) && 
					sameSuccessors((*it)->nodeSource,(*jt)->nodeSource)) {
					// move predecessors from one arc to the other (if not already there)
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcPrev.begin() ; kt != (*jt)->nodeSource->lArcPrev.end() ; ++kt) {
						bool bFound = false;
						for(LArcTempX::iterator lt = (*it)->nodeSource->lArcPrev.begin() ; lt != (*it)->nodeSource->lArcPrev.end() ; ++lt) {
							if (((*lt)->nodeSource == (*kt)->nodeSource) && (equivalentArcs(*lt,*kt))) {
								bFound = true;
								break;
							}
						}	
						if (bFound == false) {
							(*kt)->nodeDest = (*it)->nodeSource;
							(*it)->nodeSource->lArcPrev.push_back(*kt);
						} else {
							for(LArcTempX::iterator ht = (*kt)->nodeSource->lArcNext.begin() ; ht != (*kt)->nodeSource->lArcNext.end() ; ++ht) {
								if (*ht == *kt) {
									(*kt)->nodeSource->lArcNext.erase(ht);
									deleteArcTempX(*kt);
									++iArcsMerged;
									break;
								}
							}
						}
					}
					// remove links
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcNext.begin() ; kt != (*jt)->nodeSource->lArcNext.end() ; ++kt) {
						if (*kt == *jt) {
							continue;
						}
						for(LArcTempX::iterator lt = (*kt)->nodeDest->lArcPrev.begin() ; lt != (*kt)->nodeDest->lArcPrev.end() ; ++lt) {
							if (*lt == *kt) {
								ArcTempX *arcDelete = *lt;
								(*kt)->nodeDest->lArcPrev.erase(lt);
								deleteArcTempX(arcDelete);
								++iArcsMerged;
								break;
							}
						}	
					}
					// remove the redundant arc and node
					assert((*jt)->nodeSource->bWordEnd == (*it)->nodeSource->bWordEnd);
					(*it)->nodeSource->iDepth = max((*it)->nodeSource->iDepth,(*jt)->nodeSource->iDepth);		// keep highest depth
					deleteNodeTempX((*jt)->nodeSource);
					deleteArcTempX(*jt);
					jt = node->lArcPrev.erase(jt);
					++iNodesMerged;
					++iArcsMerged;
				} else {
					++jt;
				}
			}
		}
		
		// put in the queue those predecessor nodes that have all their successors processed
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			if (iNodeState[(*it)->nodeSource->iNode] != TEMPNODE_STATE_UNSEEN) {
				continue;
			}
			bool bReady = true;
			for(LArcTempX::iterator jt = (*it)->nodeSource->lArcNext.begin() ; jt != (*it)->nodeSource->lArcNext.end() ; ++jt) {	
				if (iNodeState[(*jt)->nodeDest->iNode] != TEMPNODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}
			if (bReady) {
				lNodes.push_back((*it)->nodeSource);
				iNodeState[(*it)->nodeSource->iNode] = TEMPNODE_STATE_QUEUED;
			}
		}
			
		//printf("queue: %12d nodes processed: %12d merged: %12d\n",lNodes.size(),iNodesProcessed,iNodesMerged);
	}	

	
	delete [] iNodeState;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	printf("- backward merging summary ---------------\n");
	printf("# nodes processed: %12d\n",iNodesProcessed);
	printf("# nodes merged:    %12d\n",iNodesMerged);
	printf("# arcs merged:     %12d\n",iArcsMerged);
	printf("processing time:   %12.2f seconds\n",dTimeSeconds);
}

// merge nodes backward starting from the given group of MI-nodes
void NetworkBuilderX::mergeBackwardMI(VNodeTempX &vNodeTempMI, NodeMergeInfo *fiNodeInfo) {

	int iNodesProcessed = 0;
	int iNodesMerged = 0;
	int iArcsMerged = 0;
	
	//double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// allocate a structure to keep the node states	
	//unsigned char *iNodeState = new unsigned char[m_iNodesTemp];
	//for(int i=0 ; i<m_iNodesTemp ; ++i) {
	unsigned char *iNodeState = new unsigned char[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		iNodeState[i] = TEMPNODE_STATE_UNSEEN;
	}	
	
	// start from the MI-nodes
	LNodeTempX lNodes;
	for(VNodeTempX::iterator it = vNodeTempMI.begin() ; it != vNodeTempMI.end() ; ++it) {
		assert((*it)->iType == NODE_TYPE_MI);
		lNodes.push_back(*it);
		iNodeState[(*it)->iNode] = TEMPNODE_STATE_QUEUED;	
	}
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		NodeTempX *node = lNodes.front();
		lNodes.pop_front();
		iNodeState[node->iNode] = TEMPNODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			LArcTempX::iterator jt = it;
			jt++;
			for( ; jt != node->lArcPrev.end() ; ) {
				if (equivalentNodes((*it)->nodeSource,(*jt)->nodeSource) && sameSuccessors((*it)->nodeSource,(*jt)->nodeSource)) {
					// move predecessors from one arc to the other (if not already there)
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcPrev.begin() ; kt != (*jt)->nodeSource->lArcPrev.end() ; ++kt) {
						bool bFound = false;
						for(LArcTempX::iterator lt = (*it)->nodeSource->lArcPrev.begin() ; lt != (*it)->nodeSource->lArcPrev.end() ; ++lt) {
							if (((*lt)->nodeSource == (*kt)->nodeSource) && equivalentArcs(*lt,*kt)){
								bFound = true;
								break;
							}
						}	
						if (bFound == false) {
							(*kt)->nodeDest = (*it)->nodeSource;
							(*it)->nodeSource->lArcPrev.push_back(*kt);
						} else {
							for(LArcTempX::iterator ht = (*kt)->nodeSource->lArcNext.begin() ; ht != (*kt)->nodeSource->lArcNext.end() ; ++ht) {
								if (*ht == *kt) {
									(*kt)->nodeSource->lArcNext.erase(ht);
									deleteArcTempX(*kt);
									++iArcsMerged;
									break;
								}
							}
						}
					}
					// remove links
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcNext.begin() ; kt != (*jt)->nodeSource->lArcNext.end() ; ++kt) {
						if (*kt == *jt) {
							continue;
						}
						for(LArcTempX::iterator lt = (*kt)->nodeDest->lArcPrev.begin() ; lt != (*kt)->nodeDest->lArcPrev.end() ; ++lt) {
							if (*lt == *kt) {
								ArcTempX *arcDelete = *lt;
								(*kt)->nodeDest->lArcPrev.erase(lt);
								deleteArcTempX(arcDelete);
								++iArcsMerged;
								break;
							}
						}	
					}
					// remove the redundant arc and node
					if ((*jt)->nodeSource->iType == NODE_TYPE_FI) {
						fiNodeInfo[(*jt)->nodeSource->iFI].iMerged = (*it)->nodeSource->iFI;
					}
					deleteNodeTempX((*jt)->nodeSource);
					deleteArcTempX(*jt);
					jt = node->lArcPrev.erase(jt);
					++iNodesMerged;
					++iArcsMerged;
				} else {
					++jt;
				}
			}
		}
		
		// do not add FI-nodes to the queue
		if (node->iDepth == 1) {
			continue;
		}
		
		// put in the queue those predecessor nodes that have all their successors processed
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			if (iNodeState[(*it)->nodeSource->iNode] != TEMPNODE_STATE_UNSEEN) {
				continue;
			}
			bool bReady = true;
			for(LArcTempX::iterator jt = (*it)->nodeSource->lArcNext.begin() ; jt != (*it)->nodeSource->lArcNext.end() ; ++jt) {	
				if (iNodeState[(*jt)->nodeDest->iNode] != TEMPNODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}
			if (bReady) {
				lNodes.push_back((*it)->nodeSource);
				iNodeState[(*it)->nodeSource->iNode] = TEMPNODE_STATE_QUEUED;
			}
		}
			
		//printf("queue: %12d nodes processed: %12d merged: %12d\n",lNodes.size(),iNodesProcessed,iNodesMerged);
	}	
	
	delete [] iNodeState;
	
	//double dTimeEnd = TimeUtils::getTimeMilliseconds();
	//double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	/*printf("- backward merging summary ---------------\n");
	printf("# nodes processed: %12d\n",iNodesProcessed);
	printf("# nodes merged:    %12d\n",iNodesMerged);
	printf("# arcs merged:     %12d\n",iArcsMerged);
	printf("processing time:   %12.2f seconds\n",dTimeSeconds);*/
}

// merge nodes backward starting from the given group of FI-nodes
void NetworkBuilderX::mergeBackwardFI(VNodeTempX &vNodeTempFI, NodeMergeInfo *foNodeMergeInfo) {

	int iNodesProcessed = 0;
	int iNodesMerged = 0;
	int iArcsMerged = 0;
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// allocate a structure to keep the node states	
	//unsigned char *iNodeState = new unsigned char[10*m_iNodesTemp];
	//for(int i=0 ; i<10*m_iNodesTemp ; ++i) {
	unsigned char *iNodeState = new unsigned char[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		iNodeState[i] = TEMPNODE_STATE_UNSEEN;
	}	
	
	// start from the MI-nodes
	LNodeTempX lNodes;
	for(VNodeTempX::iterator it = vNodeTempFI.begin() ; it != vNodeTempFI.end() ; ++it) {
		assert((*it)->iType == NODE_TYPE_FI);
		lNodes.push_back(*it);
		iNodeState[(*it)->iNode] = TEMPNODE_STATE_QUEUED;	
	}
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		NodeTempX *node = lNodes.front();
		lNodes.pop_front();
		iNodeState[node->iNode] = TEMPNODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			LArcTempX::iterator jt = it;
			jt++;
			for( ; jt != node->lArcPrev.end() ; ) {
				if (equivalentNodes((*it)->nodeSource,(*jt)->nodeSource) && sameSuccessors((*it)->nodeSource,(*jt)->nodeSource)) {
					// move predecessors from one arc to the other (if not already there)
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcPrev.begin() ; kt != (*jt)->nodeSource->lArcPrev.end() ; ++kt) {
						bool bFound = false;
						for(LArcTempX::iterator lt = (*it)->nodeSource->lArcPrev.begin() ; lt != (*it)->nodeSource->lArcPrev.end() ; ++lt) {
							if (((*lt)->nodeSource == (*kt)->nodeSource) && equivalentArcs(*lt,*kt)) {
								bFound = true;
								break;
							}
						}	
						if (bFound == false) {
							(*kt)->nodeDest = (*it)->nodeSource;
							(*it)->nodeSource->lArcPrev.push_back(*kt);
						} else {
							for(LArcTempX::iterator ht = (*kt)->nodeSource->lArcNext.begin() ; ht != (*kt)->nodeSource->lArcNext.end() ; ++ht) {
								if (*ht == *kt) {
									(*kt)->nodeSource->lArcNext.erase(ht);
									deleteArcTempX(*kt);
									++iArcsMerged;
									break;
								}
							}
						}
					}
					// remove links
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcNext.begin() ; kt != (*jt)->nodeSource->lArcNext.end() ; ++kt) {
						if (*kt == *jt) {
							continue;
						}
						for(LArcTempX::iterator lt = (*kt)->nodeDest->lArcPrev.begin() ; lt != (*kt)->nodeDest->lArcPrev.end() ; ++lt) {
							if (*lt == *kt) {
								ArcTempX *arcDelete = *lt;
								(*kt)->nodeDest->lArcPrev.erase(lt);
								deleteArcTempX(arcDelete);
								++iArcsMerged;
								break;
							}
						}	
					}
					assert((*jt)->nodeSource->bWordEnd == (*it)->nodeSource->bWordEnd);
					// remove the redundant arc and node
					if ((*jt)->nodeSource->iType == NODE_TYPE_FO) {
						foNodeMergeInfo[(*jt)->nodeSource->iFI].iMerged = (*it)->nodeSource->iFI;
					}
					deleteNodeTempX((*jt)->nodeSource);
					deleteArcTempX(*jt);
					jt = node->lArcPrev.erase(jt);
					++iNodesMerged;
					++iArcsMerged;
				} else {
					++jt;
				}
			}
		}
		
		// do not add FI-nodes to the queue
		if (node->iType == NODE_TYPE_FO) {
			continue;
		}
		
		// put in the queue those predecessor nodes that have all their successors processed
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			if (iNodeState[(*it)->nodeSource->iNode] != TEMPNODE_STATE_UNSEEN) {
				continue;
			}
			bool bReady = true;
			for(LArcTempX::iterator jt = (*it)->nodeSource->lArcNext.begin() ; jt != (*it)->nodeSource->lArcNext.end() ; ++jt) {	
				if (iNodeState[(*jt)->nodeDest->iNode] != TEMPNODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}
			if (bReady) {
				lNodes.push_back((*it)->nodeSource);
				iNodeState[(*it)->nodeSource->iNode] = TEMPNODE_STATE_QUEUED;
			}
		}
			
		//printf("queue: %12d nodes processed: %12d merged: %12d\n",lNodes.size(),iNodesProcessed,iNodesMerged);
	}	
	
	delete [] iNodeState;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	printf("- backward merging summary ---------------\n");
	printf("# nodes processed: %12d\n",iNodesProcessed);
	printf("# nodes merged:    %12d\n",iNodesMerged);
	printf("# arcs merged:     %12d\n",iArcsMerged);
	printf("processing time:   %12.2f seconds\n",dTimeSeconds);
}

// merge nodes backward starting from the given group of MI-nodes
void NetworkBuilderX::mergeBackwardWordNodes(VNodeTempX &vNodeWord) {

	int iNodesProcessed = 0;
	int iNodesMerged = 0;
	int iArcsMerged = 0;
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
	// allocate a structure to keep the node states	
	//unsigned char *iNodeState = new unsigned char[m_iNodesTemp];
	//for(int i=0 ; i<m_iNodesTemp ; ++i) {
	unsigned char *iNodeState = new unsigned char[m_iNodesTempID];
	for(int i=0 ; i<m_iNodesTempID ; ++i) {
		iNodeState[i] = TEMPNODE_STATE_UNSEEN;
	}	
	
	// start from the word-nodes
	LNodeTempX lNodes;
	for(VNodeTempX::iterator it = vNodeWord.begin() ; it != vNodeWord.end() ; ++it) {
		assert((*it)->iType == NODE_TYPE_WORD);
		lNodes.push_back(*it);
		iNodeState[(*it)->iNode] = TEMPNODE_STATE_QUEUED;	
	}
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		NodeTempX *node = lNodes.front();
		lNodes.pop_front();
		iNodeState[node->iNode] = TEMPNODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			LArcTempX::iterator jt = it;
			jt++;
			for( ; jt != node->lArcPrev.end() ; ) {
				if (equivalentNodes((*it)->nodeSource,(*jt)->nodeSource) && sameSuccessors((*it)->nodeSource,(*jt)->nodeSource)) {
					// move predecessors from one arc to the other (if not already there)
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcPrev.begin() ; kt != (*jt)->nodeSource->lArcPrev.end() ; ++kt) {
						bool bFound = false;
						for(LArcTempX::iterator lt = (*it)->nodeSource->lArcPrev.begin() ; lt != (*it)->nodeSource->lArcPrev.end() ; ++lt) {
							if (((*lt)->nodeSource == (*kt)->nodeSource) && equivalentArcs(*lt,*kt)){
								bFound = true;
								break;
							}
						}	
						if (bFound == false) {
							(*kt)->nodeDest = (*it)->nodeSource;
							(*it)->nodeSource->lArcPrev.push_back(*kt);
						} else {
							for(LArcTempX::iterator ht = (*kt)->nodeSource->lArcNext.begin() ; ht != (*kt)->nodeSource->lArcNext.end() ; ++ht) {
								if (*ht == *kt) {
									(*kt)->nodeSource->lArcNext.erase(ht);
									deleteArcTempX(*kt);
									++iArcsMerged;
									break;
								}
							}
						}
					}
					// remove links
					for(LArcTempX::iterator kt = (*jt)->nodeSource->lArcNext.begin() ; kt != (*jt)->nodeSource->lArcNext.end() ; ++kt) {
						if (*kt == *jt) {
							continue;
						}
						for(LArcTempX::iterator lt = (*kt)->nodeDest->lArcPrev.begin() ; lt != (*kt)->nodeDest->lArcPrev.end() ; ++lt) {
							if (*lt == *kt) {
								ArcTempX *arcDelete = *lt;
								(*kt)->nodeDest->lArcPrev.erase(lt);
								deleteArcTempX(arcDelete);
								++iArcsMerged;
								break;
							}
						}	
					}
					assert((*jt)->nodeSource->bWordEnd == (*it)->nodeSource->bWordEnd);
					// remove the redundant arc and node
					assert((*jt)->nodeSource->iType != NODE_TYPE_FI);
					deleteNodeTempX((*jt)->nodeSource);
					deleteArcTempX(*jt);
					jt = node->lArcPrev.erase(jt);
					++iNodesMerged;
					++iArcsMerged;
				} else {
					++jt;
				}
			}
		}
		
		// do not add FI-nodes to the queue
		if (node->iDepth == 3) {
			continue;
		}
		
		// put in the queue those predecessor nodes that have all their successors processed
		for(LArcTempX::iterator it = node->lArcPrev.begin() ; it != node->lArcPrev.end() ; ++it) {
			if (iNodeState[(*it)->nodeSource->iNode] != TEMPNODE_STATE_UNSEEN) {
				continue;
			}
			bool bReady = true;
			for(LArcTempX::iterator jt = (*it)->nodeSource->lArcNext.begin() ; jt != (*it)->nodeSource->lArcNext.end() ; ++jt) {	
				if (iNodeState[(*jt)->nodeDest->iNode] != TEMPNODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}
			if (bReady) {
				lNodes.push_back((*it)->nodeSource);
				iNodeState[(*it)->nodeSource->iNode] = TEMPNODE_STATE_QUEUED;
			}
		}
			
		//printf("queue: %12d nodes processed: %12d merged: %12d\n",lNodes.size(),iNodesProcessed,iNodesMerged);
	}	
	
	delete [] iNodeState;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	printf("- backward merging summary ---------------\n");
	printf("# nodes processed: %12d\n",iNodesProcessed);
	printf("# nodes merged:    %12d\n",iNodesMerged);
	printf("# arcs merged:     %12d\n",iArcsMerged);
	printf("processing time:   %12.2f seconds\n",dTimeSeconds);
}


// build the look-ahead tree from the network and transform it to a topologically sorted array
// -leaves correspond to words in the vocabulary
// -positions in the array keep the index of the parent node
int *NetworkBuilderX::buildLMLookAheadTree(int *iNodes) {

	/*word-arcs have been moved to the head of the network during puhsWordLabels, however they are always placed after
	 the first HMM-model (first three HMM-states). Collect all the word arcs and make sure there
	is at least one arc per word in the vocabulary. If the vocabulary size is V, then there must 
	
	- each node within the network should have an integer with the LMLA element in the LMLA-tree (which in
	 reality is a topologically sorted array). When a word-arc is accessed then the index within the LMLA-tree can 
	 be the lexical-unit-id, sometimes 
	
	two goals:
	
	a) efficient copy of lm-scores from n-gram tables to lm-context-dependent-lmla-tree, the arcs in the LM-state
	are sorted by lexUnitId (check this), just use the lex-unit-id to copy the score scores[lexUnitId] = score, first for the deepest backoff, last is the higher ngrams.*/

	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	// (1) collect word-arcs and initialize indices of look-ahead nodes
	MNodeTempX mNodeSeen;
	VArcTempX vArcWord;

	VNodeTempX vNodes;
	vNodes.push_back(m_nodeTempRoot);
	int iArcsInitialized = 0;
	while(vNodes.empty() == false) {
		NodeTempX *nodeAux = vNodes.back();
		vNodes.pop_back();
		for(LArcTempX::iterator it = nodeAux->lArcNext.begin() ; it != nodeAux->lArcNext.end() ; ++it) {
			// add to the queue if not seen yet
			if (mNodeSeen.find((*it)->nodeDest) == mNodeSeen.end()) {
				vNodes.push_back((*it)->nodeDest);
				mNodeSeen.insert(MNodeTempX::value_type((*it)->nodeDest,true));
			}
			if ((*it)->iType == ARC_TYPE_WORD) {
				vArcWord.push_back(*it);	
			}
			// initialize the index to the look-ahead node
			(*it)->iLANode = -1;
			++iArcsInitialized;
		}
	}
	
	cout << "arcs initialized: " << iArcsInitialized << endl;
	
	// (2) traverse the graph backwards starting at the word-arcs	
	// tree traversal is breadth-first so look-ahead indices are topologically sorted
	int iLANode = 0;
	map<int,int> mPredecessor;
	DArcTempX dArcs;
	for(VArcTempX::iterator it = vArcWord.begin() ; it != vArcWord.end() ; ++it) {
		dArcs.push_back(*it);
		(*it)->iLANode = iLANode++;
	}
	int iArcsProcessed = 0;
	while(dArcs.empty() == false) {
		ArcTempX *arcAux = dArcs.front();
		dArcs.pop_front();
		++iArcsProcessed;
		// a) one successor: same id
		if (arcAux->nodeSource->lArcNext.size() == 1) {
			for(LArcTempX::iterator it = arcAux->nodeSource->lArcPrev.begin() ; it != arcAux->nodeSource->lArcPrev.end() ; ++it) {
				assert((*it)->iLANode == -1);
				(*it)->iLANode = arcAux->iLANode;
				assert((*it)->iLANode != -1);
				if ((*it)->nodeSource->iDepth > m_iDepthMI) {
					dArcs.push_front(*it);		// deep first (successors are already processed)
				}
			}
		} 
		// b) multiple successors: different id (if not already assigned)
		else {
			if (arcAux->nodeSource->lArcPrev.front()->iLANode == -1) {
				// make sure all successors have already a look-ahead id, otherwise put back in the list
				bool bReady = true;
				for(LArcTempX::iterator it = arcAux->nodeSource->lArcNext.begin() ; it != arcAux->nodeSource->lArcNext.end() ; ++it) {
					if ((*it)->iLANode == -1) {	
						bReady = false;
						break;
					}
				}
				if (bReady == false) {
					dArcs.push_back(arcAux);		// to the back of the queue, successors need to be processed first
					continue;
				}	
				assert(arcAux->nodeSource->lArcPrev.empty() == false);
				for(LArcTempX::iterator it = arcAux->nodeSource->lArcPrev.begin() ; it != arcAux->nodeSource->lArcPrev.end() ; ++it) {
					assert((*it)->iLANode == -1);
					(*it)->iLANode = iLANode;
					if ((*it)->nodeSource->iDepth > m_iDepthMI) {
						dArcs.push_front(*it);		// deep first (successors are already processed)
					}
				}
				for(LArcTempX::iterator it = arcAux->nodeSource->lArcNext.begin() ; it != arcAux->nodeSource->lArcNext.end() ; ++it) {
					assert(mPredecessor.find((*it)->iLANode) == mPredecessor.end());
					mPredecessor[(*it)->iLANode] = iLANode;
				}
				++iLANode;
			}	
		}
	}	

	// build the array keeping a topological order of the tree
	int *iLATree = new int[iLANode];
	for(int i=0 ; i < iLANode ; ++i) {
		iLATree[i] = -1;
	}
	for(map<int,int>::iterator it = mPredecessor.begin() ; it != mPredecessor.end() ; ++it) {
		iLATree[it->first] = it->second;
	}
	
	// check
	int *iCheck = new int[iLANode];
	for(int i=0 ; i < iLANode ; ++i) {
		iCheck[i] = -1;
	}
	for(unsigned int i=0 ; i < vArcWord.size() ; ++i) {
		iCheck[i] = i;
	}
	for(int i=0 ; i < iLANode ; ++i) {
		if (iLATree[i] != -1) {
			iCheck[iLATree[i]] = max(iCheck[iLATree[i]],iCheck[i]);
		}
		assert(iCheck[i] != -1);
	}
	delete [] iCheck;
	
	cout << "predecessors: " << mPredecessor.size() << endl;
	cout << "arcs processed: " << iArcsProcessed << endl;
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;		
	
	cout << "-- building look-ahead tree --------------------------\n";
	cout << "# word arcs:     " << vArcWord.size() << " (words in vocabulary)\n";
	cout << "processing time: " << dTimeSeconds << " seconds\n";	
	cout << "# look-ahead nodes: " << iLANode << endl;	
	cout << "------------------------------------------------------\n";

	*iNodes = iLANode;
	return iLATree;
}

};	// end-of-namespace


