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


#include "Alignment.h"
#include "BestPath.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "HMMManager.h"
#include "HypothesisLattice.h"
#include "IOBase.h"
#include "LMFSM.h"
#include "LMManager.h"
#include "Mappings.h"
#include "NBestList.h"
#include "NBestListEntry.h"
#include "NBestListEntryElement.h"
#include "Numeric.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
HypothesisLattice::HypothesisLattice(PhoneSet *phoneSet, LexiconManager *lexiconManager)
{
	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
	
	m_lnodes = NULL;
	m_ledges = NULL;
	
	// lattice properties
	m_iNodes = -1;
	m_iEdges = -1;
	m_iFrames = -1;	
}

// destructor
HypothesisLattice::~HypothesisLattice()
{
	destroy();
}

// build the lattice container (reusing both the initial and final nodes)
void HypothesisLattice::buildContainer() {	

	return buildContainer(m_lnodeInitial,m_lnodeFinal);
}

// build the lattice container
void HypothesisLattice::buildContainer(LNode *lnodeInitial, LNode *lnodeFinal) {

	m_lnodeInitial = lnodeInitial;	
	m_lnodeFinal = lnodeFinal;	
	
	if (m_lnodes != NULL) {
		delete [] m_lnodes;
		delete [] m_ledges;
	}
	
	VLEdge vLEdge;
	VLNode vLNode;
	MLNode mLNode;
	
	vLNode.push_back(m_lnodeInitial);
	mLNode.insert(MLNode::value_type(m_lnodeInitial,true));
	
	while(vLNode.empty() == false) {
		
		LNode *lnode = vLNode.back();
		vLNode.pop_back();
		
		for(LEdge *ledge = lnode->edgeNext ; ledge != NULL ; ledge = ledge->edgePrev) {
		
			vLEdge.push_back(ledge);
			MLNode::iterator it = mLNode.find(ledge->nodeNext);
			if (it == mLNode.end()) {
				vLNode.push_back(ledge->nodeNext);
				mLNode.insert(MLNode::value_type(ledge->nodeNext,true));
			}	
		}	
	}
		
	// create an array of nodes	
	m_iNodes = (int)mLNode.size();
	m_lnodes = new LNode*[m_iNodes];
	int i=0;
	for(MLNode::iterator it = mLNode.begin() ; it != mLNode.end() ; ++it, ++i) {
		m_lnodes[i] = it->first;
		m_lnodes[i]->iNode = i;
	}

	// create an array of edges
	m_iEdges = (int)vLEdge.size();
	m_ledges = new LEdge*[m_iEdges];
	for(int i=0 ; i<m_iEdges ; ++i) {
		m_ledges[i] = vLEdge[i];
		m_ledges[i]->iEdge = i;
	}
	
	m_iFrames = m_lnodeFinal->iFrame+1;
	
	// set properties
	setProperty(LATTICE_PROPERTY_VERSION,LATTICE_FORMAT_VERSION);
	setProperty(LATTICE_PROPERTY_NODES,m_iNodes);
	setProperty(LATTICE_PROPERTY_EDGES,m_iEdges);
	setProperty(LATTICE_PROPERTY_FRAMES,m_iFrames);
}

// destroy the lattice
void HypothesisLattice::destroy() {

	if (m_ledges != NULL) {
		assert(m_lnodes != NULL);
		for(int i=0 ; i < m_iNodes ; ++i) {
			delete m_lnodes[i];
		}
		delete [] m_lnodes;
		for(int i=0 ; i < m_iEdges ; ++i) {
			deleteEdge(m_ledges[i]);
		}
		delete [] m_ledges;
	}
	
	m_mProperties.clear();
}

// print the lattice (using the lattice container)
void HypothesisLattice::print() {

	cout << "- nodes --------------------------------\n";
	for(int i=0 ; i < m_iNodes ; ++i) {
		print(m_lnodes[i]);
	}
	cout << "- edges --------------------------------\n";
	for(int i=0 ; i < m_iEdges ; ++i) {
		print(m_ledges[i]);
	}
	cout << "----------------------------------------\n";
}

// print the lattice (by traversing the nodes starting from the initial node)
void HypothesisLattice::printTraverse() {

	VLNode vLNode;
	MLNode mLNode;
	int iEdgesVisited = 0;
	
	vLNode.push_back(m_lnodeInitial);
	mLNode.insert(MLNode::value_type(m_lnodeInitial,true));
	
	cout << "-------------------------------------------\n";
	while(vLNode.empty() == false) {
		
		LNode *node = vLNode.back();
		vLNode.pop_back();
		
		for(LEdge *edge = node->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		
			print(edge);
			++iEdgesVisited;
			MLNode::iterator it = mLNode.find(edge->nodeNext);
			if (it == mLNode.end()) {
				vLNode.push_back(edge->nodeNext);
				mLNode.insert(MLNode::value_type(edge->nodeNext,true));
			}	
		}	
	}
	
	cout << "# nodes visited: " << mLNode.size() << endl;
	cout << "# edges visited: " << iEdgesVisited << endl;
	cout << "-------------------------------------------\n";
}

// print a lattice node
void HypothesisLattice::print(LEdge *edge) {

	string strLexUnit;
	m_lexiconManager->getStrLexUnitPronunciation(edge->lexUnit,strLexUnit);
	
	printf("[%4d] (%4d -> %4d) %4d %4d %20s %12.4f\n",edge->iEdge,edge->nodePrev->iNode,
		edge->nodeNext->iNode,edge->iFrameStart,edge->iFrameEnd,strLexUnit.c_str(),edge->fScoreAM);
}

// print a list of edges
void HypothesisLattice::print(LLEdge &lLEdge) {

	cout << "-------------------------------------\n";
	for(LLEdge::iterator it = lLEdge.begin() ; it != lLEdge.end() ; ++it) {
		print(*it);
	}
	cout << "-------------------------------------\n";
}


// print a lattice node
void HypothesisLattice::print(LNode *lnode) {

	printf("(%4d) %8d\n",lnode->iNode,lnode->iFrame);
	
	// print predecessors
	for(LEdge *edge = lnode->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		print(edge);
	}
}

// print the lattice properties
void HypothesisLattice::printProperties() {

	float fKB = ((float)(m_iNodes*sizeof(LNode)+m_iEdges*sizeof(LEdge)))/1024.0f;

	cout << "- lattice properties ---------------------\n";
	for(MProperty::iterator it = m_mProperties.begin() ; it != m_mProperties.end() ; ++it) {
		printf("%-16s %16s\n",it->first.c_str(),it->second.c_str());
	}
	cout << "size in memory: " << fKB << "KBs\n";
	cout << "------------------------------------------\n";
}

// store the lattice into the given file
void HypothesisLattice::store(const char *strFile, unsigned char iFormat) {

	if (iFormat == FILE_FORMAT_TEXT) {
		storeTextFormat(strFile);
	} else {
		assert(iFormat == FILE_FORMAT_BINARY);
		storeBinaryFormat(strFile);
	}
}

// store the lattice into the given file (text format)
void HypothesisLattice::storeTextFormat(const char *strFile) {

	ostringstream oss;

	FileOutput file(strFile,false);
	file.open();
	
	// print the lattice properties
	for(MProperty::iterator it = m_mProperties.begin() ; it != m_mProperties.end() ; ++it) {
		oss << setw(16) << left << it->first.c_str() << setw(16) << it->second.c_str() << endl;
	}
	
	// print the edges
	string strLexUnit;
	for(int i=0 ; i < m_iEdges ; ++i) {
		LEdge *ledge = m_ledges[i];
		m_lexiconManager->getStrLexUnitPronunciation(ledge->lexUnit,strLexUnit);
		oss << "(" << setw(4) << ledge->nodePrev->iNode << " " << setw(4) << ledge->nodeNext->iNode << ") " << setw(4) << ledge->iFrameStart << " " << setw(4) << ledge->iFrameEnd << " " << setw(4) << strLexUnit << " ";
		if (isProperty(LATTICE_PROPERTY_AM_PROB)) {
			oss << "am=" << setw(12) << ledge->fScoreAM;
		}
		if (isProperty(LATTICE_PROPERTY_LM_PROB)) {
			oss << "lm=" << setw(12) << ledge->fScoreLM;
		}
		if (isProperty(LATTICE_PROPERTY_INSERTION_PENALTY)) {
			oss << "ip=" << setw(12) << ledge->fInsertionPenalty;
		}
		if (isProperty(LATTICE_PROPERTY_FWD_PROB)) {
			oss << "fw=" << setw(12) << ledge->dScoreForward;
		}
		if (isProperty(LATTICE_PROPERTY_BWD_PROB)) {
			oss << "bw=" << setw(12) << ledge->dScoreBackward;
		}
		if (isProperty(LATTICE_PROPERTY_PP)) {
			oss << "pp=" << setw(12) << ledge->fPP;
		}
		if (isProperty(LATTICE_PROPERTY_CONFIDENCE)) {
			oss << "conf=" << setw(12) << ledge->fConfidence;
		}
		oss << endl;
		if (isProperty(LATTICE_PROPERTY_HMMS)) {
			for(int j=0 ; j < ledge->iPhones ; ++j) {
				oss << setw(10) << m_phoneSet->getStrPhone(ledge->phoneAlignment[j].iPhone);
				for(int k=0 ; k < NUMBER_HMM_STATES ; ++k) {
					if (isProperty(LATTICE_PROPERTY_PHONE_ALIGN)) {
						oss << " " << setw(6) << ledge->phoneAlignment[j].iStateBegin[k] << " " << setw(6) << ledge->phoneAlignment[j].iStateEnd[k] << " " << "[" << setw(6) << ledge->phoneAlignment[j].iHMMState[k] << "]";
					} else {
						oss << " [" << setw(6) << ledge->phoneAlignment[j].iHMMState[k] << "]" << endl;	
					}
				}
				oss << endl;
			}
		}
	}
	// print the nodes
	for(int i=0 ; i < m_iNodes ; ++i) {
		LNode *lnode = m_lnodes[i];
		oss << "(" << lnode->iNode << ") " << lnode->iFrame << endl;
	}	
	
	IOBase::writeString(file.getStream(),oss);

	file.close();
}

// store the lattice into the given file (binary format)
void HypothesisLattice::storeBinaryFormat(const char *strFile) {

	// create the file
	FileOutput file(strFile,true);
	file.open();
	
	// print the lattice properties
	// # properties
	int iProperties = (int)m_mProperties.size();
	IOBase::write(file.getStream(),iProperties);
	// actual properties
	for(MProperty::iterator it = m_mProperties.begin() ; it != m_mProperties.end() ; ++it) {
		// property
		IOBase::writeString(file.getStream(),it->first.c_str(),(unsigned int)it->first.length());
		// value
		IOBase::writeString(file.getStream(),it->second.c_str(),(unsigned int)it->second.length());
	}	
	
	// print the edges
	for(int i=0 ; i < m_iEdges ; ++i) {
		LEdge *ledge = m_ledges[i];
		IOBase::write(file.getStream(),ledge->nodePrev->iNode);
		IOBase::write(file.getStream(),ledge->nodeNext->iNode);
		IOBase::write(file.getStream(),ledge->lexUnit->iIndex);
		// write am-scores?
		if (isProperty(LATTICE_PROPERTY_AM_PROB)) {
			IOBase::write(file.getStream(),ledge->fScoreAM);
		}
		// write lm-scores?
		if (isProperty(LATTICE_PROPERTY_LM_PROB)) {
			IOBase::write(file.getStream(),ledge->fScoreLM);
		}
		// insertion penalty?
		if (isProperty(LATTICE_PROPERTY_INSERTION_PENALTY)) {
			IOBase::write(file.getStream(),ledge->fInsertionPenalty);
		}
		// write phone-level alignments?
		if (isProperty(LATTICE_PROPERTY_HMMS) || isProperty(LATTICE_PROPERTY_PHONE_ALIGN)) {
			IOBase::write(file.getStream(),ledge->iPhones);
			assert(ledge->phoneAlignment);
			for(int i=0 ; i < ledge->iPhones ; ++i) {
				IOBase::write(file.getStream(),ledge->phoneAlignment[i].iPhone);
				IOBase::write(file.getStream(),ledge->phoneAlignment[i].iPosition);
				for(int j=0 ; j < NUMBER_HMM_STATES ; ++j) {
					IOBase::write(file.getStream(),ledge->phoneAlignment[i].iStateBegin[j]);
					IOBase::write(file.getStream(),ledge->phoneAlignment[i].iStateEnd[j]);
					IOBase::write(file.getStream(),ledge->phoneAlignment[i].iHMMState[j]);
				}
			}
		}
		// write phone-accuracy?
		if (isProperty(LATTICE_PROPERTY_PHONE_ACCURACY)) {
			IOBase::write(file.getStream(),ledge->iPhones);
			IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(ledge->fPhoneAccuracy),
				sizeof(float)*ledge->iPhones);
		}
		// write best path?
		if (isProperty(LATTICE_PROPERTY_BEST_PATH)) {
			IOBase::write(file.getStream(),ledge->bBestPath);
		}
		// write pp?
		if (isProperty(LATTICE_PROPERTY_PP)) {
			IOBase::write(file.getStream(),ledge->fPP);
		}
	}
	// print the nodes
	for(int i=0 ; i < m_iNodes ; ++i) {
		LNode *lnode = m_lnodes[i];
		IOBase::write(file.getStream(),lnode->iNode);
		IOBase::write(file.getStream(),lnode->iFrame);
	}	
	
	// close the file
	file.close();
}


// load the lattice from the given file
void HypothesisLattice::load(const char *strFile) {

	loadBinaryFormat(strFile);
}

// load the lattice from the given file (binary format)
void HypothesisLattice::loadBinaryFormat(const char *strFile) {
		
	// open the file
	FileInput file(strFile,true);
	file.open();
	
	// load the lattice properties
	
	// # properties
	int iProperties = -1;
	
	IOBase::read(file.getStream(),&iProperties);
	char *strProperty = NULL;
	char *strValue = NULL;
	for(int i=0 ; i < iProperties ; ++i) {
		// property 
		IOBase::readString(file.getStream(),&strProperty);
		// value 
		IOBase::readString(file.getStream(),&strValue);
		// check that is not multiple defined
		if (isProperty(strProperty)) {
			BVC_ERROR << "lattice property \"" << strProperty <<  "\" defined multiple times";
		}
		setProperty(strProperty,strValue);
		delete [] strProperty;
		delete [] strValue;
	}
	
	// make sure mandatory properties are defined
	if (!isProperty(LATTICE_PROPERTY_VERSION) ||
		!isProperty(LATTICE_PROPERTY_NODES) || 
		!isProperty(LATTICE_PROPERTY_EDGES) ||	
		!isProperty(LATTICE_PROPERTY_FRAMES)) {
		BVC_ERROR << "wrong lattice format: missing properties";
	}
	
	// get # nodes and edges
	m_iNodes = atoi(getPropertyValue(LATTICE_PROPERTY_NODES));
	m_iEdges = atoi(getPropertyValue(LATTICE_PROPERTY_EDGES));
	m_iFrames = atoi(getPropertyValue(LATTICE_PROPERTY_FRAMES));

	// allocate memory for the edges
	m_ledges = new LEdge*[m_iEdges];
	int *iNodePrev = new int[m_iEdges];
	int *iNodeNext = new int[m_iEdges];
	int iLexUnitIndex = -1;
	// read the edges	
	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i] = new LEdge;
		LEdge *ledge = m_ledges[i];
		ledge->iEdge = i;
		IOBase::read(file.getStream(),&iNodePrev[i]);
		IOBase::read(file.getStream(),&iNodeNext[i]);
		ledge->iFrameStart = -1;
		ledge->iFrameEnd = -1;
		IOBase::read(file.getStream(),&iLexUnitIndex);
		ledge->lexUnit = m_lexiconManager->getLexUnitByIndex(iLexUnitIndex);
		if (!ledge->lexUnit) {
			BVC_ERROR << "wrong lattice format: unknown lexical unit, index: " << iLexUnitIndex;
		}
		assert(ledge->lexUnit != m_lexiconManager->m_lexUnitBegSentence);
		
		// context
		ledge->iContextLeft = NULL;
		ledge->iContextRight = NULL;
		
		// load am-scores?
		ledge->fScoreAM = -FLT_MAX;
		if (isProperty(LATTICE_PROPERTY_AM_PROB)) {
			IOBase::read(file.getStream(),&ledge->fScoreAM);
		}
		// language model score?
		ledge->fScoreLM = -FLT_MAX;
		if (isProperty(LATTICE_PROPERTY_LM_PROB)) {
			IOBase::read(file.getStream(),&ledge->fScoreLM);
		}
		// insertion penalty?
		ledge->fInsertionPenalty = -FLT_MAX;
		if (isProperty(LATTICE_PROPERTY_INSERTION_PENALTY)) {
			IOBase::read(file.getStream(),&ledge->fInsertionPenalty);
		}	
		// load phone-level alignments?
		ledge->iPhones = -1;
		ledge->phoneAlignment = NULL;
		if (isProperty(LATTICE_PROPERTY_HMMS) || isProperty(LATTICE_PROPERTY_PHONE_ALIGN)) {
			IOBase::read(file.getStream(),&ledge->iPhones);
			if (ledge->iPhones != (int)ledge->lexUnit->vPhones.size()) {
				BVC_ERROR << "wrong lattice format: inconsistent number of phones in edge";
			}	
			ledge->phoneAlignment = new LPhoneAlignment[ledge->iPhones];
			for(int i=0 ; i < ledge->iPhones ; ++i) {
				IOBase::read(file.getStream(),&ledge->phoneAlignment[i].iPhone);
				IOBase::read(file.getStream(),&ledge->phoneAlignment[i].iPosition);
				for(int j=0 ; j < NUMBER_HMM_STATES ; ++j) {
					IOBase::read(file.getStream(),&ledge->phoneAlignment[i].iStateBegin[j]);
					IOBase::read(file.getStream(),&ledge->phoneAlignment[i].iStateEnd[j]);
					IOBase::read(file.getStream(),&ledge->phoneAlignment[i].iHMMState[j]);
				}
			}	
		}
		// load phone accuracy
		ledge->fPhoneAccuracy = NULL;
		if (isProperty(LATTICE_PROPERTY_PHONE_ACCURACY)) {
			IOBase::read(file.getStream(),&ledge->iPhones);
			if (ledge->iPhones != (int)ledge->lexUnit->vPhones.size()) {
				BVC_ERROR << "wrong lattice format: inconsistent number of phones in edge";
			}
			ledge->fPhoneAccuracy = new float[ledge->iPhones];
			IOBase::readBytes(file.getStream(),(char*)ledge->fPhoneAccuracy,sizeof(float)*ledge->iPhones);
		}
		// load best path?
		ledge->bBestPath = false;
		if (isProperty(LATTICE_PROPERTY_BEST_PATH)) {
			IOBase::read(file.getStream(),&ledge->bBestPath);	
		}
		// load pp?
		ledge->fPP = -FLT_MAX;	
		if (isProperty(LATTICE_PROPERTY_PP)) {
			IOBase::read(file.getStream(),&ledge->fPP);
		}
	}
	
	// allocate memory for the nodes
	m_lnodes = new LNode*[m_iNodes];
	// read the nodes
	for(int i=0 ; i < m_iNodes ; ++i) {
		m_lnodes[i] = new LNode;
		LNode *lnode = m_lnodes[i];
		lnode->iNode = i;
		IOBase::read(file.getStream(),&lnode->iNode);
		IOBase::read(file.getStream(),&lnode->iFrame);
		if ((lnode->iFrame >= m_iFrames) || (lnode->iFrame < -1)) {
			BVC_ERROR << "wrong lattice format: inconsistent node time: " << lnode->iFrame;
		}
		lnode->edgePrev = NULL;
		lnode->edgeNext = NULL;
	}
	// connect edges and nodes
	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->nodePrev = m_lnodes[iNodePrev[i]];
		m_ledges[i]->nodeNext = m_lnodes[iNodeNext[i]];
		m_ledges[i]->edgeNext = m_ledges[i]->nodeNext->edgePrev;
		m_ledges[i]->nodeNext->edgePrev = m_ledges[i];
		m_ledges[i]->edgePrev = m_ledges[i]->nodePrev->edgeNext;
		m_ledges[i]->nodePrev->edgeNext = m_ledges[i];
		m_ledges[i]->iFrameStart = m_ledges[i]->nodePrev->iFrame+1;;
		m_ledges[i]->iFrameEnd = m_ledges[i]->nodeNext->iFrame;	
		assert(m_ledges[i]->iFrameEnd > m_ledges[i]->iFrameStart);
	}
	delete [] iNodePrev;
	delete [] iNodeNext;	
	
	// get the initial and final nodes 
	m_lnodeInitial = NULL;
	m_lnodeFinal = NULL;
	for(int i=0 ; i < m_iNodes ; ++i) {
		if (m_lnodes[i]->edgePrev == NULL) {
			if (m_lnodeInitial != NULL) {
				BVC_ERROR << "wrong lattice format: multiple initial nodes defined";
			}	
			m_lnodeInitial = m_lnodes[i];
		}
		if (m_lnodes[i]->edgeNext == NULL) {
			if (m_lnodeFinal != NULL) {
				BVC_ERROR << "wrong lattice format: multiple final nodes defined";
			}	
			m_lnodeFinal = m_lnodes[i];
		}	
	}
	if ((m_lnodeInitial == NULL) || (m_lnodeFinal == NULL)) {
		BVC_ERROR << "wrong lattice format: no initial/final node";
	}
	if ((m_lnodeInitial->iFrame != -1) || (m_lnodeFinal->iFrame != m_iFrames-1)) {
		BVC_ERROR << "wrong lattice format: inconsistent initial/final nodes";
	}
	if (m_lnodeInitial == m_lnodeFinal) {
		BVC_ERROR << "wrong lattice format: initial node and final node are the same";
	}
	
	// close the file
	file.close();
}




// check lattice correctness
// it performs the following checks:
// - time alignment is consistent
// - lattice is connected
void HypothesisLattice::check() {

	// (1) check that the lattice is connected (left to right)

	// mark all the nodes as not-visited
	for(int i=0 ; i<m_iNodes ; ++i) {
		m_lnodes[i]->bTouched = false;
	}
	
	VLNode vLNode;
	m_lnodeInitial->bTouched = true;
	vLNode.push_back(m_lnodeInitial);
	while(vLNode.empty() == false) {
		
		LNode *lnode = vLNode.back();
		vLNode.pop_back();
		
		for(LEdge *ledge = lnode->edgeNext ; ledge != NULL ; ledge = ledge->edgePrev) {
			// if not visited mark it as visited and put it in the vector	
			if (ledge->nodeNext->bTouched == false) {
				ledge->nodeNext->bTouched = true;
				vLNode.push_back(ledge->nodeNext);
			}	
		}
	}
	
	// check that all the nodes are visited
	for(int i=0 ; i<m_iNodes ; ++i) {
		if (m_lnodes[i]->bTouched == false) {
			BVC_ERROR << "lattice did not pass check";
		}
	}

	// (1) check that the lattice is connected (right to left)
	assert(vLNode.empty());

	// mark all the nodes as not-visited
	for(int i=0 ; i<m_iNodes ; ++i) {
		m_lnodes[i]->bTouched = false;
	}
	
	m_lnodeFinal->bTouched = true;
	vLNode.push_back(m_lnodeFinal);
	while(vLNode.empty() == false) {
		
		LNode *lnode = vLNode.back();
		vLNode.pop_back();
		
		for(LEdge *ledge = lnode->edgePrev ; ledge != NULL ; ledge = ledge->edgeNext) {
			// if not visited mark it as visited and put it in the vector	
			if (ledge->nodePrev->bTouched == false) {
				ledge->nodePrev->bTouched = true;
				vLNode.push_back(ledge->nodePrev);
			}	
		}
	}
	
	// check that all the nodes are visited
	for(int i=0 ; i<m_iNodes ; ++i) {
		if (m_lnodes[i]->bTouched == false) {
			BVC_ERROR << "lattice did not pass check";
		}
	}
}

// traverse the lattice forward and merge identical edges
// - it does not add/remove any paths (word-sequence) to the original lattice
void HypothesisLattice::forwardEdgeMerge() {
	
	int iNodesProcessed = 0;
	int iEdgesRemoved = 0;
	int iNodesRemoved = 0;
	
	// allocate a structure to keep the node states	
	unsigned char *iNodeState = new unsigned char[m_iNodes];
	for(int i=0 ; i < m_iNodes ; ++i) {
		iNodeState[i] = LATTICE_NODE_STATE_UNSEEN;
	}
	
	// start from the initial node
	LEdge *edgeAux = NULL;
	LEdge **edgePointer = NULL;
	LLNode lNodes;
	lNodes.push_back(m_lnodeInitial);
	iNodeState[m_lnodeInitial->iNode] = LATTICE_NODE_STATE_QUEUED;
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		LNode *node = lNodes.front();
		lNodes.pop_front();
		
		// delete the node if marked as deleted
		if (iNodeState[node->iNode] == LATTICE_NODE_STATE_DELETED) {
			delete node;
			continue;
		}		
		
		iNodeState[node->iNode] = LATTICE_NODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LEdge *edge1 = node->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
			LEdge *edge2Prev = NULL;
			for(LEdge *edge2 = edge1->edgePrev ; edge2 != NULL ; edge2 = edge2Prev) {
				edge2Prev = edge2->edgePrev;
				if (equivalentEdges(edge1,edge2) && equivalentPredecessors(edge1->nodeNext,edge2->nodeNext)) {
					
					if (edge1->nodeNext == edge2->nodeNext) {
						assert((edge1->nodeNext == m_lnodeFinal) && (edge2->nodeNext == m_lnodeFinal));
						// this case only happens for unigram-lattices and it is rare, should be handled separately
						continue;
					}	
					
					// (2) move predecessors from one edge to the other (if not already there)
					LEdge *edgeAdd = NULL;
					for(LEdge *edgeA = edge2->nodeNext->edgeNext ; edgeA != NULL ; ) {
						bool bFound = false;
						LEdge *edgeB = NULL;
						for(edgeB = edge1->nodeNext->edgeNext ; edgeB != NULL ; edgeB = edgeB->edgePrev) {
							if (equivalentEdges(edgeA,edgeB) && (edgeA->nodeNext == edgeB->nodeNext)) {
								bFound = true;
								break;
							}
						}
						if (bFound == false) {
							edgeAux = edgeA->edgePrev;
							edgeA->nodePrev = edge1->nodeNext;
							edgeA->edgePrev = edgeAdd;
							edgeAdd = edgeA;
							edgeA = edgeAux;
						} else {
							// extract the redundant edge
							edgePointer = &edgeA->nodeNext->edgePrev; 
							while(*edgePointer != edgeA) {
								edgePointer = &(*edgePointer)->edgeNext;
							}
							*edgePointer = edgeA->edgeNext;	
							++iEdgesRemoved;
							// move to the next predecessor
							LEdge *edgeDelete = edgeA;
							edgeA = edgeA->edgePrev;
							deleteEdge(edgeDelete);
						}
					}
					edgePointer = &edge1->nodeNext->edgeNext; 
					while(*edgePointer != NULL) {
						edgePointer = &(*edgePointer)->edgePrev;
					}
					(*edgePointer) = edgeAdd;
					// (3) extract all the successor edges
					LNode *nodeNext = edge2->nodeNext;
					for(LEdge *edge = nodeNext->edgePrev ; edge != NULL ; ) {
						edgePointer = &edge->nodePrev->edgeNext; 
						while(*edgePointer != edge) {
							assert(*edgePointer != NULL);
							edgePointer = &(*edgePointer)->edgePrev;
						}
						LEdge *edgeDelete = edge;
						*edgePointer = edge->edgePrev;
						edge = edge->edgeNext;
						if (edgeDelete == edge2Prev) {
							edge2Prev = edgeDelete->edgePrev;
						}						
						deleteEdge(edgeDelete);
						++iEdgesRemoved;	
					}
					++iNodesRemoved;
					++iEdgesRemoved;
					// delete the node now or mark it to be deleted later
					if (iNodeState[nodeNext->iNode] == LATTICE_NODE_STATE_QUEUED) {	
						iNodeState[nodeNext->iNode] = LATTICE_NODE_STATE_DELETED;
					} else {
						//printf("deleted(1): %x\n",nodeNext);
						delete nodeNext;
					}
				}
			}
		}	
		
		// put in the queue those predecessor nodes that have all their successors processed
		for(LEdge *edge = node->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (iNodeState[edge->nodeNext->iNode] != LATTICE_NODE_STATE_UNSEEN) {
				continue;
			}
			bool bReady = true;
			for(LEdge *edge1 = edge->nodeNext->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
				if (iNodeState[edge1->nodePrev->iNode] != LATTICE_NODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}	
			if (bReady) {
				lNodes.push_back(edge->nodeNext);
				iNodeState[edge->nodeNext->iNode] = LATTICE_NODE_STATE_QUEUED;
			}
		}	
	}
	
	delete [] iNodeState;
	
	// rebuild the lattice container
	buildContainer(m_lnodeInitial,m_lnodeFinal);
	
	check();
		
	float fCompressionRatioNodes = ((float)(100*iNodesRemoved))/((float)m_iNodes);
	float fCompressionRatioEdges = ((float)(100*iEdgesRemoved))/((float)m_iEdges);
	BVC_VERB << "compression ratio:";
	BVC_VERB << "# nodes removed: " << iNodesRemoved <<" (" << fCompressionRatioNodes << "%)";
	BVC_VERB << "# edges removed: " << iEdgesRemoved <<" (" << fCompressionRatioEdges << "%)";	
}

// traverse the lattice backwards and merge identical edges
// - it does not add/remove any paths (word-sequence) to the original lattice
void HypothesisLattice::backwardEdgeMerge() {
	
	int iNodesProcessed = 0;
	int iEdgesRemoved = 0;
	int iNodesRemoved = 0;
	
	// allocate a structure to keep the node states	
	unsigned char *iNodeState = new unsigned char[m_iNodes];
	for(int i=0 ; i < m_iNodes ; ++i) {
		iNodeState[i] = LATTICE_NODE_STATE_UNSEEN;
	}
	
	// start from the final node
	LEdge *edgeAux = NULL;
	LEdge **edgePointer = NULL;
	LLNode lNodes;
	lNodes.push_back(m_lnodeFinal);
	iNodeState[m_lnodeFinal->iNode] = LATTICE_NODE_STATE_QUEUED;
		
	//map<LEdge*,bool> mEdgeDeleted;
	
	// process all nodes
	while(lNodes.empty() == false) {
	
		// get a node to process
		LNode *node = lNodes.front();
		lNodes.pop_front();
		
		assert(iNodeState[node->iNode] == LATTICE_NODE_STATE_QUEUED);
		
		//printf("processing: %d %x\n",node->iNode,node);
		
		// delete the node if marked as deleted
		if (iNodeState[node->iNode] == LATTICE_NODE_STATE_DELETED) {
			//printf("deleted(2): %x\n",node);
			delete node;
			continue;
		}
		
		iNodeState[node->iNode] = LATTICE_NODE_STATE_PROCESSED;
		++iNodesProcessed;
	
		// look for equivalent arcs and merge them
		for(LEdge *edge1 = node->edgePrev ; edge1 != NULL ; edge1 = edge1->edgeNext) {
			LEdge *edge2Next = NULL;
			for(LEdge *edge2 = edge1->edgeNext ; edge2 != NULL ; edge2 = edge2Next) {
				edge2Next = edge2->edgeNext;
				//assert(mEdgeDeleted.find(edge2) == mEdgeDeleted.end());
				//assert(mEdgeDeleted.find(edge1) == mEdgeDeleted.end());
				if (equivalentEdges(edge1,edge2) && equivalentSuccessors(edge1->nodePrev,edge2->nodePrev)) {	
					
					if (edge1->nodePrev == edge2->nodePrev) {
						// this case only happens for unigram-lattices and it is rare, should be handled separately
						continue;
					}	
					
					// (2) move predecessors from one edge to the other (if not already there)
					LEdge *edgeAdd = NULL;
					for(LEdge *edgeA = edge2->nodePrev->edgePrev ; edgeA != NULL ; ) {
						bool bFound = false;
						LEdge *edgeB = NULL;
						for(edgeB = edge1->nodePrev->edgePrev ; edgeB != NULL ; edgeB = edgeB->edgeNext) {
							if (equivalentEdges(edgeA,edgeB) && (edgeA->nodePrev == edgeB->nodePrev)) {
								bFound = true;
								break;
							}
						}
						if (bFound == false) {
							edgeAux = edgeA->edgeNext;
							edgeA->nodeNext = edge1->nodePrev;
							edgeA->edgeNext = edgeAdd;
							edgeAdd = edgeA;
							edgeA = edgeAux;
						} else {
							// extract the redundant edge
							edgePointer = &edgeA->nodePrev->edgeNext; 
							while(*edgePointer != edgeA) {
								edgePointer = &(*edgePointer)->edgePrev;
							}
							*edgePointer = edgeA->edgePrev;	
							++iEdgesRemoved;
							// move to the next predecessor
							LEdge *edgeDelete = edgeA;		
							edgeA = edgeA->edgeNext;
							//mEdgeDeleted.insert(map<LEdge*,bool>::value_type(edgeDelete,true));
							deleteEdge(edgeDelete);
						}
					}
					edgePointer = &edge1->nodePrev->edgePrev; 
					while(*edgePointer != NULL) {
						edgePointer = &(*edgePointer)->edgeNext;
					}
					(*edgePointer) = edgeAdd;
					// (3) extract all the successor edges
					LNode *nodePrev = edge2->nodePrev;
					for(LEdge *edge = nodePrev->edgeNext ; edge != NULL ; ) {
						edgePointer = &edge->nodeNext->edgePrev; 
						while(*edgePointer != edge) {
							assert(*edgePointer != NULL);
							edgePointer = &(*edgePointer)->edgeNext;
						}
						LEdge *edgeDelete = edge;
						*edgePointer = edge->edgeNext;
						edge = edge->edgePrev;
						//mEdgeDeleted.insert(map<LEdge*,bool>::value_type(edgeDelete,true));
						if (edgeDelete == edge2Next) {
							edge2Next = edgeDelete->edgeNext;
						}
						deleteEdge(edgeDelete);
						++iEdgesRemoved;	
					}
					++iNodesRemoved;
					// delete the node now or mark it to be deleted later
					if (iNodeState[nodePrev->iNode] == LATTICE_NODE_STATE_QUEUED) {	
						iNodeState[nodePrev->iNode] = LATTICE_NODE_STATE_DELETED;
					} else {
						//printf("node deleted(1): %d %x\n",nodePrev->iNode,nodePrev);
						delete nodePrev;						
					}
				}
			}
		}
		
		// put in the queue those predecessor nodes that have all their successors processed
		for(LEdge *edge = node->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
			if (iNodeState[edge->nodePrev->iNode] != LATTICE_NODE_STATE_UNSEEN) {
				continue;
			}
			bool bReady = true;
			for(LEdge *edge1 = edge->nodePrev->edgeNext ; edge1 != NULL ; edge1 = edge1->edgePrev) {
				if (iNodeState[edge1->nodeNext->iNode] != LATTICE_NODE_STATE_PROCESSED) {
					bReady = false;
					break;
				}
			}	
			if (bReady) {
				//printf("node queued: %d (%x)\n",edge->nodePrev->iNode,edge->nodePrev);
				lNodes.push_back(edge->nodePrev);
				iNodeState[edge->nodePrev->iNode] = LATTICE_NODE_STATE_QUEUED;
			}
		}	
	}
	
	delete [] iNodeState;
	
	// rebuild the lattice container
	buildContainer(m_lnodeInitial,m_lnodeFinal);
	
	check();
		
	float fCompressionRatioNodes = ((float)(100*iNodesRemoved))/((float)m_iNodes);
	float fCompressionRatioEdges = ((float)(100*iEdgesRemoved))/((float)m_iEdges);
	BVC_VERB << "compression ratio:";
	BVC_VERB << "# nodes removed: " << iNodesRemoved <<" (" << fCompressionRatioNodes << "%)";
	BVC_VERB << "# edges removed: " << iEdgesRemoved <<" (" << fCompressionRatioEdges << "%)";	
}

// compute forward backward scores
void HypothesisLattice::computeForwardBackwardScores(float fScalingAM, float fScalingLM) {

	// check lattice properties
	if (!isProperty(LATTICE_PROPERTY_AM_PROB) || 
		!isProperty(LATTICE_PROPERTY_LM_PROB) ||
		!isProperty(LATTICE_PROPERTY_INSERTION_PENALTY)) {
		BVC_ERROR << "acoustic likelihood, lm likelihood and insertion penalty are needed to compute fwd-bwd scores!";
	}
			
	// keep scaling factors
	m_fAMScalingFactor = fScalingAM;
	m_fLMScalingFactor = fScalingLM;

	// compute forward scores
	untouchEdges();
	assert(m_lnodeFinal->edgePrev != NULL);
	for(LEdge *edge = m_lnodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		forward(edge);
	}
	
	// compute backward scores
	untouchEdges();
	assert(m_lnodeInitial->edgeNext != NULL);
	for(LEdge *edge = m_lnodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		backward(edge);
	}
	
	// update the lattice properties
	setProperty(LATTICE_PROPERTY_FWD_PROB,"yes");
	setProperty(LATTICE_PROPERTY_BWD_PROB,"yes");
}

// compute forward score for the edge
void HypothesisLattice::forward(LEdge *edge) {

	if (edge->nodePrev == m_lnodeInitial) {
		edge->bTouched = true;
		edge->dScoreForward = (edge->fScoreAM+edge->fInsertionPenalty)*m_fAMScalingFactor + edge->fScoreLM*m_fLMScalingFactor;
		return;
   }
   
	// (1) accumulate the forward score from predecessors
	double dScoreAcc = -DBL_MAX;
	double dAux = 0.0;
	LEdge *edgePred = edge->nodePrev->edgePrev;
	while(edgePred != NULL) {
		
		// compute the forward score of the predecessor if not available
		if (edgePred->bTouched == false) {
			forward(edgePred);	
		} 
		assert(edgePred->bTouched);
		assert(edgePred->dScoreForward != -DBL_MAX);	
		
		// accumulate forward probabilities
		dAux = edgePred->dScoreForward + edge->fScoreLM*m_fLMScalingFactor;
		dScoreAcc = Numeric::logAddition(dScoreAcc,dAux);
		
		edgePred = edgePred->edgeNext;
	}
 
	// (2) compute the forward score of the edge (the accumulated score)
	edge->dScoreForward = dScoreAcc + (edge->fScoreAM+edge->fInsertionPenalty)*m_fAMScalingFactor;
	edge->bTouched = true;	
}

// compute backward score for the node
void HypothesisLattice::backward(LEdge *edge) {

	if (edge->nodeNext == m_lnodeFinal) {
		edge->dScoreBackward = 0.0;
		edge->bTouched = true;
		return;
   }
   
	// (1) compute the backward score of all the successors (only those that are left)
	double dScoreAcc = -DBL_MAX;		// accumulated score from successors
	double dAux = 0.0;
	LEdge *edgeSucc = edge->nodeNext->edgeNext;
	while(edgeSucc != NULL) {
		
		// compute the backward score of the successor if not available
		if (edgeSucc->bTouched == false) {
			backward(edgeSucc);	
		}
		assert(edgeSucc->bTouched);
		assert(edgeSucc->dScoreBackward != -DBL_MAX);
		
		// accumulate backward probabilities
		dAux = edgeSucc->dScoreBackward + edgeSucc->fScoreLM*m_fLMScalingFactor +
			(edgeSucc->fScoreAM+edgeSucc->fInsertionPenalty)*m_fAMScalingFactor;
		dScoreAcc = Numeric::logAddition(dScoreAcc,dAux);
		
		edgeSucc = edgeSucc->edgePrev;
	}
 
	// (2) compute the backward score of the edge (the accumulated)
	edge->dScoreBackward = dScoreAcc;
	edge->bTouched = true;
}

// compute posterior probabilities from the hypothesis graph
void HypothesisLattice::computePosteriorProbabilities() {

	// make sure forward/backward scores are available
	if ((!isProperty(LATTICE_PROPERTY_FWD_PROB)) || (!isProperty(LATTICE_PROPERTY_BWD_PROB))) {
		BVC_ERROR << "forward and backward log-likelihoods are needed to compute posterior probabilities!";
	}

	// (1) compute normalization factor
	assert(m_lnodeFinal->edgePrev != NULL);
	double dNorm = -DBL_MAX;
	double dAccIni = -DBL_MAX;
	double dAccFin = -DBL_MAX;
	for(LEdge *edge = m_lnodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		dNorm = Numeric::logAddition(dNorm,edge->dScoreForward+edge->dScoreBackward);
		dAccFin = Numeric::logAddition(dAccFin,edge->dScoreForward+edge->dScoreBackward);
	}
	for(LEdge *edge = m_lnodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		dAccIni = Numeric::logAddition(dAccIni,edge->dScoreForward+edge->dScoreBackward);
	}
	
	// (2) compute posterior probabilities
	for(int i=0 ; i < m_iEdges ; ++i) {
		// note that the acoustic score is in both the fwd and bwd probabilities
		m_ledges[i]->fPP = (float)exp(m_ledges[i]->dScoreForward+m_ledges[i]->dScoreBackward-dNorm);
	}	
	
	setProperty(LATTICE_PROPERTY_PP,"yes");
}

// mark all the nodes as not touched
void HypothesisLattice::untouchNodes() {

	for(int i=0 ; i < m_iNodes ; ++i) {
		m_lnodes[i]->bTouched = false;
	}
}

// mark all the edges as not touched
void HypothesisLattice::untouchEdges() {

	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->bTouched = false;
	}
}

// reset all the auxiliar edges
void HypothesisLattice::resetAuxEdges() {

	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->edgeNextAux = NULL;
		m_ledges[i]->edgePrevAux = NULL;
	}
}



// compute posterior-probability based confidence estimates
void HypothesisLattice::computeConfidenceScore(unsigned char iConfidenceMeasure) {

	if (!isProperty(LATTICE_PROPERTY_PP)) {
		BVC_ERROR << "posterior probabilities are needed to compute confidence scores!";
	}	

	// raw posteriors
	if (iConfidenceMeasure == CONFIDENCE_MEASURE_POSTERIORS) {
	
		for(int i=0 ; i < m_iEdges ; ++i) {
			m_ledges[i]->fConfidence = m_ledges[i]->fPP;
		}
	}
	// accumulate posteriors across overlapping word instances 
	// note: confidence value is in the range [0.0 - inf] since posteriors are accumulated
	// from overlapping edges in the lattice that may not overlap between them
	else if (iConfidenceMeasure == CONFIDENCE_MEASURE_ACCUMULATED) {
	
		// (1) arrange the edges by lexical unit
		MIntLEdge mLexUnitEdges;
		for(int i=0 ; i < m_iEdges ; ++i) {
			m_ledges[i]->fConfidence = 0.0;
			MIntLEdge::iterator it = mLexUnitEdges.find(m_ledges[i]->lexUnit->iLexUnit);
			if (it != mLexUnitEdges.end()) {
				m_ledges[i]->edgeNextAux = it->second;
				it->second = m_ledges[i];
			}
			else {
				m_ledges[i]->edgeNextAux = NULL;
				mLexUnitEdges.insert(MIntLEdge::value_type(m_ledges[i]->lexUnit->iLexUnit,m_ledges[i]));	
			}
		}
		
		// (2) accumulate posteriors across overlapping lexical units
		for(MIntLEdge::iterator it = mLexUnitEdges.begin() ; it != mLexUnitEdges.end() ; ++it) {
			LEdge *edge = it->second;
			do {
				edge->fConfidence += edge->fPP;
				LEdge *edgeAux = edge->edgeNextAux;
				while(edgeAux != NULL) {
					assert(edge != edgeAux);
					if (overlap(edge,edgeAux)) {	
						edgeAux->fConfidence += edge->fPP;
						edge->fConfidence += edgeAux->fPP;
					}
					edgeAux = edgeAux->edgeNextAux; 
				}
				edge = edge->edgeNextAux;
			} while(edge != NULL);
		}	

		mLexUnitEdges.clear();
	}	
	// maximum posterior probability of an overlapping edge
	else {	
		assert(iConfidenceMeasure == CONFIDENCE_MEASURE_MAXIMUM);
		
		// (1) arrange the edges by lexical unit
		MIntLEdge mLexUnitEdges;
		for(int i=0 ; i < m_iEdges ; ++i) {
			m_ledges[i]->fConfidence = 0.0;
			MIntLEdge::iterator it = mLexUnitEdges.find(m_ledges[i]->lexUnit->iLexUnit);
			if (it != mLexUnitEdges.end()) {
				m_ledges[i]->edgeNextAux = it->second;
				it->second = m_ledges[i];
			}
			else {
				m_ledges[i]->edgeNextAux = NULL;
				mLexUnitEdges.insert(MIntLEdge::value_type(m_ledges[i]->lexUnit->iLexUnit,m_ledges[i]));	
			}
		}
				
		// (2) accumulate posteriors across all time-frames and keep the max
		float *fFrameScoreMax = new float[m_iFrames];
		for(MIntLEdge::iterator it = mLexUnitEdges.begin() ; it != mLexUnitEdges.end() ; ++it) {
			LEdge *edge = it->second;
			// accumulate posteriors at the frame level
			do {
				// initialize to the posterior probability
				for(int i=0 ; i <= (edge->iFrameEnd-edge->iFrameStart) ; ++i) {
					fFrameScoreMax[i] = edge->fPP;
				}
				LEdge *edgeAux = it->second;
				while(edgeAux != NULL) {
					// two nodes have only frames in common if they overlap somewhere
					if ((edge != edgeAux) && overlap(edge,edgeAux)) {	
						for(int i=edge->iFrameStart ; i <= edge->iFrameEnd ; ++i) {
							if (spanTimeFrame(edgeAux,i)) {
								fFrameScoreMax[i-edge->iFrameStart] += edgeAux->fPP; 
							}
						}
					}
					edgeAux = edgeAux->edgeNextAux; 
				}
				// keep the max
				float fMax = edge->fPP;		// this is the floor of the maximum value
				for(int i=0 ; i <= (edge->iFrameEnd-edge->iFrameStart) ; ++i) {
					if (fFrameScoreMax[i] > fMax) {
						fMax = fFrameScoreMax[i];
					}
				}	
				edge->fConfidence = fMax;
			
				edge = edge->edgeNextAux;
			} while(edge != NULL);
		}
		delete [] fFrameScoreMax;
		mLexUnitEdges.clear();	
	}
	
	// update the lattice properties
	setProperty(LATTICE_PROPERTY_CONFIDENCE,getStrConfidenceMeasure(iConfidenceMeasure));
}


// set all the acoustic-scores to the given value
void HypothesisLattice::setAMScores(float fValue) {

	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->fScoreAM = fValue;
	}
}

// set all the language-model-scores to the given value
void HypothesisLattice::setLMScores(float fValue) {

	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->fScoreLM = fValue;
	}
}

// compute the lattice depth given the reference text
// -> lattice depth is defined as the average number of lattice arcs containing words
// that cross every time frame
// -> lattice density is defined as the number of arcs in the lattice
// divided by the number of words in the reference
LatticeDepth *HypothesisLattice::computeDepth(bool bIgnoreNonWords) {
	
	LatticeDepth *latticeDepth = new LatticeDepth;
	
	latticeDepth->iFrames = m_iFrames;
	latticeDepth->iFramesEdges = 0;

	for(int i=0 ; i<m_iEdges ; ++i) {
		assert(m_ledges[i] != NULL);
		LexUnit *lexUnit = m_ledges[i]->lexUnit;
		if (bIgnoreNonWords == false) {
			latticeDepth->iFramesEdges += (m_ledges[i]->iFrameEnd-m_ledges[i]->iFrameStart+1);
		} else if (m_lexiconManager->isStandard(lexUnit)) {
			latticeDepth->iFramesEdges += (m_ledges[i]->iFrameEnd-m_ledges[i]->iFrameStart+1);
		}	
	}
	
	latticeDepth->fDepth = ((float)latticeDepth->iFramesEdges)/((float)latticeDepth->iFrames);
	
	return latticeDepth;
}

// -------------------------------------------------------------------------------------------
// HMM-id marking 
// -------------------------------------------------------------------------------------------

// mark each of the edges in the lattice with HMM-information, the lattice
// will be expanded as needed in order to accomodate different word-contexts
void HypothesisLattice::hmmMarking(HMMManager *hmmManager) {

	m_iContextSizeWW = hmmManager->getContextSizeHMM();
	m_iContextSizeCW = hmmManager->getContextSizeHMMCW();
	
	//m_iContextSizeWW = m_iContextSizeCW = 3;
	m_iPhoneContextPadding = m_phoneSet->size();
	int iContextSizeMax = max(m_iContextSizeCW,m_iContextSizeWW);
	
	assert(m_iContextSizeCW > 0);
	assert(m_iContextSizeWW > 0);
	
	int iNodesOriginal = m_iNodes;
	int iEdgesOriginal = m_iEdges;
	int iEdgesAdded = 0;
	int iNodesAdded = 0;	
	
	// initialize contexts
	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->iContextLeft = NULL;
		m_ledges[i]->iContextRight = NULL;
	}
	
	LLEdge lEdge;
	
	// (1) left context propagation
	
 	// create a fake edge to make the algorithm simpler
	LEdge *edgeFake = new LEdge;
	edgeFake->nodeNext = m_lnodeInitial;
	if (iContextSizeMax > 1) {
		edgeFake->iContextLeft = new unsigned char[iContextSizeMax];
		for(int i=0 ; i < iContextSizeMax ; ++i) {
			edgeFake->iContextLeft[i] = m_phoneSet->getPhoneIndexSilence();
		}			
	}
	edgeFake->lexUnit = m_lexiconManager->getLexUnitSilence();
	lEdge.push_back(edgeFake);
	
	// process the edges one by one
	while(lEdge.empty() == false) {
	
		LEdge *edgeAux = lEdge.front();
		lEdge.pop_front();	
		
		// does the edge go to the final node?
		if (edgeAux->nodeNext->edgeNext == NULL) {
			assert(edgeAux->nodeNext == m_lnodeFinal);
			continue;	
		}
		
		// build the left context that will be propagated from the edge (cross-word)
		edgeAux->iPhones = (int)edgeAux->lexUnit->vPhones.size();
		unsigned char *iContextLeftProp = newBlankContext();
		for(int i=0 ; i < min((int)m_iContextSizeCW,(int)(edgeAux->iPhones)) ; ++i) {
			iContextLeftProp[m_iContextSizeCW-i-1] = edgeAux->lexUnit->vPhones[edgeAux->iPhones-1-i];
		}
		
		// successor edges do not have a context yet
		bool bContextLeft = false;
		for(LEdge *edge = edgeAux->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
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
			for(LEdge *edge = edgeAux->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
				edge->iContextLeft = copyContext(iContextLeftProp);		
				lEdge.push_back(edge);
			}
		} 
		// successor edges already have a context
		else {
			// different context: duplicate node and edges and put the edges in the queue
			if (memcmp(iContextLeftProp,edgeAux->nodeNext->edgeNext->iContextLeft,
				iContextSizeMax*sizeof(unsigned char)) != 0) {	
				// extract the auxiliar edge from the original destination node
				LEdge **edge = &(edgeAux->nodeNext->edgePrev);
				while(*edge != edgeAux) {	
					edge = &((*edge)->edgeNext);
				}
				*edge = (*edge)->edgeNext;
				// create a new destination node and copy edges
				LNode *node = newNode(edgeAux->nodeNext->iFrame);
				++iNodesAdded;
				for(LEdge *edge = edgeAux->nodeNext->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
					if (edge == edgeAux) {
						continue;
					}
					LEdge *edgeDup = newEdge(edge);
					connectEdge(node,edgeDup,edge->nodeNext);	
					edgeDup->lexUnit = edge->lexUnit;
					edgeDup->iContextLeft = copyContext(iContextLeftProp);
					edgeDup->iContextRight = NULL;
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
		}
		delete [] iContextLeftProp;
	}
	if (iContextSizeMax > 1) {	
		delete [] edgeFake->iContextLeft;
	}
	delete edgeFake;	
	
	// (2) right context propagation
	
 	// create a fake edge to make the algorithm simpler
	edgeFake = new LEdge;
	edgeFake->nodePrev = m_lnodeFinal;
	if (iContextSizeMax > 1) {
		edgeFake->iContextRight = new unsigned char[iContextSizeMax];
		for(int i=0 ; i < iContextSizeMax ; ++i) {
			edgeFake->iContextRight[i] = m_phoneSet->getPhoneIndexSilence();
		}			
	}
	edgeFake->lexUnit = m_lexiconManager->getLexUnitSilence();
	lEdge.push_back(edgeFake);
	
	// process the edges one by one
	while(lEdge.empty() == false) {
	
		LEdge *edgeAux = lEdge.front();
		lEdge.pop_front();	
		
		// does the node come from the inital node?
		if (edgeAux->nodePrev->edgePrev == NULL) {
			assert(edgeAux->nodePrev == m_lnodeInitial);
			continue;	
		}
		
		// build the right context that will be propagated from the edge (cross-word)
		edgeAux->iPhones = (int)edgeAux->lexUnit->vPhones.size();
		unsigned char *iContextRightProp = newBlankContext();
		for(int i=0 ; i < min((int)m_iContextSizeCW,(int)(edgeAux->iPhones)) ; ++i) {
			iContextRightProp[i] = edgeAux->lexUnit->vPhones[i];
		}
		
		// predecessor edges do not have a context yet
		bool bContextRight = false;
		for(LEdge *edge = edgeAux->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
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
			for(LEdge *edge = edgeAux->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
				edge->iContextRight = copyContext(iContextRightProp);
				lEdge.push_back(edge);
			}
		} 
		// successor edges already have a context
		else {
			// different context: duplicate node and edges and put the edges in the queue
			if (memcmp(iContextRightProp,edgeAux->nodePrev->edgePrev->iContextRight,
				iContextSizeMax*sizeof(unsigned char)) != 0) {	
				// extract the auxiliar edge from the original source node
				LEdge **edge = &(edgeAux->nodePrev->edgeNext);
				while(*edge != edgeAux) {	
					edge = &((*edge)->edgePrev);
				}
				*edge = (*edge)->edgePrev;
				// create a new source node and copy edges
				LNode *node = newNode(edgeAux->nodePrev->iFrame);
				++iNodesAdded;
				for(LEdge *edge = edgeAux->nodePrev->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
					if (edge == edgeAux) {
						continue;
					}
					LEdge *edgeDup = newEdge(edge);
					connectEdge(edge->nodePrev,edgeDup,node);	
					edgeDup->lexUnit = edge->lexUnit;
					edgeDup->iContextLeft = copyContext(edge->iContextLeft);
					edgeDup->iContextRight = copyContext(iContextRightProp);
					assert((edgeDup->edgeNext != NULL) || (edgeDup->edgePrev != NULL));
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
		}
		delete [] iContextRightProp;
	}
	if (iContextSizeMax > 1) {	
		delete [] edgeFake->iContextRight;
	}
	delete edgeFake;	
	
	// rebuild the lattice container
	buildContainer();	
	
	// get the hmm-states for each edge
	for(int i=0 ; i < m_iEdges ; ++i) {
		
		LEdge *edge = m_ledges[i];
		
		//m_lexiconManager->print(edge->lexUnit);
		
		edge->iPhones = (int)edge->lexUnit->vPhones.size();
		if (edge->iPhones == 0) {
			m_lexiconManager->print(edge->lexUnit);
		}
		assert(edge->iPhones > 0);
		edge->phoneAlignment = new LPhoneAlignment[edge->iPhones];
		unsigned char *iContextLeft = new unsigned char[iContextSizeMax];
		unsigned char *iContextRight = new unsigned char[iContextSizeMax];
		unsigned char iPosition;
		for(int i=0 ; i < edge->iPhones ; ++i) {
		
			unsigned char iPhone = edge->lexUnit->vPhones[i];
			
			// build left context
			if (i == 0) {
				memcpy(iContextLeft,edge->iContextLeft,iContextSizeMax*sizeof(unsigned char));
			} else {
				for(int j=0 ; j < min(i,(int)m_iContextSizeWW); ++j) {
					iContextLeft[m_iContextSizeWW-j-1] = edge->lexUnit->vPhones[i-1-j];
				}
				for(int j=min(i,(int)m_iContextSizeWW) ; j < iContextSizeMax ; ++j) {
					iContextLeft[m_iContextSizeWW-j-1] = m_iPhoneContextPadding;
				}
			}
			//printContext(iContextLeft);
			
			// build right context
			if (i == edge->iPhones-1) {
				memcpy(iContextRight,edge->iContextRight,iContextSizeMax*sizeof(unsigned char));
			} else {
				for(int j=0 ; j < min(edge->iPhones-i-1,(int)m_iContextSizeWW); ++j) {
					iContextRight[j] = edge->lexUnit->vPhones[i+j+1];
				}
				for(int j=min(edge->iPhones-i-1,(int)m_iContextSizeWW) ; j < iContextSizeMax ; ++j) {
					iContextRight[j] = m_iPhoneContextPadding;
				}
			}
			//printContext(iContextRight);	
			
			// get the within-word position
			if (i == 0) {
				iPosition = (edge->iPhones == 1) ? WITHIN_WORD_POSITION_MONOPHONE : WITHIN_WORD_POSITION_START;
			} else if (i == edge->iPhones-1) {
				iPosition = WITHIN_WORD_POSITION_END;
			} else {
				iPosition = WITHIN_WORD_POSITION_INTERNAL;
			}
			// keep phone
			edge->phoneAlignment[i].iPhone = iPhone;
			edge->phoneAlignment[i].iPosition = iPosition;
			// get hmm-states	
			for(int iState=0 ; iState < NUMBER_HMM_STATES ; ++iState) {
				HMMStateDecoding *hmmState = hmmManager->getHMMStateDecoding(iContextLeft,iPhone,
					iContextRight,iPosition,iState);
				edge->phoneAlignment[i].iHMMState[iState] = hmmState->getId();
				//printf("hmm-id: %d\n",hmmState->getId());
				// within-phone first and last frames are still unknown
				edge->phoneAlignment[i].iStateBegin[iState] = -1;
				edge->phoneAlignment[i].iStateEnd[iState] = -1;	
			}	
		}
		delete [] iContextLeft;
		delete [] iContextRight;
	}
	
	// remove contexts
	for(int i=0 ; i < m_iEdges ; ++i) {
		assert(m_ledges[i]->iContextLeft != NULL);
		assert(m_ledges[i]->iContextRight != NULL);
	}	
	
	float fExpansionRatioNodes = ((float)(100*iNodesAdded))/((float)iNodesOriginal);
	float fExpansionRatioEdges = ((float)(100*iEdgesAdded))/((float)iEdgesOriginal);
	BVC_VERB << "# nodes added: " << iNodesAdded <<" (" << fExpansionRatioNodes << "%)";
	BVC_VERB << "# edges added: " << iEdgesAdded <<" (" << fExpansionRatioEdges << "%)";	
	
	// update the lattice properties
	setProperty(LATTICE_PROPERTY_HMMS,"yes");	
}

// WER computation ---------------------------------------------------------------------------

// compute the Lattice Word Error Rate (also known as oracle)
// - perfect match: prune all paths with edit errors (the goal is to get WER = 0 or failure)
// note: the alignment is not time synchronous so pruning based on best-partial-score is 
//       not possible (unless we are looking for a perfect-match)
LatticeWER *HypothesisLattice::computeWER(VLexUnit &vLexUnitReference, BestPath **bestPath, 
	Mappings *mappings, bool bPerfectMatch, int iBeamSize) {

	// make sure reference is not empty
	if (vLexUnitReference.empty()) {
		return NULL;
	}
	
	if (bPerfectMatch) {
		iBeamSize = 0;
	}
	
	// make sure none of the lexical units in the reference are fillers
	for(VLexUnit::iterator it = vLexUnitReference.begin() ; it != vLexUnitReference.end() ; ++it) {
		if (m_lexiconManager->isFiller(*it)) {
			return NULL;
		}
	}
	
	// count the number of OOV words
	int iOOV = 0;
	for(VLexUnit::iterator it = vLexUnitReference.begin() ; it != vLexUnitReference.end() ; ++it) {
		if (m_lexiconManager->isUnknown(*it)) {
			++iOOV;
		}
	}	
	
 	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	// create an array of active tokens
	m_tokenCurrent = new LWERToken[LATTICE_WER_MAX_TOKENS];
	m_tokenNext = new LWERToken[LATTICE_WER_MAX_TOKENS];
	m_iActiveTokensCurrent = 0;	
	
	// create the table of nodes used to keep the active tokens at each node
	m_iTokenNode = new int[m_iNodes];
	for(int i=0 ; i < m_iNodes ; ++i) {
		m_iTokenNode[i] = -1;
	}	
	
	// allocate memory for the history items and create the linked list of available items
	m_iHistoryItems = LATTICE_WER_HISTORY_ITEMS_START;
	m_historyItems = new LWERHistoryItem[m_iHistoryItems];
	for(int i=0 ; i < m_iHistoryItems-1 ; ++i) {
		m_historyItems[i].iPrev = i+1;
	}
	m_historyItems[m_iHistoryItems-1].iPrev = -1;
	m_iHistoryItemAvailable = 0;	
	
	// propagate tokens from the starting node
	m_tokenCurrent[0].iScore = 0;
	m_tokenCurrent[0].iNode = m_lnodeInitial->iNode;
	m_tokenCurrent[0].iReference = -1;
	m_tokenCurrent[0].iHistoryItem = -1;
	m_iActiveTokensCurrent=1;
	
	// propagate tokens until no token is active
	m_iActiveTokensNext = 0;	
	m_tokenBest.iNode = -1;
	m_tokenBest.iScore = INT_MAX;
	while(m_iActiveTokensCurrent > 0) {
	
		// propagate active tokens at the current time-frame
		m_iActiveTokensNext = 0;
		int iTokensToProcess = m_iActiveTokensCurrent;
		for(int i=0 ; ((i < LATTICE_WER_MAX_TOKENS) && (iTokensToProcess > 0)) ; ++i) {
	
			// skip inactive tokens (pruned)
			if (m_tokenCurrent[i].iNode == -1) {
				continue;
			}	
			
			--iTokensToProcess;
			
			// is the token at the final node and all the words in the reference seen?
			if ((m_tokenCurrent[i].iNode == m_lnodeFinal->iNode) && 
				(m_tokenCurrent[i].iReference == (int)vLexUnitReference.size()-1)) {
				// best scoring token so far?
				if (m_tokenCurrent[i].iScore < m_tokenBest.iScore) {
					memcpy(&m_tokenBest,&m_tokenCurrent[i],sizeof(LWERToken));
					//printf("best final token: %8.2f %4d\n",m_tokenBest.fScore,m_tokenBest.iReference);
				} 
				// first token with perfect score? if so, finish the search
				if (m_tokenBest.iScore == 0) {
					break;
				}
			}
			// propagate tokens through the successor nodes
			else {
				for(LEdge *edge = m_lnodes[m_tokenCurrent[i].iNode]->edgeNext ; edge != NULL ; edge = edge->edgePrev) {	
					if (m_tokenCurrent[i].iNode != m_lnodeFinal->iNode) {
						// skip hypothesis
						if (m_lexiconManager->isFiller(edge->lexUnit)) {
							activateToken(&m_tokenCurrent[i],edge,ALIGNMENT_EVENT_SKIP_HYP);
						} 
						// insertion
						else if (!bPerfectMatch) {
							activateToken(&m_tokenCurrent[i],edge,ALIGNMENT_EVENT_INSERTION);	
						}
					}
					if (m_tokenCurrent[i].iReference+1 < (int)vLexUnitReference.size()) {
						// skip reference
						if (m_lexiconManager->isFiller(vLexUnitReference[m_tokenCurrent[i].iReference+1])) {
							activateToken(&m_tokenCurrent[i],NULL,ALIGNMENT_EVENT_SKIP_REF);
						}
						// deletion
						else if (!bPerfectMatch) {
							activateToken(&m_tokenCurrent[i],NULL,ALIGNMENT_EVENT_DELETION);
						}
					}
					// substitution/correct
					if ((m_tokenCurrent[i].iNode != m_lnodeFinal->iNode) && 
						(m_tokenCurrent[i].iReference+1 < (int)vLexUnitReference.size()) &&
						(m_lexiconManager->isFiller(edge->lexUnit) == false) && 
						(m_lexiconManager->isFiller(vLexUnitReference[m_tokenCurrent[i].iReference+1]) == false)) {
						// no mappings
						if (mappings == NULL) {	
							if (edge->lexUnit->iLexUnit == vLexUnitReference[m_tokenCurrent[i].iReference+1]->iLexUnit) {
								activateToken(&m_tokenCurrent[i],edge,ALIGNMENT_EVENT_CORRECT);
							} else if (!bPerfectMatch) {
								activateToken(&m_tokenCurrent[i],edge,ALIGNMENT_EVENT_SUBSTITUTION);
							}
						}
						// apply mappings to hypothesis and reference word
						else {
							const char *strHyp = (*mappings)[m_lexiconManager->getStrLexUnit(edge->lexUnit->iLexUnit)];
							const char *strRef = (*mappings)[m_lexiconManager->getStrLexUnit(
								vLexUnitReference[m_tokenCurrent[i].iReference+1]->iLexUnit)];
							if (strcmp(strHyp,strRef) == 0) {
								activateToken(&m_tokenCurrent[i],edge,ALIGNMENT_EVENT_CORRECT);	
							} else if (!bPerfectMatch) {
								activateToken(&m_tokenCurrent[i],edge,ALIGNMENT_EVENT_SUBSTITUTION);
							}
						}
					}
				}	
			}
		}
		if (m_tokenBest.iScore == 0) {
			break;
		}
		
		// pruning ------------------------------------
		
		// (1) get the best score
		int iScoreBest = INT_MAX;
		if (m_tokenBest.iNode != -1) {
			iScoreBest = m_tokenBest.iScore;
		}	
		for(int i=0 ; i < m_iActiveTokensNext ; ++i) {
			if (m_tokenNext[i].iScore < iScoreBest) {
				iScoreBest = m_tokenNext[i].iScore;
			}
			//printf("active: %6d (score: %6.2f) reference: %6d node: %x event: %d\n",m_tokenNext[i].node->iNode,
			//	m_tokenNext[i].fScore,m_tokenNext[i].iReference,m_tokenNext[i].node,m_tokenNext[i].historyItem->iAlignmentEvent);
			m_iTokenNode[m_tokenNext[i].iNode] = -1;	
		}
		// build histogram
		int *iHistogram = new int[iBeamSize+1];
		for(int i=0 ; i < iBeamSize+1 ; ++i) {
			iHistogram[i] = 0;
		}	
		for(int i=0 ; i < m_iActiveTokensNext ; ++i) {
			int iBin = min(m_tokenNext[i].iScore-iScoreBest,iBeamSize);
			assert((iBin >= 0) && (iBin <= iBeamSize));
			iHistogram[iBin]++;
		}	
		int iThreshold = 0;
		int iSeen = 0;
		for(int i=0 ; i < iBeamSize+1 ; ++i) {
			iSeen += iHistogram[i];
			if (iSeen > LATTICE_WER_MAX_ACTIVE_TOKENS) {
				break;
			}
			iThreshold += 1;
		}
		delete [] iHistogram;	
		// (2) prune active tokens outside the beam
		int iPruned = 0;
		for(int i=0 ; i < m_iActiveTokensNext ; ++i) {	
			if (m_tokenNext[i].iScore-iScoreBest > iThreshold) {
				m_tokenNext[i].iNode = -1;
				++iPruned;
			}
		}
		m_iActiveTokensNext -= iPruned;	
		
		// sanity check
		for(int i=0 ; i < m_iNodes ; ++i) {
			assert(m_iTokenNode[i] == -1);
		}
		
		//printf("bestScore: %5.2f (reference: %4d) active states: %6d pruned: %6d\n",
		//	fScoreBest,iReferenceBest,m_iActiveTokensNext,iPruned);
		
		// swap
		m_iActiveTokensCurrent = m_iActiveTokensNext;
		LWERToken *tokenAux = m_tokenCurrent;
		m_tokenCurrent = m_tokenNext;
		m_tokenNext = tokenAux;
	}
	
	// unable to find a perfect match? (WER > 0)
	if (m_tokenBest.iNode == -1) {
		assert(bPerfectMatch);
		// clean-up
		delete [] m_tokenCurrent;
		delete [] m_tokenNext;
		delete [] m_iTokenNode;	
		delete [] m_historyItems;	
		return NULL;
	}	
	
	// if the best path is already marked clear it
	if (isProperty(LATTICE_PROPERTY_BEST_PATH)) {
		clearBestPath();
	}
	
	// do back-tracing on the best token to compute the Lattice WER
	LatticeWER *latticeWER = new LatticeWER;
	latticeWER->iWordsReference = (int)vLexUnitReference.size();
	latticeWER->iInsertions = 0;
	latticeWER->iDeletions = 0;
	latticeWER->iSubstitutions = 0;
	latticeWER->iCorrect = 0;
	latticeWER->iErrors = 0;	
	
	if (bestPath) {
		*bestPath = new BestPath(m_lexiconManager,0.0);
	}
	
	assert(m_tokenBest.iNode != -1);
	int iHistoryItem = m_tokenBest.iHistoryItem;
	while(iHistoryItem != -1) {
	
		LWERHistoryItem *historyItem = m_historyItems+iHistoryItem;	
	
		// keep and mark the best path in the lattice
		if (historyItem->iAlignmentEvent != ALIGNMENT_EVENT_DELETION) {
			latticeWER->vLexUnitBest.push_back(m_ledges[historyItem->iEdge]->lexUnit);	
			m_ledges[historyItem->iEdge]->bBestPath = true;
		}	
			
		// correct
		if (historyItem->iAlignmentEvent == ALIGNMENT_EVENT_CORRECT) {
			//cout << "correct:	" << m_lexiconManager->getStrLexUnit(historyItem->edge->lexUnit->iLexUnit) <<
			//	" " << m_lexiconManager->getStrLexUnit(vLexUnitReference[historyItem->iReference]->iLexUnit) << endl;			
			++latticeWER->iCorrect;
		} 
		// insertion	
		else if (historyItem->iAlignmentEvent == ALIGNMENT_EVENT_INSERTION) {
			//cout << "insertion: " << m_lexiconManager->getStrLexUnit(historyItem->edge->lexUnit->iLexUnit) << endl;
			++latticeWER->iInsertions;
		}
		// deletion
		else if (historyItem->iAlignmentEvent == ALIGNMENT_EVENT_DELETION) {
			//cout << "deletion: " << m_lexiconManager->getStrLexUnit(vLexUnitReference[historyItem->iReference]->iLexUnit) 
			//	<< endl;
			++latticeWER->iDeletions;
		}
		// substitution
		else if (historyItem->iAlignmentEvent == ALIGNMENT_EVENT_SUBSTITUTION) {
			//cout << "substitution: " << m_lexiconManager->getStrLexUnit(historyItem->edge->lexUnit->iLexUnit) <<
			//	m_lexiconManager->getStrLexUnit(vLexUnitReference[historyItem->iReference]->iLexUnit) << endl;
			++latticeWER->iSubstitutions;
		}
		// skip hypothesis
		else if (historyItem->iAlignmentEvent == ALIGNMENT_EVENT_SKIP_HYP) {
			//cout << "skip: " << m_lexiconManager->getStrLexUnit(historyItem->edge->lexUnit->iLexUnit) << endl;
		}
		// skip reference
		else {
			assert(historyItem->iAlignmentEvent == ALIGNMENT_EVENT_SKIP_REF);
			//cout << "skip: " << m_lexiconManager->getStrLexUnit(vLexUnitReference[historyItem->iReference]->iLexUnit) << endl;
		}
				
		// add it to the best path
		if ((bestPath) && (historyItem->iAlignmentEvent != ALIGNMENT_EVENT_DELETION) &&
			(historyItem->iAlignmentEvent != ALIGNMENT_EVENT_SKIP_REF)) {
			float fScore = 0.0;
			float fScoreAM = isProperty(LATTICE_PROPERTY_AM_PROB) ? m_ledges[historyItem->iEdge]->fScoreAM : 0.0;
			float fScoreLM = isProperty(LATTICE_PROPERTY_LM_PROB) ? m_ledges[historyItem->iEdge]->fScoreLM : 0.0;
			float fConfidence = isProperty(LATTICE_PROPERTY_CONFIDENCE) ? m_ledges[historyItem->iEdge]->fConfidence : 0.0;
			float fInsertionPenalty = 
				isProperty(LATTICE_PROPERTY_INSERTION_PENALTY) ? m_ledges[historyItem->iEdge]->fInsertionPenalty : 0.0;
			(*bestPath)->newElementFront(m_ledges[historyItem->iEdge]->iFrameStart,m_ledges[historyItem->iEdge]->iFrameEnd,
				fScore,fScoreAM,fScoreLM,fConfidence,m_ledges[historyItem->iEdge]->lexUnit,fInsertionPenalty);
		}
				
		iHistoryItem = (m_historyItems+iHistoryItem)->iPrev;
	}
	
	// compute WER from alignment errors
	latticeWER->iErrors = latticeWER->iInsertions+latticeWER->iDeletions+latticeWER->iSubstitutions;
	latticeWER->fWER = ((float)latticeWER->iErrors*100.0f)/((float)latticeWER->iWordsReference);	
	latticeWER->iOOV = iOOV;
	
	// set the best-path property
	setProperty(LATTICE_PROPERTY_BEST_PATH,"yes");
		
	// clean-up
	delete [] m_tokenCurrent;
	delete [] m_tokenNext;
	delete [] m_iTokenNode;	
	delete [] m_historyItems;
	
 	double dTimeEnd = TimeUtils::getTimeMilliseconds();
 	double dTime = (dTimeEnd-dTimeBegin)/1000.0;
 	double dRTF = dTime/(m_iFrames/100.0);
 	
	printf("lattice WER: %5.2f computation time: %.4fs (RTF= %5.4f) sub: %4d del: %4d ins: %4d correct: %4d\n",
		latticeWER->fWER,dTime,dRTF,latticeWER->iSubstitutions,latticeWER->iDeletions,latticeWER->iInsertions,
		latticeWER->iCorrect);
	
	return latticeWER;
}

// history item garbage collection: free-up unused history items
void HypothesisLattice::historyItemGarbageCollection() {

	// (1) initialize all the history items as unused
	for(int i=0 ; i < m_iHistoryItems ; ++i) {
		m_historyItems[i].bUsed = false;
	}

	// (2) mark history items in use as used
	int iItemsUsed = 0;
	
	// current tokens
	int iTokensToProcess = m_iActiveTokensCurrent;
	for(int i=0 ; ((i < LATTICE_WER_MAX_TOKENS) && (iTokensToProcess > 0)) ; ++i) {

		// skip inactive tokens (pruned)
		if (m_tokenCurrent[i].iNode == -1) {
			continue;
		}	
		
		--iTokensToProcess;
		
		for(int iHistoryItem = m_tokenCurrent[i].iHistoryItem ; iHistoryItem != -1 ; iHistoryItem = (m_historyItems+iHistoryItem)->iPrev) {
			LWERHistoryItem *historyItem = m_historyItems+iHistoryItem;
			if (historyItem->bUsed) {
				break;
			}
			historyItem->bUsed = true;
			++iItemsUsed;
		}
	}

	// next tokens
	for(int i=0 ; i < m_iActiveTokensNext ; ++i) {
		
		for(int iHistoryItem = m_tokenNext[i].iHistoryItem ; iHistoryItem != -1 ; iHistoryItem = (m_historyItems+iHistoryItem)->iPrev) {
			LWERHistoryItem *historyItem = m_historyItems+iHistoryItem;
			if (historyItem->bUsed) {
				break;
			}
			historyItem->bUsed = true;
			++iItemsUsed;
		}
	}
	
	// best token
	if (m_tokenBest.iNode != -1) {
		for(int iHistoryItem = m_tokenBest.iHistoryItem ; iHistoryItem != -1 ; iHistoryItem = (m_historyItems+iHistoryItem)->iPrev) {
			LWERHistoryItem *historyItem = m_historyItems+iHistoryItem;
			if (historyItem->bUsed) {
				break;
			}
			historyItem->bUsed = true;
			++iItemsUsed;
		}	
	}

	assert(iItemsUsed <= m_iHistoryItems);
	if (iItemsUsed > (0.2*m_iHistoryItems)) {

		//printf("garbageCollector: reallocating to %d elements\n",m_iHistoryItems*2);
		
		// allocate a new data structure with double capacity		
		LWERHistoryItem *historyItems = new LWERHistoryItem[m_iHistoryItems*2];
		
		// copy the active items from the old data structure
		for(int i=0 ; i < m_iHistoryItems ; ++i) {
			historyItems[i].iPrev = m_historyItems[i].iPrev;
			historyItems[i].iEdge = m_historyItems[i].iEdge;
			historyItems[i].iReference = m_historyItems[i].iReference;
			historyItems[i].iAlignmentEvent = m_historyItems[i].iAlignmentEvent;
			historyItems[i].bUsed = m_historyItems[i].bUsed;
		}
		
		// create the linked list of available items
		for(int i=m_iHistoryItems ; i < (2*m_iHistoryItems)-1 ; ++i) {
			historyItems[i].iPrev = i+1;
			historyItems[i].bUsed = false;	
		}
		historyItems[(2*m_iHistoryItems)-1].iPrev = -1;
		historyItems[(2*m_iHistoryItems)-1].bUsed = false;
		
		delete [] m_historyItems;
		m_historyItems = historyItems;
		m_iHistoryItemAvailable = m_iHistoryItems;
		m_iHistoryItems *= 2;		
		
	} 
	// just recycle unused items
	else {
		//printf("garbageCollector: recycling (%d items used)\n",iItemsUsed);
		
		// put unused items in the list of available items
		int *iHistoryItemAux = &m_iHistoryItemAvailable;
		for(int i=0 ; i < m_iHistoryItems ; ++i) {
			if (m_historyItems[i].bUsed == false) {
				*iHistoryItemAux = i;
				iHistoryItemAux = &(m_historyItems[i].iPrev);
			}
		}
		*iHistoryItemAux = -1;
		assert(m_iHistoryItemAvailable != -1);
	}	
}

// rescore the lattice using the given rescoring method (Dijkstra)
// - all edges in the lattice contain the cost of moving from the source node to the 
//   destination node, that cost is expressed in terms of log-likelihood (like in 
//   standard decoding) or a posterior probability
// - the algorithm works by computing the minimum distance from each node to the final
//   node using Dijkstra 
BestPath *HypothesisLattice::rescore(const char *strRescoringMethod) {

	// determine the function to use for getting the edge log-weight
	bool bLikelihood = false;
	if (strcmp(strRescoringMethod,RESCORING_METHOD_LIKELIHOOD) == 0) {
		bLikelihood = true;
		// check lattice properties
		if (!isProperty(LATTICE_PROPERTY_AM_PROB) || 
			!isProperty(LATTICE_PROPERTY_LM_PROB) ||
			!isProperty(LATTICE_PROPERTY_INSERTION_PENALTY)) {
			BVC_ERROR << "likelihood-based rescoring needs acoustic and lm likelihood and insertion penalty";
		}
	} else {
		assert(strcmp(strRescoringMethod,RESCORING_METHOD_POSTERIORS) == 0);
		bLikelihood = false;
		// check lattice properties
		if (!isProperty(LATTICE_PROPERTY_PP)) {
			BVC_ERROR << "pp-based rescoring needs posterior probabilities";
		}
	} 
	
	// allocate a structure to keep the node states	
	// initialize distances and traceback pointers
	unsigned char *iNodeState = new unsigned char[m_iNodes];
	double *dNodeDistance = new double[m_iNodes];
	LEdge **edgeNext = new LEdge*[m_iNodes];
	for(int i=0 ; i < m_iNodes ; ++i) {
		iNodeState[i] = LATTICE_NODE_STATE_UNSEEN;
		dNodeDistance[i] = -DBL_MAX;
		edgeNext[i] = NULL;
	}
	
	// (1) compute minimum distance from each node to the final node
	
	// start from the final node
	LLNode lNodes;
	lNodes.push_back(m_lnodeFinal);
	iNodeState[m_lnodeFinal->iNode] = LATTICE_NODE_STATE_QUEUED;
	edgeNext[m_lnodeFinal->iNode] = NULL;
	dNodeDistance[m_lnodeFinal->iNode] = 0.0;
	
	// process all nodes
	while(!lNodes.empty()) {
	
		// get a node to process
		LNode *node = lNodes.front();
		lNodes.pop_front();
		
		iNodeState[node->iNode] = LATTICE_NODE_STATE_PROCESSED;

		// compute the minimum distance from each node to the final node
		for(LEdge *edge = node->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		
			// relax?
			double dAcc;
			if (bLikelihood) {
				dAcc = edgeLogLikelihood(edge);
			} else {
				dAcc = edgeLogPP(edge);
			}
			dAcc += dNodeDistance[node->iNode];
			if (dNodeDistance[edge->nodePrev->iNode] < dAcc) {
				dNodeDistance[edge->nodePrev->iNode] = dAcc;
				edgeNext[edge->nodePrev->iNode] = edge;
				if (iNodeState[edge->nodePrev->iNode] != LATTICE_NODE_STATE_QUEUED) {
					lNodes.push_back(edge->nodePrev);
					iNodeState[edge->nodePrev->iNode] = LATTICE_NODE_STATE_QUEUED;
				}
			}
		}
	}
	
	assert(edgeNext[m_lnodeInitial->iNode] != NULL);
	
	// mark all edges as not in the best path	
	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->bBestPath = false;
	}	
	
	// (2) retrieve best path using traceback pointers
	BestPath *bestPath = new BestPath(m_lexiconManager,0.0);
	LEdge *edgeBest = edgeNext[m_lnodeInitial->iNode];
	do {
		// mark best path within the lattice
		edgeBest->bBestPath = true;
		// keep best path
		bestPath->newElementBack(edgeBest->iFrameStart,edgeBest->iFrameEnd,0.0,
			edgeBest->fScoreAM,edgeBest->fScoreLM,edgeBest->fConfidence,
			edgeBest->lexUnit,edgeBest->fInsertionPenalty);
		edgeBest = edgeNext[edgeBest->nodeNext->iNode];	
	} while(edgeBest != NULL);
	
	//bestPath->print();
	
	delete [] dNodeDistance;
	delete [] edgeNext;
	delete [] iNodeState;
	
	setProperty(LATTICE_PROPERTY_BEST_PATH,"yes");

	return bestPath;	
}



// return the best path in the lattice (must be already marked)
BestPath *HypothesisLattice::getBestPath() {

	if (isProperty(LATTICE_PROPERTY_BEST_PATH) == false) {
		return NULL;
	}

	BestPath *bestPath = new BestPath(m_lexiconManager,0.0);
	LNode *node = m_lnodeInitial;
	do {
		LEdge *edgeBest = NULL;
		for(LEdge *edge = node->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (edge->bBestPath) {
				edgeBest = edge;
				break;
			}
		}
		assert(edgeBest != NULL);
		
		// keep best path
		bestPath->newElementBack(edgeBest->iFrameStart,edgeBest->iFrameEnd,0.0,
			edgeBest->fScoreAM,edgeBest->fScoreLM,edgeBest->fConfidence,
			edgeBest->lexUnit,edgeBest->fInsertionPenalty);	
		
		node = edgeBest->nodeNext;
	} while(node != m_lnodeFinal);

	return bestPath;
}

// return the alignment of the best path in the lattice (must be already marked)
VLPhoneAlignment *HypothesisLattice::getBestPathAlignment() {
	
	if ((!isProperty(LATTICE_PROPERTY_BEST_PATH)) || (!isProperty(LATTICE_PROPERTY_PHONE_ALIGN))) {
		BVC_ERROR << "best-path is not marked in the lattice and/or phone alignment is not present!";
	}

	VLPhoneAlignment *vLPhoneAlignment = new VLPhoneAlignment;
	
	LNode *node = m_lnodeInitial;
	do {
		LEdge *edgeBest = NULL;
		for(LEdge *edge = node->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
			if (edge->bBestPath) {
				edgeBest = edge;
				break;
			}
		}
		assert(edgeBest != NULL);	
		
		for(int i=0 ; i < edgeBest->iPhones ; ++i) {
			vLPhoneAlignment->push_back(&edgeBest->phoneAlignment[i]);
		}
		
		node = edgeBest->nodeNext;
	} while(node != m_lnodeFinal);	

	return vLPhoneAlignment;
}


// attach lm-probabilities to the edges of the lattice
// note: lattice expansion may be needed in order for each edge to have a unique
//       word context
void HypothesisLattice::attachLMProbabilities(LMFSM *lmFSM) {

	LLEdge lEdge;
	int iNodesAdded = 0;
	int iEdgesAdded = 0;
	int iEdgesProcessed = 0;
	
 	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	// reset lm-state and score for all edges
	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->iLMState = -1;
		m_ledges[i]->iLMStatePrev = -1;
		m_ledges[i]->fScoreLM = 0.0;
	}
	
	// get the initial lm-state	
	int iLMStateInitial = lmFSM->getInitialState();
	
	// put initial edges in the queue
	for(LEdge *edge = m_lnodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		// regular word: update lm-state 
		if (m_lexiconManager->isStandard(edge->lexUnit)) {		
			edge->iLMState = lmFSM->updateLMState(iLMStateInitial,edge->lexUnit->iLexUnit,&edge->fScoreLM);
		} else {
			edge->iLMState = iLMStateInitial;
			edge->fScoreLM = 0.0;
		}
		lEdge.push_back(edge);
	}
	
	// process edges one by one
	while(lEdge.empty() == false) {
	
		LEdge *edge = lEdge.front();
		lEdge.pop_front();
		++iEdgesProcessed;
		
		// edge goes to final node: add the lm-score resulting from transitioning to the final state
		if (edge->nodeNext == m_lnodeFinal) {
			edge->fScoreLM += lmFSM->toFinalState(edge->iLMState);
			continue;
		}
		
		//either all the next edges have a lm-state or none of them have, if they have then check that all the destination //edges have the right word-context, otherwise make a copy of all of them and create a new destination node too
		
		assert(edge->nodeNext->edgeNext != NULL);
		// successor edges do not have word-context:	attach lm-state and score to them
		if (edge->nodeNext->edgeNext->iLMStatePrev == -1) {
			// attach lm-state and score to successor edges
			for(LEdge *edge2 = edge->nodeNext->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
				// regular word: update lm-state 
				if (m_lexiconManager->isStandard(edge2->lexUnit)) {
					edge2->iLMState = lmFSM->updateLMState(edge->iLMState,edge2->lexUnit->iLexUnit,&edge2->fScoreLM);	
				} 
				// silence/filler: lm state does not change
				else {
					edge2->iLMState = edge->iLMState;
					edge2->fScoreLM = 0.0;
				}
				edge2->iLMStatePrev = edge->iLMState;
				lEdge.push_back(edge2);
			}	
		}
		// successor edges have a word context
		// (a) same previous lm-state: do nothing
		else if (edge->nodeNext->edgeNext->iLMStatePrev == edge->iLMState) {
			// sanity check
			for(LEdge *edge2 = edge->nodeNext->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
				assert(edge2->iLMStatePrev == edge->iLMState);
			}
		} 
		// (b) different previous lm-state: create a new destination node and copy successors
		else {
					
			// extract the auxiliar edge from the original destination node
			LEdge **edgeAux = &(edge->nodeNext->edgePrev);
			while(*edgeAux != edge) {	
				edgeAux = &((*edgeAux)->edgeNext);
			}
			*edgeAux = (*edgeAux)->edgeNext;	
			// copy successors
			LNode *node = newNode(edge->nodeNext->iFrame);
			++iNodesAdded;
			for(LEdge *edge2 = edge->nodeNext->edgeNext ; edge2 != NULL ; edge2 = edge2->edgePrev) {
				assert(edge2 != edge);
				LEdge *edgeDup = newEdge(edge2);
				connectEdge(node,edgeDup,edge2->nodeNext);	
				edgeDup->lexUnit = edge2->lexUnit;
				// attach lm-state and lm-score
				// - regular word: update lm-state 
				if (m_lexiconManager->isStandard(edgeDup->lexUnit)) {
					edgeDup->iLMState = lmFSM->updateLMState(edge->iLMState,edgeDup->lexUnit->iLexUnit,&edgeDup->fScoreLM);	
				} 
				// - silence/filler: lm state does not change
				else {
					edgeDup->iLMState = edge->iLMState;
					edgeDup->fScoreLM = 0.0;
				}	
				assert((edgeDup->edgePrev != NULL) || (edgeDup->edgeNext != NULL));
				lEdge.push_back(edgeDup);
				++iEdgesAdded;
			}
			edge->nodeNext = node;
			node->edgePrev = edge;
			edge->edgeNext = NULL;	
		}
	}
	
	BVC_VERB << "edges processed: " << iEdgesProcessed;
	float fExpansionRatioNodes = ((float)(100*iNodesAdded))/((float)m_iNodes);
	float fExpansionRatioEdges = ((float)(100*iEdgesAdded))/((float)m_iEdges);
	BVC_VERB << "expansion ratio:";
	BVC_VERB << "# nodes added: " << iNodesAdded <<" (" << fExpansionRatioNodes << "%)";
	BVC_VERB << "# edges added: " << iEdgesAdded <<" (" << fExpansionRatioEdges << "%)";
	
	// update lattice properties
	setProperty(LATTICE_PROPERTY_LM_PROB,"yes");
	setProperty(LATTICE_PROPERTY_NGRAM,LMManager::getStrNGram(lmFSM->getNGramOrder()));

	// rebuild the lattice container
	buildContainer(m_lnodeInitial,m_lnodeFinal);
	
	check();
	
	// compute processing time
 	double dTimeEnd = TimeUtils::getTimeMilliseconds();
 	double dTime = (dTimeEnd-dTimeBegin)/1000.0;
 	double dRTF = dTime/(m_iFrames/100.0);	
 	BVC_VERB << "processing time: " << dTime << " seconds, RTF: " << dRTF;
}

// attach insertion penalties
void HypothesisLattice::attachInsertionPenalty(LexiconManager *lexiconManager) {

	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->fInsertionPenalty = m_ledges[i]->lexUnit->fInsertionPenalty;
	}
	
	// update properties
	setProperty(LATTICE_PROPERTY_INSERTION_PENALTY,"yes");
}

// clear the best path in the lattice
void HypothesisLattice::clearBestPath() {

	assert(isProperty(LATTICE_PROPERTY_BEST_PATH));
	
	// this is a safer method
	for(int i=0 ; i < m_iEdges ; ++i) {
		m_ledges[i]->bBestPath = false;
	}	
	
	removeProperty(LATTICE_PROPERTY_BEST_PATH);
}

// add a path to the lattice
void HypothesisLattice::addPath(Alignment *alignment, bool bIsBest) {

	// check consistency
	if (alignment->getFrames() != (unsigned int)m_iFrames) {
		BVC_ERROR << "inconsistent number of frames: lattice has " << m_iFrames << ", path to insert has " << alignment->getFrames();
	}

	// if the added path will be marked as best path then the original best path needs to be cleared
	if (bIsBest) {
		if (isProperty(LATTICE_PROPERTY_BEST_PATH)) {
			clearBestPath();
		}	
	}

	LNode *node = m_lnodeInitial;
	VWordAlignment *vWordAlignment = alignment->getWordAlignment();
	for(VWordAlignment::iterator it = vWordAlignment->begin() ; it != vWordAlignment->end() ; ++it) {
		LEdge *edge = newEdge((*it)->iFrameBegin,(*it)->iFrameEnd,
			m_lexiconManager->getLexUnitPron((*it)->iLexUnitPron),-FLT_MAX,-FLT_MAX,-FLT_MAX);
		edge->nodePrev = node;
		edge->edgePrev = node->edgeNext;
		node->edgeNext = edge;
		VWordAlignment::iterator jt = it;
		// link to final node
		if (++jt == vWordAlignment->end()) {
			edge->nodeNext = m_lnodeFinal;
			edge->edgeNext = m_lnodeFinal->edgePrev;
			m_lnodeFinal->edgePrev = edge;
		} 
		// link to auxiliar node
		else {
			node = newNode(edge->iFrameEnd);
			edge->nodeNext = node;
			node->edgePrev = edge;
			edge->edgeNext = NULL;
		}	
		if (bIsBest) {
			edge->bBestPath = true;
		} else {
			edge->bBestPath = false;
		}
	}
	
	// rebuild the lattice container
	buildContainer(m_lnodeInitial,m_lnodeFinal);
	
	check();
	
	// TODO: do forward/backward compacting in order to merge equivalent edges if possible
	
	// update properties
	if (bIsBest) {
		setProperty(LATTICE_PROPERTY_BEST_PATH,"yes");
	}
}

// compute lattice likelihood
double HypothesisLattice::getLikelihood() {

	if (!isProperty(LATTICE_PROPERTY_FWD_PROB)) {
		BVC_ERROR << "forward scores are needed in order to compute lattice likelihood!";
	}

	// add probability from terminal edges
	double dLikelihood = -DBL_MAX;
	for(LEdge *edge = m_lnodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		dLikelihood = Numeric::logAddition(dLikelihood,edge->dScoreForward);	
	}
	
	return dLikelihood;
}

// compute the phone accuracy for phones within each arc given the reference word sequence
void HypothesisLattice::computePhoneAccuracy(VLPhoneAlignment &vLPhoneAlignment, bool bSetSilenceToZero, 
	bool bSetFillersToZero) {

	if (isProperty(LATTICE_PROPERTY_PHONE_ACCURACY)) {
		for(int i=0 ; i < m_iEdges ; ++i) {
			assert(m_ledges[i]->fPhoneAccuracy != NULL);
			delete [] m_ledges[i]->fPhoneAccuracy;
			m_ledges[i]->fPhoneAccuracy = NULL;
		}
	}

	// for each edge
	for(int i=0 ; i < m_iEdges ; ++i) {
		LEdge *edge = m_ledges[i];
		edge->fPhoneAccuracy = new float[edge->iPhones];
		for(int iPhone = 0 ; iPhone < edge->iPhones ; ++iPhone) {
			if ((bSetSilenceToZero && (edge->lexUnit == m_lexiconManager->getLexUnitSilence())) || 
				(bSetFillersToZero && (m_lexiconManager->isFiller(edge->lexUnit)))) {
				edge->fPhoneAccuracy[iPhone] = 0.0;
				continue;
			}
			LPhoneAlignment *phoneAlignment = &edge->phoneAlignment[iPhone];	
			// for each phoneme in the reference that overlap with this phoneme
			float fAccuracyMax = -FLT_MAX;
			for(VLPhoneAlignment::iterator it = vLPhoneAlignment.begin() ; it != vLPhoneAlignment.end() ; ++it) {
				float fOverlap = getOverlapProportion(*it,phoneAlignment);
				assert(fOverlap >= 0.0);
				assert(fOverlap <= 1.1);
				if (fOverlap > 0.0) {
					float fAccuracy;
					// same phone
					if (phoneAlignment->iPhone == (*it)->iPhone) {
						fAccuracy = -1.0f + 2.0f*fOverlap;
					}
					// different phone
					else {
						fAccuracy = -1.0f + fOverlap;	
					}
					if (fAccuracy > fAccuracyMax) {
						fAccuracyMax = fAccuracy;
					}
				}
			}
			assert(fAccuracyMax != -FLT_MAX);
			edge->fPhoneAccuracy[iPhone] = fAccuracyMax;
		}	
	}
	
	setProperty(LATTICE_PROPERTY_PHONE_ACCURACY,"yes");
}

// n-best list generation -----------------------------------------------

// create a n-best list from the lattice
NBestList *HypothesisLattice::createNBestList(int iN, const char *strRescoringMethod) {

	// determine the function to use for getting the edge log-weight
	bool bLikelihood = false;
	if (strcmp(strRescoringMethod,RESCORING_METHOD_LIKELIHOOD) == 0) {
		bLikelihood = true;
		// check lattice properties
		if (!isProperty(LATTICE_PROPERTY_AM_PROB) || 
			!isProperty(LATTICE_PROPERTY_LM_PROB) ||
			!isProperty(LATTICE_PROPERTY_INSERTION_PENALTY)) {
			BVC_ERROR << "likelihood-based rescoring needs acoustic and lm likelihood and insertion penalty";
		}
	} else {
		assert(strcmp(strRescoringMethod,RESCORING_METHOD_POSTERIORS) == 0);
		bLikelihood = false;
		// check lattice properties
		if (!isProperty(LATTICE_PROPERTY_PP)) {
			BVC_ERROR << "pp-based rescoring needs posterior probabilities";
		}
	} 
	
	resetAuxEdges();
	
	// compute Viterbi paths and scores
	untouchEdges();
	assert(m_lnodeFinal->edgePrev != NULL);
	for(LEdge *edge = m_lnodeFinal->edgePrev ; edge != NULL ; edge = edge->edgeNext) {
		viterbi(edge);
	}	

	// compute reverse Viterbi paths and scores
	untouchEdges();
	assert(m_lnodeInitial->edgeNext != NULL);
	for(LEdge *edge = m_lnodeInitial->edgeNext ; edge != NULL ; edge = edge->edgePrev) {
		viterbiReverse(edge);
	}	
	
	// create the n-best
	untouchEdges();
	// sort edges by path score
	LLEdge lEdge;
	for(int i=0 ; i<m_iEdges ; ++i) {
		assert((m_ledges[i]->edgePrevAux) || (m_ledges[i]->nodePrev == m_lnodeInitial));
		assert((m_ledges[i]->edgeNextAux) || (m_ledges[i]->nodeNext == m_lnodeFinal));
		lEdge.push_back(m_ledges[i]);
	}
	lEdge.sort(comparePathScore);
		
	// each entry in the n-best list is generated by finding the unmarked 
	// edge in the lattice with the highest score and then getting its path
	NBestList *nBestList = new NBestList();
	for(int i = 0 ; i < iN ; ++i) {
	
		// get the best unmarked edge in the lattice (it is not in any n-best entry yet)
		LEdge *edgeBest = NULL;
		for(LLEdge::iterator it = lEdge.begin() ; it != lEdge.end() ; ) {
			if ((*it)->lexUnit->iType != LEX_UNIT_TYPE_STANDARD) {
				++it;
				continue;
			}
			if ((*it)->bTouched == true) {
				it = lEdge.erase(it);
			} else {
				edgeBest = *it;
				break;
			}
		}
		if (!edgeBest) {
			break;
		}
			
		// recover path using the trace-back and trace-forward pointers
		LLEdge lPath;
		LEdge *edgeAux = edgeBest;
		while(edgeAux->nodePrev != m_lnodeInitial) {
			edgeAux->edgePrevAux->bTouched = true;
			lPath.push_front(edgeAux->edgePrevAux);
			edgeAux = edgeAux->edgePrevAux;	
		}
		lPath.push_back(edgeBest);
		edgeBest->bTouched = true;
		edgeAux = edgeBest;
		while(edgeAux->nodeNext != m_lnodeFinal) {
			edgeAux->edgeNextAux->bTouched = true;
			lPath.push_back(edgeAux->edgeNextAux);
			edgeAux = edgeAux->edgeNextAux;	
		}
		
		// add the path to the n-best list
		double dLikelihoodPath = edgeBest->dScoreViterbi + edgeBest->dScoreViterbiReverse;
		NBestListEntry *nBestListEntry = new NBestListEntry(m_lexiconManager);
		for(LLEdge::iterator jt = lPath.begin() ; jt != lPath.end() ; ++jt) {
			double dLikelihoodAM = isProperty(LATTICE_PROPERTY_AM_PROB) ? (*jt)->fScoreAM : 0.0;
			double dLikelihoodLM = isProperty(LATTICE_PROPERTY_LM_PROB) ? (*jt)->fScoreLM : 0.0;
			double dIP = (*jt)->lexUnit->fInsertionPenalty;
			double dPP = isProperty(LATTICE_PROPERTY_PP) ? (*jt)->fPP : 0.0;
			nBestListEntry->add(new NBestListEntryElement((*jt)->iFrameStart,(*jt)->iFrameEnd,
				(*jt)->lexUnit,dLikelihoodAM,dLikelihoodLM,dIP,dPP));	
		}
		nBestListEntry->setLikelihood(dLikelihoodPath);
		nBestList->add(nBestListEntry);
	}
	
	return nBestList;
}
		
// viterbi: for each edge keep predecessor edge and path score up to the edge
void HypothesisLattice::viterbi(LEdge *edge) {

	if (edge->nodePrev == m_lnodeInitial) {
		edge->bTouched = true;
		edge->dScoreViterbi = (edge->fScoreAM+edge->fInsertionPenalty)*m_fAMScalingFactor + edge->fScoreLM*m_fLMScalingFactor;
		edge->edgePrevAux = NULL;
		return;
   }
   
	// (1) get the Viterbi score from predecessors
	double dScoreMax = -DBL_MAX;
	double dAux = 0.0;	
	LEdge *edgePred = edge->nodePrev->edgePrev;
	while(edgePred) {
		
		// compute the Viterbi score of the predecessor if not available
		if (edgePred->bTouched == false) {
			viterbi(edgePred);	
		} 
		assert(edgePred->bTouched);
		assert(edgePred->dScoreViterbi != -DBL_MAX);	
		
		// keep best score and best predecessor edge
		dAux = edgePred->dScoreViterbi + edge->fScoreLM*m_fLMScalingFactor;
		if (dAux > dScoreMax) {
			dScoreMax = dAux;
			edge->edgePrevAux = edgePred;
		}
		
		edgePred = edgePred->edgeNext;
	}
	assert(edge->edgePrevAux);
 
	// (2) compute the viterbi score of the edge (best score)
	edge->dScoreViterbi = dScoreMax + (edge->fScoreAM+edge->fInsertionPenalty)*m_fAMScalingFactor;
	edge->bTouched = true;	
}


// reverse viterbi: for each edge keep successor edge and path score from the edge
void HypothesisLattice::viterbiReverse(LEdge *edge) {

	if (edge->nodeNext == m_lnodeFinal) {
		edge->dScoreViterbiReverse = 0.0;
		edge->bTouched = true;
		edge->edgeNextAux = NULL;
		return;
   }
   
	// (1) compute the reverse Viterbi score of all the successors (only those that are left)
	double dScoreMax = -DBL_MAX;		// accumulated score from successors
	double dAux = 0.0;
	LEdge *edgeSucc = edge->nodeNext->edgeNext;
	while(edgeSucc) {
		
		// compute the reverse Viterbi score of the successor if not available
		if (edgeSucc->bTouched == false) {
			viterbiReverse(edgeSucc);	
		}
		assert(edgeSucc->bTouched);
		assert(edgeSucc->dScoreViterbiReverse != -DBL_MAX);
		
		// keep best score and best successor edge
		dAux = edgeSucc->dScoreViterbiReverse + edgeSucc->fScoreLM*m_fLMScalingFactor +
			(edgeSucc->fScoreAM+edgeSucc->fInsertionPenalty)*m_fAMScalingFactor;
		if (dAux > dScoreMax) {
			dScoreMax = dAux;
			edge->edgeNextAux = edgeSucc;
		}
		
		edgeSucc = edgeSucc->edgePrev;
	}
	assert(edge->edgeNextAux);	
 
	// (2) compute the reverse Viterbi score of the edge (best score)
	edge->dScoreViterbiReverse = dScoreMax;
	edge->bTouched = true;
}

};	// end-of-namespace
