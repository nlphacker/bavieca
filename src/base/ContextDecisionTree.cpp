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

#include "ContextDecisionTree.h"

// constructor
ContextDecisionTree::ContextDecisionTree(int iFeatureDimensionality, int iCovarianceModelling, PhoneSet *phoneSet, unsigned char iContextModelingOrder, Accumulator *accumulator, Rule **rules, int iRules) {
	
	m_iFeatureDimensionality = iFeatureDimensionality;
	m_iCovarianceModelling = iCovarianceModelling;
	m_phoneSet = phoneSet;
	m_iContextModelingOrder = iContextModelingOrder;
	m_dtnodes = NULL;
	m_iNodes = 0;
	m_iLeaves = 0;
	
	// create the root node
	m_dtnodeRoot = new DTNode;
	m_dtnodeRoot->dtnodeYes = NULL;
	m_dtnodeRoot->dtnodeNo = NULL;
	m_dtnodeRoot->dtnodeMerged = NULL;
	m_dtnodeRoot->bIsMerged = false;
	m_dtnodeRoot->accumulator = accumulator;
	m_dtnodeRoot->rule = NULL;
	m_dtnodeRoot->rules = new Rule*[iRules];
	m_dtnodeRoot->iRules = iRules;
	for(int i=0 ; i<iRules ; ++i) {
		m_dtnodeRoot->rules[i] = rules[i];
	}	
	m_dtnodeRoot->iHMMState = -1;	
	m_iNodes++;
	
	// compute the occupation
	m_dOccupation = 0.0;
	Accumulator *accumulatorAux = m_dtnodeRoot->accumulator;
	while(accumulatorAux != NULL) {	
		m_dOccupation += accumulatorAux->getOccupation();
		accumulatorAux = accumulatorAux->getNext();
	}
	
	// allocate memory for the estimation statistics
	m_bEstimationDataAllocated = true;
	
	m_dMeanCluster = new double[m_iFeatureDimensionality];
	m_dCovarianceCluster = new double[m_iFeatureDimensionality];
	m_dObservationCluster = new double[m_iFeatureDimensionality];
	m_dObservationSquareCluster = new double[m_iFeatureDimensionality];
	
	m_dMeanClusterNo = new double[m_iFeatureDimensionality];
	m_dCovarianceClusterNo = new double[m_iFeatureDimensionality];
	m_dObservationClusterNo = new double[m_iFeatureDimensionality];
	m_dObservationSquareClusterNo = new double[m_iFeatureDimensionality];
	
	return;
}

// constructor
ContextDecisionTree::ContextDecisionTree(int iFeatureDimensionality, int iCovarianceModelling, PhoneSet *phoneSet, unsigned char iContextModelingOrder) {

	m_iFeatureDimensionality = iFeatureDimensionality;
	m_iCovarianceModelling = iCovarianceModelling;	
	m_phoneSet = phoneSet;
	m_iContextModelingOrder = iContextModelingOrder;
	m_dtnodeRoot = NULL;
	m_dtnodes = NULL;
	m_iNodes = 0;	
	m_iLeaves = 0;
	m_iFeatureDimensionality = iFeatureDimensionality;	
	
	m_bEstimationDataAllocated = false;
	
	return;
}

// constructor
ContextDecisionTree::ContextDecisionTree() {

	m_iFeatureDimensionality = -1;
	m_phoneSet = NULL;
	m_iContextModelingOrder = 0;
	m_dtnodeRoot = NULL;
	m_dtnodes = NULL;
	m_iNodes = 0;	
	m_iLeaves = 0;
	
	m_bEstimationDataAllocated = false;
	
	return;
}

// copy constructor (creates a replica of the object, for tree traversal only)
ContextDecisionTree::ContextDecisionTree(ContextDecisionTree *contextDecisionTreeSeed) {
	
	m_iFeatureDimensionality = contextDecisionTreeSeed->m_iFeatureDimensionality;
	m_iCovarianceModelling = contextDecisionTreeSeed->m_iCovarianceModelling;
	m_phoneSet = contextDecisionTreeSeed->m_phoneSet;
	m_iContextModelingOrder = contextDecisionTreeSeed->m_iContextModelingOrder;

	// get the number of nodes
	m_iNodes = contextDecisionTreeSeed->m_iNodes;
	m_iLeaves = contextDecisionTreeSeed->m_iLeaves;
	
	// allocate memory for the nodes
	m_dtnodes = new DTNode[m_iNodes];
	memcpy(m_dtnodes,contextDecisionTreeSeed->m_dtnodes,m_iNodes*sizeof(DTNode));
		
	// reset unused fields
	for(int i=0 ; i < m_iNodes ; ++i) {
		m_dtnodes[i].rules = NULL;
		m_dtnodes[i].iRules = -1;
		m_dtnodes[i].dOccupation = -1;
		m_dtnodes[i].accumulator = NULL;
		m_dtnodes[i].dtnodeMerged = NULL;
		m_dtnodes[i].bIsMerged = false;
	}
	
	// replicate the rule objects
	
	// recreate the trees
	for(int i=0 ; i < m_iNodes ; ++i) {
		
	}
	
	
	// TODO this needs to be completed...
	
	m_bEstimationDataAllocated = false;

	return;
}

// destructor
ContextDecisionTree::~ContextDecisionTree() {

	// tree loaded from disk
	if (m_dtnodes != NULL) {
		for(int i=0 ; i < m_iNodes ; ++i) {
			delete m_dtnodes[i].rule;
		}
		delete [] m_dtnodes;
	} 
	// tree created from scratch
	else {
		destroyTree(m_dtnodeRoot);
	}
	
	// release memory
	if (m_bEstimationDataAllocated == true) {
	
		delete [] m_dMeanCluster;
		delete [] m_dCovarianceCluster;
		delete [] m_dObservationCluster;
		delete [] m_dObservationSquareCluster;
		
		delete [] m_dMeanClusterNo;
		delete [] m_dCovarianceClusterNo;
		delete [] m_dObservationClusterNo;
		delete [] m_dObservationSquareClusterNo;	
	}
	
	return;	
}

// load from file
ContextDecisionTree *ContextDecisionTree::load(FileInput &file, int iFeatureDimensionality, 
	int iCovarianceModelling, PhoneSet *phoneSet, PhoneticRulesManager *phoneticRulesManager, 
	unsigned char iContextModelingOrder) {

	VPhoneticRule *vPhoneticRule = phoneticRulesManager->getRules();

	ContextDecisionTree *contextDecisionTree = new ContextDecisionTree(iFeatureDimensionality,iCovarianceModelling,phoneSet,iContextModelingOrder);

	// read the context modeling order
	unsigned char iOrder;
 	IOBase::read(file.getStream(),&iOrder);	
	assert(iContextModelingOrder == iOrder);

	// read the number of nodes
 	IOBase::read(file.getStream(),&contextDecisionTree->m_iNodes);	
	assert(contextDecisionTree->m_iNodes >= 1);
	
	contextDecisionTree->m_iLeaves = 0;
	
	// allocate memory for the nodes
	contextDecisionTree->m_dtnodes = new DTNode[contextDecisionTree->m_iNodes];
	
	// read the nodes one by one
	int iIndex = 0;
	int iIndexChildYes = -1;
	int iIndexChildNo = -1;	
	int iHMMState = -1;
	for(int i=0 ; i < contextDecisionTree->m_iNodes ; ++i) {
		
		// load the children indices
	 	IOBase::read(file.getStream(),&iIndexChildYes);	
 		IOBase::read(file.getStream(),&iIndexChildNo);	
		// leaf node
		if (iIndexChildYes == -1) {
			assert(iIndexChildNo == -1);
			
			// load the HMM-state index
	 		IOBase::read(file.getStream(),&iHMMState);	
	 		assert(iHMMState >= 0);
			
			contextDecisionTree->m_dtnodes[i].dtnodeYes = NULL;	
			contextDecisionTree->m_dtnodes[i].dtnodeNo = NULL;
			contextDecisionTree->m_dtnodes[i].iHMMState = iHMMState;
			contextDecisionTree->m_dtnodes[i].rule = NULL;
			
			contextDecisionTree->m_iLeaves++;
		} 
		// internal node
		else {

			assert(iIndexChildNo != -1);
			assert((iIndexChildYes >= 0) && (iIndexChildYes < contextDecisionTree->m_iNodes));
			assert((iIndexChildNo >= 0) && (iIndexChildNo < contextDecisionTree->m_iNodes));
			
			contextDecisionTree->m_dtnodes[i].dtnodeYes = &contextDecisionTree->m_dtnodes[iIndexChildYes];	
			contextDecisionTree->m_dtnodes[i].dtnodeNo = &contextDecisionTree->m_dtnodes[iIndexChildNo];
			contextDecisionTree->m_dtnodes[i].iHMMState = -1; 
			contextDecisionTree->m_dtnodes[i].rule = new Rule;
			
			// load the rule	
			
			// rule type
			IOBase::read(file.getStream(),&contextDecisionTree->m_dtnodes[i].rule->iType);
			// rule identity
			switch(contextDecisionTree->m_dtnodes[i].rule->iType) {
				// left context
				case RULE_TYPE_PHONETIC_PHONE_LEFT:
				// central phone
				case RULE_TYPE_PHONETIC_PHONE: 
				// right context
				case RULE_TYPE_PHONETIC_PHONE_RIGHT: {
					if (contextDecisionTree->m_dtnodes[i].rule->iType != RULE_TYPE_PHONETIC_PHONE) {
						// read the context position
						IOBase::read(file.getStream(),&contextDecisionTree->m_dtnodes[i].rule->iContextPosition);
					}
					// read the phonetic rule index
					int iRule = -1;
					IOBase::read(file.getStream(),&iRule);
					assert((iRule >= 0) && (iRule < phoneticRulesManager->getRules()->size()));
					// get the rule
					contextDecisionTree->m_dtnodes[i].rule->phoneticRule = phoneticRulesManager->getRule(iRule);
					if (contextDecisionTree->m_dtnodes[i].rule->phoneticRule == NULL) {
						return NULL;
					}
					break;
				}
				// within-word position
				case RULE_TYPE_POSITION: {
					contextDecisionTree->m_dtnodes[i].rule->phoneticRule = NULL;
					// read the position
					IOBase::read(file.getStream(),&contextDecisionTree->m_dtnodes[i].rule->iPosition);
					assert(HMMState::isPositionValid(contextDecisionTree->m_dtnodes[i].rule->iPosition));
					break;
				}
				// HMM-state
				case RULE_TYPE_STATE: {
					contextDecisionTree->m_dtnodes[i].rule->phoneticRule = NULL;
					// read the state
					IOBase::read(file.getStream(),&contextDecisionTree->m_dtnodes[i].rule->iState);
					assert(HMMState::isStateValid(contextDecisionTree->m_dtnodes[i].rule->iState));
					break;
				}
				default: {
					assert(0);
				}
			}			
		}
	}	
	// root node?
	contextDecisionTree->m_dtnodeRoot = &contextDecisionTree->m_dtnodes[0];

	return contextDecisionTree;
}

// store to file
void ContextDecisionTree::store(FileOutput &file) {

	IOBase::write(file.getStream(),m_iContextModelingOrder);	
	IOBase::write(file.getStream(),m_iNodes);	
	
	LDTNode lDTNode;	
	lDTNode.push_back(m_dtnodeRoot);
	int iIndex = 0;
	int iIndexChildYes = -1;
	int iIndexChildNo = -1;
	while(lDTNode.empty() == false) {
	
		DTNode *node = lDTNode.front();
		lDTNode.pop_front();
		
		if (node->dtnodeYes != NULL) {
			lDTNode.push_back(node->dtnodeYes);
			++iIndex;
			iIndexChildYes = iIndex;
		} else {
			iIndexChildYes = -1;
		}
		if (node->dtnodeNo != NULL) {
			lDTNode.push_back(node->dtnodeNo);
			++iIndex;
			iIndexChildNo = iIndex;
		} else {
			iIndexChildNo = -1;
		}
		
		// write the children index
		IOBase::write(file.getStream(),iIndexChildYes);	
		IOBase::write(file.getStream(),iIndexChildNo);	
		
		// internal node: write the phonetic rule 
		if (node->rule != NULL) {
			// rule type
			IOBase::write(file.getStream(),node->rule->iType);	
			// rule identity
			switch(node->rule->iType) {
				// left context
				case RULE_TYPE_PHONETIC_PHONE_LEFT: 
				// central phone
				case RULE_TYPE_PHONETIC_PHONE: 
				// right context
				case RULE_TYPE_PHONETIC_PHONE_RIGHT: {
					if (node->rule->iType != RULE_TYPE_PHONETIC_PHONE) {
						// write the context position
						IOBase::write(file.getStream(),node->rule->iContextPosition);	
					}
					// write the phonetic rule index
					IOBase::write(file.getStream(),node->rule->phoneticRule->iRule);	
					break;
				}
				// within-word position
				case RULE_TYPE_POSITION: {
					// write the position
					IOBase::write(file.getStream(),node->rule->iPosition);
					break;
				}
				// HMM-state
				case RULE_TYPE_STATE: {
					// write the state
					IOBase::write(file.getStream(),node->rule->iState);
					break;
				}
				default: {
					assert(0);
				}
			}
		}
		// leaf node: write the HMM-state
		else {
			assert(node->iHMMState != -1);	
			IOBase::write(file.getStream(),node->iHMMState);
		}
	}	
	assert(iIndex+1 == m_iNodes);
}

// return the HMM-state associated to the given context-dependent phone
int ContextDecisionTree::getHMMIndex(unsigned char *iPhoneLeft, unsigned char iPhone, 
	unsigned char *iPhoneRight, unsigned char iPosition, 
	unsigned char iState) {
		
	DTNode *node = m_dtnodeRoot;
	while(node->rule != NULL) {
		bool bAnswer = true;	
		switch(node->rule->iType) {	
			case RULE_TYPE_PHONETIC_PHONE_LEFT: {
				assert(iPhoneLeft[node->rule->iContextPosition] != UCHAR_MAX);
				bAnswer = node->rule->phoneticRule->bPhone[iPhoneLeft[node->rule->iContextPosition]];
				break;
			}
			case RULE_TYPE_PHONETIC_PHONE: {
				bAnswer = node->rule->phoneticRule->bPhone[iPhone];
				break;
			}
			case RULE_TYPE_PHONETIC_PHONE_RIGHT: {
				assert(iPhoneRight[node->rule->iContextPosition] != UCHAR_MAX);
				bAnswer = node->rule->phoneticRule->bPhone[iPhoneRight[node->rule->iContextPosition]];
				break;
			}
			case RULE_TYPE_POSITION: {	
				bAnswer = (node->rule->iPosition == iPosition);
				break;
			}
			case RULE_TYPE_STATE: {
				bAnswer = (node->rule->iState == iState);
				break;
			}
			default: {
				assert(0);	
				break;
			}
		}
		if (bAnswer == true) {
			node = node->dtnodeYes;
		} else {
			node = node->dtnodeNo;
		}	
	}
	
	// make sure we got to a leaf
	assert(node->iHMMState != -1);
	
	return node->iHMMState;	
}

// cluster a node (if possible)
void ContextDecisionTree::clusterNode(DTNode *node, float fMinimumClusterOccupation, float fMinimumLikelihoodGainClustering) {

	// iterate through the phonetic rules and pick the one that produces the higher likelihood gain
	bool bRuleApplied = false;
	double dLikelihoodRuleYes;
	double dLikelihoodRuleNo;
	double dOccupationRuleYes;
	double dOccupationRuleNo;
	double dOccupationRuleBestYes = -1.0;
	double dOccupationRuleBestNo = -1.0;
	double dLikelihoodRuleBest = -DBL_MAX;
	double dLikelihoodRuleBestYes = -DBL_MAX;
	double dLikelihoodRuleBestNo = -DBL_MAX;
	Rule *ruleBest = NULL;
	int iRuleBest = -1;
	
	// (1) apply all the available rules and keep the best rule
	for(int iRule = 0 ; iRule < node->iRules ; ++iRule) {
		if (node->rules[iRule] == NULL) {
			continue;
		}
		// positive answer
		bool bYes = computeLikelihoodRule(node,node->rules[iRule],true,&dLikelihoodRuleYes,&dOccupationRuleYes,fMinimumClusterOccupation);
		if (bYes == false) {
			continue;
		}
		// negative answer
		bool bNo = computeLikelihoodRule(node,node->rules[iRule],false,&dLikelihoodRuleNo,&dOccupationRuleNo,fMinimumClusterOccupation);
		if (bNo == false) {
			continue;
		}
		/*bool b = computeLikelihoodRule(node,node->rules[iRule],&fLikelihoodRuleYes,&fLikelihoodRuleNo,&fOccupationRuleYes,&fOccupationRuleNo,fMinimumClusterOccupation,fCovarianceFlooringRatio,fCovarianceGlobal);
		if (b == false) {
			continue;
		}*/
		//printf("rule likelihood: %.4f (%.4f %.4f)\n",node->fLikelihood,fLikelihoodRuleYes,fLikelihoodRuleNo);
		// keep the best rule
		if ((dLikelihoodRuleYes + dLikelihoodRuleNo) > dLikelihoodRuleBest) {	
			dLikelihoodRuleBest = (dLikelihoodRuleYes + dLikelihoodRuleNo);
			dLikelihoodRuleBestYes = dLikelihoodRuleYes;
			dLikelihoodRuleBestNo = dLikelihoodRuleNo;
			dOccupationRuleBestYes = dOccupationRuleYes;
			dOccupationRuleBestNo = dOccupationRuleNo;
			ruleBest = node->rules[iRule];
			iRuleBest = iRule;
		}
	}
	
	// (2) apply the best split if it is good enough and continue with the clustering, otherwise end it
	
	// if no rule could be applied, end the clustering
	if (ruleBest == NULL) {
		return;
	}
	// compute the likelihood gain resulting from applying the best rule
	double dLikelihoodGain = dLikelihoodRuleBest-node->dLikelihood;
	if (dLikelihoodGain < 0.0) {
		printf("Unexpected: likelihoodGain= %f (%f %f)\n",dLikelihoodGain,dLikelihoodRuleBest,node->dLikelihood);
	}
	//assert(fLikelihoodGain >= 0.0);
	
	// the best rule is applied only if the likelihood gain is high enough
	if (dLikelihoodGain >= fMinimumLikelihoodGainClustering) {
	
		// keep the applied rule
		node->rule = ruleBest;
				
		// create the child nodes
		
		// "yes" child
		node->dtnodeYes = new DTNode;
		node->dtnodeYes->dLikelihood = dLikelihoodRuleBestYes;
		node->dtnodeYes->dOccupation = dOccupationRuleBestYes;		
		node->dtnodeYes->dtnodeYes = NULL;
		node->dtnodeYes->dtnodeNo = NULL;
		node->dtnodeYes->dtnodeMerged = NULL;
		node->dtnodeYes->bIsMerged = false;
		node->dtnodeYes->rule = NULL;
		node->dtnodeYes->rules = new Rule*[node->iRules];
		node->dtnodeYes->iRules = node->iRules;
		node->dtnodeYes->iHMMState = -1;
		m_iNodes++;
		
		// copy the parent rules except the one that was just applied
		for(int i=0 ; i < node->iRules ; ++i) {
			node->dtnodeYes->rules[i] = node->rules[i];
		}
		node->dtnodeYes->rules[iRuleBest] = NULL;
		
		// "no" child
		node->dtnodeNo = new DTNode;
		node->dtnodeNo->dLikelihood = dLikelihoodRuleBestNo;
		node->dtnodeNo->dOccupation = dOccupationRuleBestNo;
		node->dtnodeNo->dtnodeYes = NULL;
		node->dtnodeNo->dtnodeNo = NULL;
		node->dtnodeNo->dtnodeMerged = NULL;
		node->dtnodeNo->bIsMerged = false;
		node->dtnodeNo->rule = NULL;
		node->dtnodeNo->rules = new Rule*[node->iRules];
		node->dtnodeNo->iRules = node->iRules;
		node->dtnodeNo->iHMMState = -1;
		m_iNodes++;
		
		// copy the parent rules except the one that was just applied
		for(int i=0 ; i < node->iRules ; ++i) {
			node->dtnodeNo->rules[i] = node->rules[i];
		}
		node->dtnodeNo->rules[iRuleBest] = NULL;
		
		// distribute the parent's accumulators among the child nodes
		node->dtnodeYes->accumulator = NULL;
		node->dtnodeNo->accumulator = NULL;
		Accumulator *accumulator = node->accumulator;
		Accumulator *accumulatorAux = NULL;
		while(accumulator != NULL) {	
			accumulatorAux = accumulator->getNext();
			if (question(ruleBest,accumulator) == true) {	
				accumulator->setNext(node->dtnodeYes->accumulator);
				node->dtnodeYes->accumulator = accumulator;
			} else {
				accumulator->setNext(node->dtnodeNo->accumulator);
				node->dtnodeNo->accumulator = accumulator;
			}
			accumulator = accumulatorAux;
		}	
		node->accumulator = NULL;	
		
		// continue with the clustering
		clusterNode(node->dtnodeYes,fMinimumClusterOccupation,fMinimumLikelihoodGainClustering);
		clusterNode(node->dtnodeNo,fMinimumClusterOccupation,fMinimumLikelihoodGainClustering);
	}	
	
	//important: two complimentary ways to stop the clustering:
	// (1) when no cluster has enough occupation to be splitted
	// (2) when the likelihood gain of applying rules to any leaf node is too small (i.e. further clustering is unproductive)
	// (3) additionally, when the number of leaves reaches a certain maximum value

	return;
}

// compute the likelihood of a cluster of context dependent units 
bool ContextDecisionTree::computeLikelihoodCluster(DTNode *node, double *dLikelihoodCluster, double *dOccupationCluster, float fCovarianceFlooringRatio, float *fCovarianceGlobal) {

	// (1) compute the mean and covariance of the new cluster
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = 0.0;
		m_dCovarianceCluster[i] = 0.0;
		m_dObservationCluster[i] = 0.0;
		m_dObservationSquareCluster[i] = 0.0;
	}
	m_dOccupationCluster = 0.0;
	
	// collect data from all the accumulators in the node
	Accumulator *accumulator = node->accumulator;
	while(accumulator != NULL) {	
		for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
			m_dObservationCluster[i] += accumulator->getObservation()[i];
			m_dObservationSquareCluster[i] += accumulator->getObservationSquare()[i];
		}
		m_dOccupationCluster += accumulator->getOccupation();
		accumulator = accumulator->getNext();
	}
	// keep the occupation
	node->dOccupation = m_dOccupationCluster;	
	
	// if all the data goes to either side of the rule, it cannot be applied
	assert(m_dOccupationCluster > 0);
	
	// (1) compute the mean
	for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = (float)(m_dObservationCluster[i]/m_dOccupationCluster);
	}
	
	// (2) compute the covariance (using the already computed mean)
	for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dCovarianceCluster[i] = (float)(((m_dOccupationCluster*m_dMeanCluster[i]*m_dMeanCluster[i])+m_dObservationSquareCluster[i]-(2*m_dObservationCluster[i]*m_dMeanCluster[i]))/m_dOccupationCluster);	
	}
	
	// (3) compute the likelihood of the new cluster
	double dDeterminant = 0.0;
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		assert(m_dCovarianceCluster[i] > 0.0);
		dDeterminant += log(m_dCovarianceCluster[i]);
	}
	assert(finite(dDeterminant));
	*dLikelihoodCluster = ((-0.5*m_dOccupationCluster)*(dDeterminant+(m_iFeatureDimensionality*(log(2.0*PI_NUMBER)+1.0))));	
	assert(finite(*dLikelihoodCluster));
	
	*dOccupationCluster = m_dOccupationCluster;
	
	return true;
}

// compute the cluster likelihood
double ContextDecisionTree::computeLikelihoodCluster(double dOccupation, double *dCovariance) {

	double dDeterminant = 0.0;
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		assert(dCovariance[i] > 0.0);
		dDeterminant += log(dCovariance[i]);
	}
	assert(finite(dDeterminant));
	double dLikelihood = -0.5*dOccupation*(dDeterminant+(m_iFeatureDimensionality*(log(2.0*PI_NUMBER)+1.0)));	
	assert(finite(dLikelihood));

	return dLikelihood;
}

// compute the cluster likelihood
double ContextDecisionTree::computeLikelihoodCluster(double dOccupation, float *fCovariance) {

	double dDeterminant = 0.0;
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		assert(fCovariance[i] > 0.0);
		dDeterminant += log(fCovariance[i]);
	}
	assert(finite(dDeterminant));
	double dLikelihood = -0.5*dOccupation*(dDeterminant+(m_iFeatureDimensionality*(log(2.0*PI_NUMBER)+1.0)));	
	assert(finite(dLikelihood));

	return dLikelihood;
}


// compute the likelihood of a cluster of nphones given the parent cluster, the rule and the decission made
// note: return false if after applying the rule:
// - there is not enough data to robustly estimate the parameters (minimum occupation count)
// - there is too much data resulting from the answer to the question, so the opposite answer wont produce enough data
bool ContextDecisionTree::computeLikelihoodRule(DTNode *node, Rule *rule, bool bAnswer, double *dLikelihoodCluster, double *dOccupationCluster, float fMinimumClusterOccupation) {

	// (1) accumulate data for the new cluster
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dObservationCluster[i] = 0.0;
		m_dObservationSquareCluster[i] = 0.0;
	}
	m_dOccupationCluster = 0.0;
	
	// rule: only accumulators for which the rule produces the given answer
	assert(rule != NULL);
	Accumulator *accumulator = node->accumulator;
	while(accumulator != NULL) {	
		if (question(rule,accumulator) == bAnswer) {
			for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
				m_dObservationCluster[i] += accumulator->getObservation()[i];
				m_dObservationSquareCluster[i] += accumulator->getObservationSquare()[i];
			}
			m_dOccupationCluster += accumulator->getOccupation();
		}	
		accumulator = accumulator->getNext();
	}	
	*dOccupationCluster = m_dOccupationCluster;
	// check minimum and maximum occupation 
	if ((m_dOccupationCluster < fMinimumClusterOccupation) || 
		(m_dOccupationCluster > (node->dOccupation-fMinimumClusterOccupation))) {
		return false;
	}
	
	// (2) compute the mean and covariance of the new cluster
		
	// initialize variables
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = 0.0;
		m_dCovarianceCluster[i] = 0.0;
	}	
	
	// compute the mean
	for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = m_dObservationCluster[i]/m_dOccupationCluster;
	}
	
	// compute the covariance (using the already computed mean)
	for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dCovarianceCluster[i] = ((m_dOccupationCluster*m_dMeanCluster[i]*m_dMeanCluster[i])+m_dObservationSquareCluster[i]-(2*m_dObservationCluster[i]*m_dMeanCluster[i]))/m_dOccupationCluster;	
	}
	
	// (3) compute the likelihood of the new cluster
	*dLikelihoodCluster = computeLikelihoodCluster(m_dOccupationCluster,m_dCovarianceCluster);
	
	return true;
}


// compute the likelihood of the clusters resulting from applying a question to the given cluster
// note: return false if after applying the rule:
// - there is not enough data to robustly estimate the parameters (minimum occupation count)
// - there is too much data resulting from the answer to the question, so the opposite answer wont produce enough data
bool ContextDecisionTree::computeLikelihoodRule(DTNode *node, Rule *rule, double *dLikelihoodClusterYes, double *dLikelihoodClusterNo, double *dOccupationClusterYes, double *dOccupationClusterNo, float fMinimumClusterOccupation) {

	// initialize variables
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dObservationCluster[i] = 0.0;
		m_dObservationSquareCluster[i] = 0.0;
	}	
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dObservationClusterNo[i] = 0.0;
		m_dObservationSquareClusterNo[i] = 0.0;
	}	
	m_dOccupationCluster = 0.0;
	m_dOccupationClusterNo = 0.0;
	
	// rule: only accumulators for which the rule produces the given answer
	assert(rule != NULL);
	Accumulator *accumulator = node->accumulator;
	while(accumulator != NULL) {	
		// positive answer
		if (question(rule,accumulator) == true) {
			for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
				m_dObservationCluster[i] += accumulator->getObservation()[i];
				m_dObservationSquareCluster[i] += accumulator->getObservationSquare()[i];
			}
			m_dOccupationCluster += accumulator->getOccupation();
		} 
		// negative answer
		else {
			for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
				m_dObservationClusterNo[i] += accumulator->getObservation()[i];
				m_dObservationSquareClusterNo[i] += accumulator->getObservationSquare()[i];
			}
			m_dOccupationClusterNo += accumulator->getOccupation();
		}
		accumulator = accumulator->getNext();
	}	
	*dOccupationClusterYes = m_dOccupationCluster;
	*dOccupationClusterNo = m_dOccupationClusterNo;
	// check minimum occupation
	if ((m_dOccupationCluster < fMinimumClusterOccupation) || (m_dOccupationClusterNo < fMinimumClusterOccupation)) {
		return false;
	}	
	
	// if all the data goes to either side of the rule, it cannot be applied
	assert(m_dOccupationCluster >= fMinimumClusterOccupation);
	assert(m_dOccupationClusterNo >= fMinimumClusterOccupation);
	
	// initialize variables
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = 0.0;
		m_dCovarianceCluster[i] = 0.0;
	}	
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dMeanClusterNo[i] = 0.0;
		m_dCovarianceClusterNo[i] = 0.0;
	}	
	
	// (1) compute the mean
	for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = m_dObservationCluster[i]/m_dOccupationCluster;
		m_dMeanClusterNo[i] = m_dObservationClusterNo[i]/m_dOccupationClusterNo;
	}
	
	// (2) compute the covariance (using the already computed mean)
	for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dCovarianceCluster[i] = ((m_dOccupationCluster*m_dMeanCluster[i]*m_dMeanCluster[i])+m_dObservationSquareCluster[i]-(2*m_dObservationCluster[i]*m_dMeanCluster[i]))/m_dOccupationCluster;	
		m_dCovarianceClusterNo[i] = ((m_dOccupationClusterNo*m_dMeanClusterNo[i]*m_dMeanClusterNo[i])+m_dObservationSquareClusterNo[i]-(2*m_dObservationClusterNo[i]*m_dMeanClusterNo[i]))/m_dOccupationClusterNo;	
	}
	
	// (3) compute the likelihood of the new cluster
	*dLikelihoodClusterYes = computeLikelihoodCluster(m_dOccupationCluster,m_dCovarianceCluster);
	*dLikelihoodClusterNo = computeLikelihoodCluster(m_dOccupationClusterNo,m_dCovarianceClusterNo);
	
	return true;
}


// compute the likelihood of a decission tree by summing up the likelihood in the leaves
double ContextDecisionTree::computeTreeLikelihood(DTNode *node) {
	
	if (node->dtnodeYes == NULL) {
		assert(node->dtnodeNo == NULL);
		return node->dLikelihood;
	} else {
		return (computeTreeLikelihood(node->dtnodeYes) + computeTreeLikelihood(node->dtnodeNo));
	}
}

// return the number of leaves in the tree (including those leaves merged to another leaf)
int ContextDecisionTree::countTreeLeaves(DTNode *node) {
	
	if (node->dtnodeYes == NULL) {
		assert(node->dtnodeNo == NULL);
		return 1;
	} else {
		return (countTreeLeaves(node->dtnodeYes) + countTreeLeaves(node->dtnodeNo));
	}
}

// return the number of leaves in the tree (excluding those leaves merged to another leaf)
int ContextDecisionTree::countTreeLeavesUnique(DTNode *node) {
	
	if (node->rule == NULL) {
		if (node->dtnodeMerged != NULL) {	
			return 1;
		} else {
			return 0;
		}
	} else {
		return (countTreeLeaves(node->dtnodeYes) + countTreeLeaves(node->dtnodeNo));
	}
}


// get the tree leaves
void ContextDecisionTree::getTreeLeaves(DTNode *node, LDTNode &lLeaves) {

	if (node->dtnodeYes == NULL) {
		assert(node->dtnodeNo == NULL);
		lLeaves.push_back(node);
	} else {
		getTreeLeaves(node->dtnodeYes,lLeaves);
		getTreeLeaves(node->dtnodeNo,lLeaves);
	}

	return;
}

// print tree leaves
void ContextDecisionTree::printTreeLeaves(DTNode *node) {
	
	if (node->rule != NULL) {
		printTreeLeaves(node->dtnodeYes);	
		printTreeLeaves(node->dtnodeNo);	
	} else {
		printf("leaf: %d\n",node->iHMMState);
	}
}

// check tree consistency
bool ContextDecisionTree::checkConsistency(DTNode *node) {

	if (node->rule != NULL) {
		if (node->rule->iType > RULE_TYPE_STATE) {
			return false;
		}
		if (node->iHMMState != -1) {
			return false;
		}
		if (checkConsistency(node->dtnodeYes) == false) {
			return false;
		}
		if (checkConsistency(node->dtnodeNo) == false) {
			return false;
		}
	} else {
		if (node->iHMMState == -1) {
			return false;
		}
	}
	
	return true;
}

// compute the likelihood of a decission tree by summing up the likelihood in the leaves
int ContextDecisionTree::countTreeLeavesOccupancy(DTNode *node) {
	
	if (node->dtnodeYes == NULL) {
		assert(node->dtnodeNo == NULL);
		return 1;
	} else {
		return (countTreeLeaves(node->dtnodeYes) + countTreeLeaves(node->dtnodeNo));
	}
}

// destroy a decission tree
void ContextDecisionTree::destroyTree(DTNode *node) {
	
	if (node->rule != NULL) {
		destroyTree(node->dtnodeYes);
		destroyTree(node->dtnodeNo);
	} 
	delete [] node->rules;
	delete node;
}

// compute the likelihood of the cluster resulting from mergining two leaves
double ContextDecisionTree::computeLikelihoodFromMerging(DTNode *nodeA, DTNode *nodeB, float fCovarianceFlooringRatio, float *fCovarianceGlobal) {

	// initialize the accumulators used to compute the mean and covariance of the new cluster
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = 0.0;
		m_dCovarianceCluster[i] = 0.0;
		m_dObservationCluster[i] = 0.0;
		m_dObservationSquareCluster[i] = 0.0;
	}
	m_dOccupationCluster = 0.0;
	
	// accumulate data from nodeA	
	Accumulator *accumulator = nodeA->accumulator;
	while(accumulator) {	
		for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
			m_dObservationCluster[i] += accumulator->getObservation()[i];
			m_dObservationSquareCluster[i] += accumulator->getObservationSquare()[i];
		}
		m_dOccupationCluster += accumulator->getOccupation();
		accumulator = accumulator->getNext();
	}
	// accumulate data from nodeB	
	accumulator = nodeB->accumulator;
	while(accumulator) {	
		for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
			m_dObservationCluster[i] += accumulator->getObservation()[i];
			m_dObservationSquareCluster[i] += accumulator->getObservationSquare()[i];
		}
		m_dOccupationCluster += accumulator->getOccupation();
		accumulator = accumulator->getNext();
	}
	
	// (1) compute the mean
	for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dMeanCluster[i] = m_dObservationCluster[i]/m_dOccupationCluster;
	}
	
	// (2) compute the covariance (using the already computed mean)
	for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
		m_dCovarianceCluster[i] = ((m_dOccupationCluster*m_dMeanCluster[i]*m_dMeanCluster[i])+m_dObservationSquareCluster[i]-(2*m_dObservationCluster[i]*m_dMeanCluster[i]))/m_dOccupationCluster;	
	}
	
	// (2) compute the likelihood of the new cluster
	double dLikelihoodCluster = computeLikelihoodCluster(m_dOccupationCluster,m_dCovarianceCluster);

	return dLikelihoodCluster;
}

// compact the tree leaves by merging those leaves that when merged the likelihood decrease is below the minimum splitting gain used for clustering (return the resulting number of leaves)
int ContextDecisionTree::compactTreeLeaves(float fMinimumLikelihoodGainClustering, float fCovarianceFlooringRatio, float *fCovarianceGlobal) {
	
	// get an array with the tree leaves
	LDTNode lLeaves;
	getTreeLeaves(m_dtnodeRoot,lLeaves);
	
	int iLeavesStart = lLeaves.size();
	
	// for each leaf in the tree check if there is any other leaf to merge it with
	for(LDTNode::iterator it = lLeaves.begin() ; it != lLeaves.end() ; ) {
		// skip the leaf node if it was already merged 
		if ((*it)->dtnodeMerged != NULL) {
			continue;
		}
		bool bRestart = false;
		for(LDTNode::iterator jt = lLeaves.begin() ; jt != lLeaves.end() ; ++jt) {
			// skip the self node or nodes that are already merged
			if ((it == jt) || ((*jt)->dtnodeMerged != NULL)) {
				continue;
			}
			// compute the likelihood of the cluster resulting from merging both leaves
			double dLikelihoodCluster = computeLikelihoodFromMerging(*it,*jt,fCovarianceFlooringRatio,fCovarianceGlobal);
			double dLikelihoodDecrease = ((*it)->dLikelihood + (*jt)->dLikelihood) - dLikelihoodCluster;
			//assert(fLikelihoodDecrease >= 0.0);
			// cluster the leaves together if the likelihood decrease is smaller enough
			if (dLikelihoodDecrease < fMinimumLikelihoodGainClustering) {
				//printf("compacting is possible\n");
				// mark one leaf as merged and move all its accumulators to the other leaf
				(*jt)->dtnodeMerged = (*it);
				assert((*it)->dtnodeMerged == NULL); 
				(*it)->bIsMerged = true;
				(*jt)->bIsMerged = true;
				(*jt)->dLikelihood = 0.0;
				(*it)->dLikelihood = dLikelihoodCluster;
				// find the last accumulator and append the other node accumulators to that one
				Accumulator *accumulatorAux = (*it)->accumulator;
				while(accumulatorAux->getNext() != NULL) {
					accumulatorAux = accumulatorAux->getNext();
				}
				accumulatorAux->setNext((*jt)->accumulator);
				(*jt)->accumulator = NULL;
				lLeaves.erase(jt);
				// since the original leave node has change it is necessary to start from the beginning
				bRestart = true;
				break;
			}	
		}		
		// check if it is necessary to restart 
		if (bRestart == false) {
			++it;
		}
	}	
	
	// sanity checks
	//printf("sanity checks\n");
	int iLeaf1 = 0;
	for(LDTNode::iterator it = lLeaves.begin() ; it != lLeaves.end() ; ++it, ++iLeaf1) {
		// skip the leaf node if it was already merged 
		if ((*it)->dtnodeMerged != NULL) {
			continue;
		}
		int iLeaf2 = 0;
		for(LDTNode::iterator jt = lLeaves.begin() ; jt != lLeaves.end() ; ++jt, ++iLeaf2) {
			// skip the self node or nodes that are already merged
			if ((it == jt) || ((*jt)->dtnodeMerged != NULL)) {
				continue;
			}
			double dLikelihoodCluster = computeLikelihoodFromMerging(*it,*jt,fCovarianceFlooringRatio,fCovarianceGlobal);
			double dLikelihoodDecrease = ((*it)->dLikelihood + (*jt)->dLikelihood) - dLikelihoodCluster;
			assert(dLikelihoodDecrease >= fMinimumLikelihoodGainClustering);
			//printf("(%3d,%3d) -> %10.4f\n",iLeaf1,iLeaf2,fLikelihoodDecrease);
		}
	}	
	
	int iLeavesEnd = lLeaves.size();	
	assert(iLeavesStart >= iLeavesEnd);

	return iLeavesEnd;
}

// create an HMM-state when there is no data associated to it (default to globla distribution)
void ContextDecisionTree::createHMMStateNoData(unsigned char iPhone, unsigned char iState, HMMState **hmmStates, int *iLeavesSeen, int *iUnitsCreated, float *fMeanGlobal, float *fCovarianceGlobal) {
	
	// create the HMM-state	
	m_dtnodeRoot->iHMMState = *iUnitsCreated;
	hmmStates[*iUnitsCreated] = new HMMState(m_iFeatureDimensionality,m_iCovarianceModelling,m_phoneSet,iPhone,iState,-1,1,*iUnitsCreated);
	Gaussian *gaussian = hmmStates[m_dtnodeRoot->iHMMState]->getGaussian(0);
	
	// default mean and covariance to global distribution
	for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
		gaussian->fMean[i] = fMeanGlobal[i];
		gaussian->fCovariance[i] = fCovarianceGlobal[i];
	}	
	gaussian->fWeight = 1.0;
	
	// sanity check (compute the likelihood of the HMM-state)
	double dDeterminant = 0.0;
	for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
		dDeterminant += log(gaussian->fCovariance[i]);
	}
	assert(finite(dDeterminant));
	
	// increment the number of leaves and accumulate likelihood
	(*iUnitsCreated)++;
	(*iLeavesSeen)++;
	
	return;
}

// create the new set of HMM-states from the leaves of the tree (return the number of physical nphones created)
void ContextDecisionTree::createClusteredContextDependentUnitsFromLeaves(DTNode *node, HMMState **hmmStates,
	int *iLeavesSeen, int *iUnitsCreated, double *dLikelihood, float fCovarianceFlooringRatio) {

	// if it is not a leaf, continue
	if (node->rule != NULL) {
		createClusteredContextDependentUnitsFromLeaves(node->dtnodeYes,hmmStates,iLeavesSeen,iUnitsCreated,dLikelihood,fCovarianceFlooringRatio);
		createClusteredContextDependentUnitsFromLeaves(node->dtnodeNo,hmmStates,iLeavesSeen,iUnitsCreated,dLikelihood,fCovarianceFlooringRatio);
		return;
	}	
	// leaf node
	else {
	
		// check if the node was merged and the HMM-state already exists 
		if ((node->bIsMerged == true) && (node->iHMMState != -1)) {
			// increment the number of leaves and accumulate likelihood
			(*iLeavesSeen)++;	
			return;
		}
	
		// create the new HMM-state from accumulators in the node (or using the HMM-state of a merged node)	
		assert(node->iHMMState == -1);
		Gaussian *gaussian = NULL;
		
		// check if the leaf node was merged to another leaf node	
		if (node->dtnodeMerged != NULL) {
			
			DTNode *dtnodeAux = node->dtnodeMerged;
			while(dtnodeAux->dtnodeMerged != NULL) {
				dtnodeAux = dtnodeAux->dtnodeMerged;
			}
			// check if the HMM-state of the merged leaf already exist
			if (dtnodeAux->iHMMState != -1) {
				node->iHMMState = dtnodeAux->iHMMState;
				// increment the number of leaves (likelihood is only accumulated when a HMM-state is created)
				(*iLeavesSeen)++;
				return;
			} 
			// otherwise create it
			else {
				// create the HMM-state
				unsigned char iPhone = node->accumulator->getPhone();
				unsigned char iState = node->accumulator->getState();
				node->iHMMState = *iUnitsCreated;
				hmmStates[*iUnitsCreated] = new HMMState(m_iFeatureDimensionality,m_iCovarianceModelling,m_phoneSet,iPhone,iState,-1,1,*iUnitsCreated);	
				dtnodeAux->iHMMState = node->iHMMState;
				gaussian = hmmStates[node->iHMMState]->getGaussian(0);	
			}	
		} else {
			// create the HMM-state	
			unsigned char iPhone = node->accumulator->getPhone();
			unsigned char iState = node->accumulator->getState();
			node->iHMMState = *iUnitsCreated;
			hmmStates[*iUnitsCreated] = new HMMState(m_iFeatureDimensionality,m_iCovarianceModelling,m_phoneSet,iPhone,iState,-1,1,*iUnitsCreated);
			gaussian = hmmStates[node->iHMMState]->getGaussian(0);
		}
		assert(node->iHMMState != -1);
		
		// (1) reset accumulators
		for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
			m_dMeanCluster[i] = 0.0;
			m_dCovarianceCluster[i] = 0.0;
		}
		for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
			m_dObservationCluster[i] = 0.0;
			m_dObservationSquareCluster[i] = 0.0;
		}
		m_dOccupationCluster = 0.0;	
		
		// (2) accumulate data
		Accumulator *accumulator = node->accumulator;
		while(accumulator != NULL) {	
			for(int i=0 ; i<m_iFeatureDimensionality ; ++i) {
				m_dObservationCluster[i] += accumulator->getObservation()[i];
				m_dObservationSquareCluster[i] += accumulator->getObservationSquare()[i];
			}
			m_dOccupationCluster += accumulator->getOccupation();
			accumulator = accumulator->getNext();
		}	
	
		// (3) compute the mean
		for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
			gaussian->fMean[i] = (float)(m_dObservationCluster[i]/m_dOccupationCluster);
		}
		
		// (4) compute the covariance (using the already computed mean)
		for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
			gaussian->fCovariance[i] = (float)(((m_dOccupationCluster*gaussian->fMean[i]*gaussian->fMean[i])+m_dObservationSquareCluster[i]-(2*m_dObservationCluster[i]*gaussian->fMean[i]))/m_dOccupationCluster);	
		}
		// set the weight
		gaussian->fWeight = 1.0;
		
		// sanity check (compute the likelihood of the HMM-state)
		double dLikelihoodCluster = computeLikelihoodCluster(m_dOccupationCluster,gaussian->fCovariance);
		
		// increment the number of leaves and accumulate likelihood
		(*iUnitsCreated)++;
		(*iLeavesSeen)++;
		(*dLikelihood) += dLikelihoodCluster;	
		
		return;
	}
}


