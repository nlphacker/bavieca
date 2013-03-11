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


#include "ContextModeling.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"

#include <assert.h>

namespace Bavieca {

// constructor
ContextModeling::ContextModeling(int iFeatureDimensionality, int iCovarianceModelling, PhoneSet *phoneSet, int iContextModelingOrderWW, int iContextModelingOrderCW, const char *strFilePhoneticRules, bool bGlobalClustering, bool bBottomUpMerging, float fMinimumClusterOccupation, float fMinimumLikelihoodGain)
{
	m_iDim = iFeatureDimensionality;
	m_iCovarianceModelling = iCovarianceModelling;
	m_phoneSet = phoneSet;
	m_strFilePhoneticRules = strFilePhoneticRules;
	m_bGlobalClustering = bGlobalClustering;
	m_bBottomUpMerging = bBottomUpMerging;
	m_fMinimumClusterOccupation = fMinimumClusterOccupation;
	m_fMinimumLikelihoodGainClustering = fMinimumLikelihoodGain;	
	m_iPhones = m_phoneSet->size();
	
	// initialize
	m_hmmStatesClustered = NULL;
	m_iContextModelingOrderCW = iContextModelingOrderCW;
	m_iContextModelingOrderWW = iContextModelingOrderWW;
	assert(m_iContextModelingOrderCW%2 == 1);
	assert(m_iContextModelingOrderWW%2 == 1);
	m_iContextSizeCW = (iContextModelingOrderCW-1)/2;
	m_iContextSizeWW = (iContextModelingOrderWW-1)/2;
	assert(m_iContextSizeCW >= 1);
	assert(m_iContextSizeWW >= 1);
	assert(m_iContextSizeWW >= m_iContextSizeCW);
	
	m_vMeanGlobal = NULL;
	m_vCovarianceGlobal = NULL;
}

// destructor
ContextModeling::~ContextModeling()
{
	if (m_vMeanGlobal) {
		delete m_vMeanGlobal;
	}
	if (m_vCovarianceGlobal) {
		delete m_vCovarianceGlobal;
	}
}

// tree-based n-phone clustering (top-down clustering using phonetic rules)
// the following information is needed:
// - estimation counts from all the observed units in the training data
// - a set of phonetic rules
bool ContextModeling::clusterContextDependentUnits(MAccumulatorLogical &mAccumulators, HMMManager *hmmManager) {

	int iOffset;

	// (1) get the logical accumulators
	int m_iAccumulators = -1;
	//MAccumulator &mAccumulators = m_hmmManager->getAccumulators();
	m_iAccumulators = (int)mAccumulators.size(); 
	
	// (2) create a list containing all the aplicable rules

	// load the phonetic rules
	PhoneticRulesManager *phoneticRulesManager = new PhoneticRulesManager(m_strFilePhoneticRules.c_str(),m_phoneSet);
	phoneticRulesManager->load();
	
	//phoneticRulesManager->print();
	VPhoneticRule *vPhoneticRule = phoneticRulesManager->getRules();
	if (vPhoneticRule->empty()) {
		return false;
	}
	int iPhoneticRules = (int)vPhoneticRule->size();
	
	// calculate the total number of rules
	int iRules = m_iContextSizeWW*2*iPhoneticRules + NUMBER_HMM_STATE_POSITIONS; // left/right context rules and within-word position rules
	if (m_bGlobalClustering) {
		iRules += iPhoneticRules + NUMBER_HMM_STATES;		// rules for basephone and HMM-state will be added
	}
	Rule **rules = new Rule*[iRules];
	// add the phonetic rules to the list
	for(int h = 0 ; h < m_iContextSizeWW ; ++h) {
		// left-contexts
		for(int j=0 ; j<iPhoneticRules ; ++j) {
			rules[j+iPhoneticRules*h] = new Rule;
			rules[j+iPhoneticRules*h]->iType = RULE_TYPE_PHONETIC_PHONE_LEFT;
			rules[j+iPhoneticRules*h]->iContextPosition = h;
			rules[j+iPhoneticRules*h]->phoneticRule = (*vPhoneticRule)[j];	
		}	
		// right-contexts
		for(int j=0 ; j<iPhoneticRules ; ++j) {
			rules[j+iPhoneticRules*(h+m_iContextSizeWW)] = new Rule;
			rules[j+iPhoneticRules*(h+m_iContextSizeWW)]->iType = RULE_TYPE_PHONETIC_PHONE_RIGHT;
			rules[j+iPhoneticRules*(h+m_iContextSizeWW)]->iContextPosition = h;
			rules[j+iPhoneticRules*(h+m_iContextSizeWW)]->phoneticRule = (*vPhoneticRule)[j];	
		}	
	}
	// within-word position
	iOffset = m_iContextSizeWW*2*iPhoneticRules;
	rules[iOffset] = new Rule;
	rules[iOffset]->iType = RULE_TYPE_POSITION;
	rules[iOffset]->iPosition = WITHIN_WORD_POSITION_START;
	rules[iOffset+1] = new Rule;
	rules[iOffset+1]->iType = RULE_TYPE_POSITION;
	rules[iOffset+1]->iPosition = WITHIN_WORD_POSITION_INTERNAL;
	rules[iOffset+2] = new Rule;
	rules[iOffset+2]->iType = RULE_TYPE_POSITION;
	rules[iOffset+2]->iPosition = WITHIN_WORD_POSITION_END;
	rules[iOffset+3] = new Rule;
	rules[iOffset+3]->iType = RULE_TYPE_POSITION;
	rules[iOffset+3]->iPosition = WITHIN_WORD_POSITION_MONOPHONE;
	
	if (m_bGlobalClustering) {
		// basephone
		iOffset = m_iContextSizeWW*2*iPhoneticRules+NUMBER_HMM_STATE_POSITIONS; 
		for(int i=0 ; i<iPhoneticRules ; ++i) {
			rules[iOffset+i] = new Rule;
			rules[iOffset+i]->iType = RULE_TYPE_PHONETIC_PHONE;
			rules[iOffset+i]->phoneticRule = (*vPhoneticRule)[i];	
		}	
		// HMM-state
		iOffset = (m_iContextSizeWW*2+1)*iPhoneticRules+NUMBER_HMM_STATE_POSITIONS; 
		for(int i=0 ; i<NUMBER_HMM_STATES ; ++i) {
			rules[iOffset+i] = new Rule;
			rules[iOffset+i]->iType = RULE_TYPE_STATE;
			rules[iOffset+i]->iState = i;
		}
	}
	
	// (2) do the actual clustering

	// local clustering: each HMM-state is clustered independently
	if (m_bGlobalClustering == false) {
	
		double dLikelihoodAccumulated = 0.0;
		double dLikelihoodAccumulatedClustering = 0.0;
		int iLeavesAccumulated = 0;
		int iLeavesAccumulatedData = 0;
		int iLeavesAccumulatedWithoutData = 0;
		
		int iContextDecisionTrees = m_phoneSet->size()*NUMBER_HMM_STATES;
		ContextDecisionTree **contextDecisionTrees = new ContextDecisionTree*[iContextDecisionTrees];
		
		// for each state
		for(unsigned int iBasePhone = 0 ; iBasePhone < m_phoneSet->size() ; ++iBasePhone) {
			for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			
				// (2.1) get a linked list with the basephone accumulators
				double dOccupation = 0.0;
				Accumulator *accumulator = NULL;
				for(MAccumulatorLogical::iterator it = mAccumulators.begin() ; it != mAccumulators.end() ; ++it) {	
					if ((it->second->getPhone() == iBasePhone) &&
						 (it->second->getState() == iState)) {
						dOccupation += it->second->getOccupation();
						it->second->setNext(accumulator);
						accumulator = it->second;
					}
				}
				
				// create the context decision tree
				ContextDecisionTree *contextDecisionTree = new ContextDecisionTree(m_iDim,
					m_iCovarianceModelling,m_phoneSet,m_iContextModelingOrderWW,accumulator,rules,iRules);
				contextDecisionTrees[iBasePhone*NUMBER_HMM_STATES+iState] = contextDecisionTree;
				
				// check whether there is data to cluster the state
				if (accumulator == NULL) {
					iLeavesAccumulated++;
					iLeavesAccumulatedWithoutData++;
					// log message
					BVC_WARNING << "no data available to estimate the parameters of the HMM-state: " << m_phoneSet->getStrPhone(iBasePhone) << "(" << iState <<"), defaulting to global distribution";
					continue;
				}
				
				// compute the likelihood of the root node (context independent phone)
				double dLikelihoodRoot = 0.0;
				double dOccupationRoot = 0.0;
				contextDecisionTree->computeLikelihoodRoot(&dLikelihoodRoot,&dOccupationRoot,*m_vCovarianceGlobal);
				dLikelihoodAccumulated += dLikelihoodRoot;
					
				// clustering: only phones that need context dependency
				double dLikelihoodAfterClustering = dLikelihoodRoot;
				int iLeaves = 1;
				int iLeavesData = 1;
				if (m_phoneSet->isPhoneContextModeled(iBasePhone)) {
	
					// to make sure the user knows what is being done
					if (((int)iBasePhone == m_phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE)) && (iState == 0)) {
						BVC_WARNING << "modeling context for phonetic symbol: " << PHONETIC_SYMBOL_SILENCE;
					}
					
					contextDecisionTree->clusterRoot(m_fMinimumClusterOccupation,m_fMinimumLikelihoodGainClustering);
					// get the likelihood from the leaves of the node
					dLikelihoodAfterClustering = contextDecisionTree->computeTreeLikelihood();
					dLikelihoodAccumulatedClustering += dLikelihoodAfterClustering;
					iLeaves = contextDecisionTree->countTreeLeaves();
					iLeavesData = contextDecisionTree->countTreeLeavesOccupancy();
					iLeavesAccumulated += iLeaves;
					iLeavesAccumulatedData += iLeavesData;	
				}
				// do not cluster: for example SIL, breath, noise, etc
				else {
				
					dLikelihoodAccumulatedClustering += dLikelihoodRoot;
					iLeavesAccumulated++;
					iLeavesAccumulatedData++;
				}
				double dPercentIncrease = 100.0*((dLikelihoodRoot-dLikelihoodAfterClustering)/dLikelihoodRoot);
				char strInformation[1024+1];
				sprintf(strInformation,"%s(%d) occ: %12.2f, likelihood: %16.4f -> %16.4f (%5.2f%%) %d %s",
					m_phoneSet->getStrPhone(iBasePhone),iState,dOccupationRoot,dLikelihoodRoot,
					dLikelihoodAfterClustering,dPercentIncrease,iLeavesData,
					Accumulator::getContextModelingOrder(m_iContextModelingOrderWW));
				BVC_VERB << strInformation;
			}
		}
		
		assert(iLeavesAccumulated == (iLeavesAccumulatedData+iLeavesAccumulatedWithoutData));
		
		double dLikelihoodDecreasePercentClustering = (dLikelihoodAccumulatedClustering/dLikelihoodAccumulated)*100;
		int iPhysicalNphones = 0;
		
		printf("likelihood before clustering:   %20.4f %6.2f%% (%8d observed %s, %d monophones)\n",
			dLikelihoodAccumulated,100.0,m_iAccumulators,Accumulator::getContextModelingOrder(m_iContextModelingOrderWW),
			m_iPhones*NUMBER_HMM_STATES);	
		printf("likelihood top-down clustering: %20.4f %6.2f%% (%8d physical %s)\n",
			dLikelihoodAccumulatedClustering,dLikelihoodDecreasePercentClustering,iLeavesAccumulated,
			Accumulator::getContextModelingOrder(m_iContextModelingOrderWW));
		
		// (3) bottom-up merging (optional)
		int iLeavesAfterMerging = 0;
		if (m_bBottomUpMerging) {
			
			double dLikelihoodAfterMerging = 0;
			for(int iHMM = 0 ; iHMM < m_iPhones*NUMBER_HMM_STATES ; ++iHMM) {
				// if there is no data associated to the HMM-state continue
				if (contextDecisionTrees[iHMM]->getOccupation() == 0.0) {
					iLeavesAfterMerging++;
					continue;
				}
				iLeavesAfterMerging += contextDecisionTrees[iHMM]->compactTreeLeaves(m_fMinimumLikelihoodGainClustering
					,*m_vCovarianceGlobal);
				dLikelihoodAfterMerging += contextDecisionTrees[iHMM]->computeTreeLikelihood();
			}	
		
			double dLikelihoodDecreasePercentClusteringMerging = (dLikelihoodAfterMerging/dLikelihoodAccumulated)*100;
		
			printf("likelihood bottom-up merging:   %20.4f %6.2f%% (%8d physical %s)\n",
				dLikelihoodAfterMerging,dLikelihoodDecreasePercentClusteringMerging,iLeavesAfterMerging,Accumulator::getContextModelingOrder(m_iContextModelingOrderWW));
			
			iPhysicalNphones = iLeavesAfterMerging;
		} else {
			iPhysicalNphones = iLeavesAccumulatedData+iLeavesAccumulatedWithoutData;
		}
		
		// create the new set of HMM-states from the leaves of the trees
		m_hmmStatesClustered = new HMMState*[iPhysicalNphones];
		for(int iState = 0 ; iState < iPhysicalNphones ; ++iState) {
			m_hmmStatesClustered[iState] = NULL; 
		}
		int iLeavesCheck = 0;
		double dLikelihoodCheck = 0.0;
		int iNphonesCreated = 0;
		int iOffset = 0;
		for(int iHMM = 0 ; iHMM < m_iPhones*NUMBER_HMM_STATES ; ++iHMM) {
			// if there is no data associated to the HMM-state continue
			if (contextDecisionTrees[iHMM]->getOccupation() == 0.0) {
				if ((m_vMeanGlobal == NULL) || (m_vCovarianceGlobal == NULL)) {
					computeGlobalDistribution(mAccumulators);
				}	
				unsigned char iPhoneAux = iHMM/NUMBER_HMM_STATES;
				unsigned char iStateAux = iHMM%NUMBER_HMM_STATES;
				contextDecisionTrees[iHMM]->createHMMStateNoData(iPhoneAux,iStateAux,
					m_hmmStatesClustered,&iLeavesCheck,	&iNphonesCreated,*m_vMeanGlobal,*m_vCovarianceGlobal);
				iOffset = iNphonesCreated;
				continue;
			}
			// create the nphones from the leaves of the tree
			contextDecisionTrees[iHMM]->createClusteredContextDependentUnitsFromLeaves(m_hmmStatesClustered,&iLeavesCheck,
				&iNphonesCreated,&dLikelihoodCheck);
			contextDecisionTrees[iHMM]->checkConsistency();	
			
			assert(iNphonesCreated > 0);
			iOffset = iNphonesCreated;
		}
		assert(iOffset == iLeavesAfterMerging);
		if (m_bBottomUpMerging) {
			assert(iLeavesAfterMerging == iNphonesCreated);
		} else {
			assert(iPhysicalNphones == iLeavesCheck);
		}
		for(int iState = 0 ; iState < iPhysicalNphones ; ++iState) {
			assert(m_hmmStatesClustered[iState] != NULL); 
		}
		
		// convert the HMMs to context-dependent
		if (hmmManager != NULL) {
			hmmManager->toContextDependentUnits(m_hmmStatesClustered,iPhysicalNphones,m_iContextModelingOrderWW,
				m_iContextModelingOrderCW,phoneticRulesManager,contextDecisionTrees,iContextDecisionTrees);
		}
		
		// write the trees to disk 
		/*if (m_strFileDecisionTrees != NULL) {
		
			if (storeContextDecisionTrees(m_strFileDecisionTrees,contextDecisionTrees,iContextDecisionTrees) == false) {
				BVC_ERROR << "unable to write the decision tree(s) to the file: " << m_strFileDecisionTrees;
			}
		}*/
	} 
	// global clustering: all the HMM-states are clustered together
	else {	
	
		float fLikelihood = 0.0;
		double dLikelihoodClustering = 0.0;
		int iLeaves = 0;
		int iLeavesData = 0;		
		
		// create the single root node
		DTNode *nodeRoot = new DTNode;
		nodeRoot->dtnodeYes = NULL;
		nodeRoot->dtnodeNo = NULL;
		nodeRoot->dtnodeMerged = NULL;
		nodeRoot->bIsMerged = false;
		nodeRoot->rule = NULL;
		nodeRoot->rules = new Rule*[iRules];
		nodeRoot->iRules = iRules;
		for(int i=0 ; i<iRules ; ++i) {
			nodeRoot->rules[i] = rules[i];
		}	
		nodeRoot->iHMMState = -1;
		
		// create a linked list with the basephone accumulators
		double dOccupation = 0.0;
		Accumulator *accumulator = NULL;
		for(MAccumulatorLogical::iterator it = mAccumulators.begin() ; it != mAccumulators.end() ; ++it) {
			dOccupation += it->second->getOccupation();
			it->second->setNext(accumulator);
			accumulator = it->second;
		}	
		
		// check whether there is enough data for clustering
		if (accumulator == NULL) {
			BVC_ERROR << "no data available to carry out context dependent clustering";
		}
		
		// create the context decision tree
		ContextDecisionTree *contextDecisionTree = new ContextDecisionTree(m_iDim,m_iCovarianceModelling,m_phoneSet,m_iContextModelingOrderWW,
			accumulator,rules,iRules);
				
		// compute the likelihood of the root node (context independent phone)
		double dLikelihoodRoot = 0.0;
		double dOccupationRoot = 0.0;
		contextDecisionTree->computeLikelihoodRoot(&dLikelihoodRoot,
			&dOccupationRoot,*m_vCovarianceGlobal);
		
		// do the clustering
		contextDecisionTree->clusterRoot(m_fMinimumClusterOccupation,m_fMinimumLikelihoodGainClustering);
		// get the likelihood from the leaves of the node
		dLikelihoodClustering += contextDecisionTree->computeTreeLikelihood();
		iLeaves = contextDecisionTree->countTreeLeaves();
		iLeavesData = contextDecisionTree->countTreeLeavesOccupancy();
	
		
		// compute the likelihood of the root node (context independent phone)
		/*bool bReturnValue = computeLikelihoodCluster(nodeRoot,NULL,true,&(nodeRoot->fLikelihood),NULL);
		assert(bReturnValue == true);
		fLikelihood = nodeRoot->fLikelihood;
		// cluster the node
		clusterNode(nodeRoot);
		// get the likelihood from the leaves of the node
		fLikelihoodClustering = computeTreeLikelihood(nodeRoot);
		iLeaves = countTreeLeaves(nodeRoot);
		iLeavesData = countTreeLeavesOccupancy(nodeRoot);*/
					
		if (iLeaves != iLeavesData) {
			printf("warning: leaves: %d leaves with data: %d\n",iLeaves,iLeavesData);
		}
		
		//printf("(top-down triphone clustering) likelihood decrease: %.4f %.4f (logical/physical triphones: %d %d)\n",fLikelihoodAccumulated,fLikelihoodAccumulatedClustering,m_iAccumulatorsAllocated,iLeavesAccumulated);
		
		double dLikelihoodDecreasePercentClustering = (dLikelihoodClustering/fLikelihood)*100;
		
		printf("likelihood before clustering:   %20.4f %6.2f%% (%8d logical triphones, %d monophones)\n",fLikelihood,100.0,m_iAccumulators,m_iPhones*NUMBER_HMM_STATES);	
		printf("likelihood top-down clustering: %20.4f %6.2f%% (%8d physical triphones)\n",dLikelihoodClustering,dLikelihoodDecreasePercentClustering,iLeaves);
		
		// (3) bottom-up merging (optional)		
		/*int iLeavesAfterMerging = 0;
		if (m_bBottomUpMerging == true) {
		
			iLeavesAfterMerging = compactTreeLeaves(nodeRoot);
			float fLikelihoodAfterMerging = computeTreeLikelihood(nodeRoot);
			float fLikelihoodDecreasePercentClusteringMerging = (fLikelihoodAfterMerging/fLikelihood)*100;
		
			printf("likelihood bottom-up merging:   %20.4f %6.2f%% (%8d physical triphones)\n",fLikelihoodAfterMerging,fLikelihoodDecreasePercentClusteringMerging,iLeavesAfterMerging);
			
			iPhysicalNphones = iLeavesAfterMerging;
		} else {
			iPhysicalNphones = iLeavesData;
		}
		
		// create the new set of HMM-states from the leaves of the tree
		//m_hmmManagerTriphones = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
		m_hmmStatesClustered = new HMMState*[iPhysicalNphones];
		for(int iState = 0 ; iState < iPhysicalNphones ; ++iState) {
			m_hmmStatesClustered[iState] = NULL; 
		}
		int iLeavesCheck = 0;
		float fLikelihoodCheck = 0.0;
		int iNphonesCreated = 0;
		createClusteredTriphonesFromLeaves(nodeRoot,&iLeavesCheck,&iNphonesCreated,&fLikelihoodCheck);
		
		// asserts
		if (m_bBottomUpMerging == true) {
			assert(iNphonesCreated == iLeavesAfterMerging);
		} else {
			assert(iPhysicalNphones == iLeavesCheck);
		}
		for(int iState = 0 ; iState < iPhysicalNphones ; ++iState) {
			assert(m_hmmStatesClustered[iState] != NULL); 
		}	
		
		// create the triphone table
		HMMState **hmmStateTriphoneTable = createTriphoneTableGlobal(nodeRoot);
		
		int iHMMStatesLogical = m_iPhones*m_iPhones*m_iPhones*NUMBER_HMM_STATE_POSITIONS*NUMBER_HMM_STATES;
		m_hmmManager->toContextDependentTriphones(m_hmmStatesClustered,iPhysicalNphones,HMM_CONTEXT_MODELING_TRIPHONES,hmmStateTriphoneTable,iHMMStatesLogical);
		
		// clean-up	
		destroyTree(nodeRoot);*/
	}		
	
	// clean-up
	/*for(int i=0 ; i<iRules ; ++i) {
		delete rules[i];
	}
	delete [] rules;
	delete phoneticRulesManager;*/
	
	// clean-up
	delete [] rules;
	
	return true;
}

// store context decision trees to a file
void ContextModeling::storeContextDecisionTrees(const char *strFile, ContextDecisionTree **contextDecisionTrees, int iContextDecisionTrees) {

	FileOutput file(strFile,true);
	file.open();
	
	// write the clustering method
	unsigned char iMethod;
	m_bGlobalClustering ? iMethod = CLUSTERING_METHOD_GLOBAL : iMethod = CLUSTERING_METHOD_LOCAL;
	IOBase::write(file.getStream(),m_iContextModelingOrderWW);	
	IOBase::write(file.getStream(),m_iPhones);	
	IOBase::write(file.getStream(),iContextDecisionTrees);	
	
	// write the decision trees one by one
	for(int i=0 ; i < iContextDecisionTrees ; ++i) {
		contextDecisionTrees[i]->store(file);
	}
	
	file.close();
}

// load context decision trees from a file
ContextDecisionTree **ContextModeling::loadContextDecisionTrees(const char *strFile, int iFeatureDimensionality, int iCovarianceModelling, PhoneticRulesManager *phoneticRulesManager, int *iContextDecisionTrees, PhoneSet *phoneSet, unsigned char iContextModelingOrder, unsigned char iClusteringMethod) {

	FileInput file(strFile,true);
	file.open();	

	int iClusteringMethodAux = -1;
	IOBase::read(file.getStream(),&iClusteringMethodAux);
	assert(iClusteringMethodAux == iClusteringMethod);
	
	// read the context modeling order
	unsigned char iContextModelingOrderAux;
	IOBase::read(file.getStream(),&iContextModelingOrderAux);
	assert(iContextModelingOrderAux == iContextModelingOrder);
	
	// read the number of phones
	unsigned char iPhonesAux;
	IOBase::read(file.getStream(),&iPhonesAux);
	assert(iPhonesAux == phoneSet->size());
	
	// read the number of trees
	IOBase::read(file.getStream(),iContextDecisionTrees);
	assert(*iContextDecisionTrees > 0);
	
	// allocate memory for the array of trees
	ContextDecisionTree **contextDecisionTrees = new ContextDecisionTree*[*iContextDecisionTrees];

	// read the decision trees one by one
	for(int i=0 ; i < *iContextDecisionTrees ; ++i) {
		contextDecisionTrees[i] = ContextDecisionTree::load(file,iFeatureDimensionality,iCovarianceModelling,
			phoneSet,phoneticRulesManager,iContextModelingOrder);
		assert(contextDecisionTrees[i]);
	}
	
	file.close();

	return contextDecisionTrees;
}

// compute the global distribution given the accumulators
void ContextModeling::computeGlobalDistribution(MAccumulatorLogical mAccumulators) {

	HMMManager::getCovarianceElements(m_iDim,m_iCovarianceModelling);

	// compute the global distribution
	Vector<double> vObservation(m_iDim);	
	Vector<double> vObservationSquare(m_iDim);	
	m_vMeanGlobal = new Vector<float>(m_iDim);
	m_vCovarianceGlobal = new Vector<float>(m_iDim);	
	double dOccupationTotal = 0.0;
	
	// compute the global mean
	for(MAccumulatorLogical::iterator it = mAccumulators.begin() ; it != mAccumulators.end() ; ++it) {	
		Accumulator *accumulator = it->second;
		vObservation.add(accumulator->getObservation());
		dOccupationTotal += accumulator->getOccupation();
	}
	m_vMeanGlobal->mul((float)(1.0/dOccupationTotal),vObservation);	 
	
	// compute global covariance	
	for(MAccumulatorLogical::iterator it = mAccumulators.begin() ; it != mAccumulators.end() ; ++it) {	
		Accumulator *accumulator = it->second;
		vObservationSquare.add(accumulator->getObservationSquareDiag());
	}	
		
	// covariance is diagonal
	m_vCovarianceGlobal->mul((float)(1.0/dOccupationTotal),vObservationSquare);
	m_vCovarianceGlobal->addSquare(-1.0,*m_vMeanGlobal);
}

};	// end-of-namespace
