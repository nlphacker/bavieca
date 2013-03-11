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


#ifndef CONTEXTDECISIONTREE_H
#define CONTEXTDECISIONTREE_H

using namespace std;

#include <list>

#include "Accumulator.h"
#include "Global.h"
#include "HMMState.h"
#include "PhoneSet.h"
#include "PhoneticRulesManager.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "Vector.h"

namespace Bavieca {

// types of rule used for clustering
#define RULE_TYPE_PHONETIC_PHONE_LEFT			0
#define RULE_TYPE_PHONETIC_PHONE					1				// only for "global" clustering
#define RULE_TYPE_PHONETIC_PHONE_RIGHT			2
#define RULE_TYPE_POSITION							3
#define RULE_TYPE_STATE								4				// only for "global" clustering

// rule used for clustering context dependent units
typedef struct _Rule {
	unsigned char iType;
	unsigned char iContextPosition;
	union {
		PhoneticRule *phoneticRule;
		unsigned char iState;	
		unsigned char iPosition;
	};
} Rule;

// node in the decision tree
typedef struct _DTNode {
	Rule **rules;							// array of rules (phonetic, state, etc) to be applied (the last position is NULL)
	int iRules;								// number of rules that are left to be applied
	double dOccupation;					// number of feature vectors aligned to nphones in this node
	double dLikelihood;					// likelihood
	Accumulator *accumulator;			// list of accumulators (each accumulator corresponds to a logical state)
	Rule *rule;								// rule applied
	_DTNode *dtnodeYes;					// yes child (corresponds to a positive answer to the phonetic question)
	_DTNode *dtnodeNo;					//	no child (corresponds to a negative answer to the phonetic question)
	_DTNode *dtnodeMerged;				// leaf with which the node was merged during bottom-up merging (or NULL if it wasn't)
	bool bIsMerged;						// whether the leaf was merged to another leaf or another leaf was merged with the leaf
	int iHMMState;							// HMM-state index
} DTNode;

typedef list<DTNode*> LDTNode;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ContextDecisionTree {

	private:
	
		int m_iDim;												// feature dimensionality
		int m_iCovarianceModelling;						// covariance modelling type
		PhoneSet *m_phoneSet;								// phonetic symbol set
		unsigned char m_iContextModelingOrder;			// context modeling order
		DTNode *m_dtnodeRoot;								// tree root
		DTNode *m_dtnodes;									// tree nodes
		int m_iNodes;											// # tree nodes
		int m_iLeaves;											// # tree leaves
		double m_dOccupation;								// tree occupation
		
		// estimation statistics
		bool m_bEstimationDataAllocated;		 
		
		Vector<double> *m_vMeanCluster;
		Vector<double> *m_vCovarianceCluster;
		double m_dOccupationCluster;
		Vector<double> *m_vObservationCluster;
		Vector<double> *m_vObservationSquareCluster;
		
		Vector<double> *m_vMeanClusterNo;
		Vector<double> *m_vCovarianceClusterNo;
		double m_dOccupationClusterNo;
		Vector<double> *m_vObservationClusterNo;
		Vector<double> *m_vObservationSquareClusterNo;			
		
	private:
	
		// compute the likelihood of a decission tree by summing up the likelihood in the leaves
		double computeTreeLikelihood(DTNode *node);
		
		// return the number of leaves in the tree (including those leaves merged to another leaf)
		int countTreeLeaves(DTNode *node);
		
		// return the number of leaves in the tree (excluding those leaves merged to another leaf)
		int countTreeLeavesUnique(DTNode *node);
		
		// compute the likelihood of a decission tree by summing up the likelihood in the leaves
		int countTreeLeavesOccupancy(DTNode *node);			
		
		// create the new set of HMM-states from the leaves of the tree (return the number of physical nphones created)
		void createClusteredContextDependentUnitsFromLeaves(DTNode *dtnode, HMMState **hmmStates, int *iLeaves, int *iUnitsCreated, double *dLikelihood);
		
		// destroy a decission tree
		void destroyTree(DTNode *node);
		
		// check tree consistency
		bool checkConsistency(DTNode *node);
		
	public:

		// constructor
		ContextDecisionTree(int iFeatureDimensionality, int iCovarianceModelling, PhoneSet *phoneSet, unsigned char iContextModelingOrder, Accumulator *accumulator, Rule **rules, int iRules);

		// constructor
		ContextDecisionTree(int iFeatureDimensionality, int iCovarianceModelling, PhoneSet *phoneSet, unsigned char iContextModelingOrder);

		// constructor
		ContextDecisionTree();
		
		// copy constructor (creates a replica of the object, for tree traversal only)
		ContextDecisionTree(ContextDecisionTree *contextDecisionTreeSeed);

		// destructor
		~ContextDecisionTree();
		
		// load from file
		static ContextDecisionTree *load(FileInput &file, int iFeatureDimensionality, int iCovarianceModelling, 
			PhoneSet *phoneSet, PhoneticRulesManager *phoneticRulesManager, unsigned char iContextModelingOrder);
		
		// store to file
		void store(FileOutput &file);
		
		// cluster a node (if possible)
		void clusterNode(DTNode *node, float fMinimumClusterOccupation, float fMinimumLikelihoodGainClustering);	
		
		// cluster a node (if possible)
		inline void clusterRoot(float fMinimumClusterOccupation, float fMinimumLikelihoodGainClustering) {
		
			return clusterNode(m_dtnodeRoot,fMinimumClusterOccupation,fMinimumLikelihoodGainClustering);
		}
		
		// compute the likelihood of a cluster of context dependent units 
		void computeLikelihoodCluster(DTNode *node, double *dLikelihoodCluster, double *dOccupationCluster, Vector<float> &vCovarianceGlobal);		
		
		// compute the likelihood of the root node
		inline void computeLikelihoodRoot(double *dLikelihoodCluster, double *dOccupationCluster, Vector<float> &vCovarianceGlobal) {
		
			computeLikelihoodCluster(m_dtnodeRoot,dLikelihoodCluster,dOccupationCluster,vCovarianceGlobal);
			m_dtnodeRoot->dLikelihood = *dLikelihoodCluster;
		}
		
		// compute the likelihood of a cluster of nphones given the parent cluster, the rule and the decission made
		// note: return false if after applying the rule:
		// - there is not enough data to robustly estimate the parameters (minimum occupation count)
		// - there is too much data resulting from the answer to the question, so the opposite answer will produce not enough data
		bool computeLikelihoodRule(DTNode *node, Rule *rule, bool bAnswer, double *dLikelihoodCluster, double *dOccupationCluster, float fMinimumClusterOccupation);		
		
		// compute the likelihood of the clusters resulting from applying a question to the given cluster
		// note: return false if after applying the rule:
		// - there is not enough data to robustly estimate the parameters (minimum occupation count)
		// - there is too much data resulting from the answer to the question, so the opposite answer wont produce enough data
		bool computeLikelihoodRule(DTNode *node, Rule *rule, double *dLikelihoodClusterYes, double *dLikelihoodClusterNo, double *dOccupationClusterYes, double *dOccupationClusterNo, float fMinimumClusterOccupation);
		
		// compute the likelihood of a decission tree by summing up the likelihood in the leaves
		inline double computeTreeLikelihood() {
		
			return computeTreeLikelihood(m_dtnodeRoot);
		}
		
		// return the number of leaves in the tree (including those leaves merged to another leaf)
		inline int countTreeLeaves() {
		
			return countTreeLeaves(m_dtnodeRoot);
		}
		
		// compute the likelihood of a decission tree by summing up the likelihood in the leaves
		inline int countTreeLeavesOccupancy() {
		
			return countTreeLeaves(m_dtnodeRoot);
		}
		
		// return the number of leaves in the tree (excluding those leaves merged to another leaf)
		inline int countTreeLeavesUnique() {
		
			return countTreeLeavesUnique(m_dtnodeRoot);
		}
		
		// get the tree leaves
		void getTreeLeaves(DTNode *node, LDTNode &lLeaves);
		
		// print tree leaves
		void printTreeLeaves(DTNode *node);
		
		// check tree consistency
		inline bool checkConsistency() {
			
			return checkConsistency(m_dtnodeRoot);
		}
		
		// compute the likelihood of the cluster resulting from mergining two leaves
		double computeLikelihoodFromMerging(DTNode *nodeA, DTNode *nodeB, Vector<float> &vCovarianceGlobal);	
		
		// compact the tree leaves by merging those leaves that when merged the likelihood decrease is below the minimum splitting gain used for clustering (return the resulting number of leaves)
		int compactTreeLeaves(float fMinimumLikelihoodGainClustering, Vector<float> &vCovarianceGlobal);
		
		// create the new set of HMM-states from the leaves of the tree (return the number of physical nphones created)
		inline void createClusteredContextDependentUnitsFromLeaves(HMMState **hmmStates, int *iLeavesSeen, int *iUnitsCreated, double *dLikelihood) {
		
			return createClusteredContextDependentUnitsFromLeaves(m_dtnodeRoot,hmmStates,iLeavesSeen,iUnitsCreated,dLikelihood);
		}
		
		// create an HMM-state when there is no data associated to it (default to globla distribution)
		void createHMMStateNoData(unsigned char iPhone, unsigned char iState, HMMState **hmmStates, int *iLeavesSeen, int *iUnitsCreated, Vector<float> &vMeanGlobal, Vector<float> &vCovarianceGlobal);
		
		// return the tree occupation
		inline double getOccupation() {
		
			return m_dOccupation;
		}
		
		// compute the cluster likelihood
		double computeLikelihoodCluster(double dOccupation, Vector<double> &vCovariance);
		
		// compute the cluster likelihood
		double computeLikelihoodCluster(double dOccupation, Vector<float> &vCovariance);
		
		// return the HMM-state associated to the given context-dependent phone
		int getHMMIndex(unsigned char *iPhoneLeft, unsigned char iPhone, 
			unsigned char *iPhoneRight, unsigned char iPosition, 
			unsigned char iState);
			
		bool question(Rule *rule, Accumulator *accumulator) {
		
			switch(rule->iType) {
			
				case RULE_TYPE_PHONETIC_PHONE_LEFT: {
					unsigned char iPhone = accumulator->getLeftPhone(rule->iContextPosition);
					return rule->phoneticRule->bPhone[iPhone];
				}
				case RULE_TYPE_PHONETIC_PHONE: {
					unsigned char iPhone = accumulator->getPhone();
					return rule->phoneticRule->bPhone[iPhone];
				}
				case RULE_TYPE_PHONETIC_PHONE_RIGHT: {
					unsigned char iPhone = accumulator->getRightPhone(rule->iContextPosition);
					return rule->phoneticRule->bPhone[iPhone];
				}
				case RULE_TYPE_POSITION: {	
					unsigned char iPosition = accumulator->getPosition();
					return (rule->iPosition == iPosition);
				}
				case RULE_TYPE_STATE: {
					unsigned char iState = accumulator->getState();
					return (rule->iState == iState);
				}
				default: {
					assert(0);	
					break;
				}
			}
			
			return false;
		}
};

};	// end-of-namespace

#endif
