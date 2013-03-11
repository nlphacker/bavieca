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


#ifndef REGRESSIONTREEFILE_H
#define REGRESSIONTREEFILE_H

using namespace std;

#include <vector>
#include <queue>
#include <list>

#include "HMMManager.h"

#include "Matrix.h"
#include "MatrixStatic.h"
#include "Vector.h"
#include "VectorStatic.h"

namespace Bavieca {

// clustering method
#define CLUSTERING_METHOD_KMEANS				0
#define CLUSTERING_METHOD_EM					1

// transform type 
#define TRANSFORM_MLLR_MEAN						0
#define TRANSFORM_MLLR_MEAN_COVARIANCE			1

// general clustering parameters
#define MINIMUM_CLUSTER_SIZE_DEFAULT			50		// minimum # Gaussian comp per cluster (to prevent singularities)

// EM clustering parameters
#define EM_CONVERGENCE_EPSILON      			0.001			// convergence for the EM-algorithm
#define EM_MAXIMUM_ITERATIONS					100				// maximum number of iterations
#define EM_MINIMUM_LIKELIHOOD		   			FLT_MIN/2.0		// minimum likelihood value (to prevent singularities)
#define EM_MAXIMUM_LIKELIHOOD		   			FLT_MAX/2.0		// minimum likelihood value (to prevent singularities)

// node type according to its position in the tree
#define NODE_POSITION_INTERNAL					0
#define NODE_POSITION_LEAF						1

typedef struct {
	GaussianDecoding *gaussian;
	double dOccupation;
	Vector<double> *vObservation;
} GaussianStats;

typedef vector<GaussianStats*> VGaussianStats;

typedef struct _MLLRTransform MLLRTransform;

typedef struct _RTNode {
	int iId;									// node identifier
	int iBaseClass;								// base class
	float fOccupation;							// Gaussian occupation of the subtree whose root is this node
	_RTNode *left;								// left child
	_RTNode *right;								// right child
	// fields needed to create the tree (data-driven Gaussian clustering)
	float fDistortion;							// distortion of elements in the tree-node
	int iDepth;									// tree-node depth
	bool bTouched;								// general purpose flag
	VGaussianDecoding	vGaussian;				// gaussians at this node (keeps gaussians clustered together or gaussians with adaptation data)
	VGaussianStats vGaussianStats;
	Matrix<double> *matrixObservationTranspose;	// stats at the base-class level that are needed to compute the covariance transform
	MLLRTransform *transform;					// MLLR tranform associated to this node
} RTNode;

typedef queue<RTNode*> QRTNode;
typedef list<RTNode*> LRTNode;

typedef struct _MLLRTransform {
	int iId;								// transform id
	RTNode *node;							// tree node 
	double dOccupation;						// number of frames used to compute the transform
	Matrix<float> *matrixK; 				// stats to compute mean transform
	Matrix<float> **matrixG; 				// stats to compute mean transform
	Matrix<float> *matrixMean; 				// mean transform
	Matrix<double> *matrixCovariance; 		// covariance transform	
} MLLRTransform;

typedef vector<MLLRTransform*> VMLLRTransform;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class RegressionTree {

	private:
	
		HMMManager *m_hmmManager;
		int m_iDim;
		int m_iHMMStates;
		int m_iGaussians;
		HMMStateDecoding *m_hmmStates;
		RTNode **m_baseClasses;
		int m_iBaseClasses;
		RTNode *m_nodeRoot;					// tree-node	
		int m_iNodes;						// total number of nodes in the tree
		int *m_iGaussianBaseClass;
		
		// transform related fields
		MLLRTransform *m_transformGlobal;	// global transform
		MLLRTransform **m_transformMap;		// maps base-classes to transforms
		VMLLRTransform m_vMLLRTransform;	// MLLR transforms computed from the adaptation data
		
		// auxiliar variables 
		int m_iId;							// used to build the regression tree
		float *m_fCholeskyFactor;			// cholesky factor
		
		// return a random number from a range
		int getRandomNumber(int iStart, int iEnd) {
			
			assert(iStart<=iEnd);
		
			#if defined __linux__ || defined __APPLE__	
			// on linux 32 bits RAND_MAX = 2147483647 (this is the maximum signed integer)
			assert(RAND_MAX > iEnd-iStart);			
			// this is a very weak random number
			int iRandom = iStart + rand()%(iEnd-iStart+1);	
			assert((iRandom >= iStart) && (iRandom <= iEnd));
			return iRandom;
			#endif
		
			#ifdef _WIN32
			assert(0);
			#endif
			
			/*assert((double)(RAND_MAX*(RAND_MAX+2)) >= (iEnd-iStart));
		
			float fRandom;
			do {
				fRandom = rand()*(RAND_MAX+1)+rand();					// random number in [0,(RAND_MAX*(RAND_MAX+2))]
			} while(fRandom == (RAND_MAX*(RAND_MAX+2)));
			fRandom /= (RAND_MAX*(RAND_MAX+2));							// random number in [0,1)
			fRandom *= iEnd-iStart+1;										
			fRandom += iStart;												// random number in [iStart,iEnd+1)
			int iRandom = floorf(fRandom);								// integer random number in [iStart,iEnd]]
			assert((iRandom >= iStart) && (iRandom <= iEnd));
			
			return iRandom; */
			return 0;
		}
		
		// return the number of Gaussian components observed by leaf nodes hanging from teh given node
		int getGaussianComponentsObserved(RTNode *node) {
		
			if ((node->left == NULL) && (node->right == NULL)) {
				return (int)node->vGaussianStats.size();
			} else {
				assert((node->left != NULL) && (node->right != NULL));
				return getGaussianComponentsObserved(node->left)+getGaussianComponentsObserved(node->right);
			}	
		}
		
		// return the total occupation of the leaf nodes hanging from the given node
		double getOccupation(RTNode *node) {
		
			if ((node->left == NULL) && (node->right == NULL)) {
				double dOccupation = 0.0;
				for(VGaussianStats::iterator it = node->vGaussianStats.begin() ; it != node->vGaussianStats.end() ; ++it) {
					dOccupation += (*it)->dOccupation;
				}	
				return dOccupation;
			} else {
				assert((node->left != NULL) && (node->right != NULL));
				return getOccupation(node->left)+getOccupation(node->right);
			}	
		}		
		
		// comparison function
		static bool compareSize(const RTNode *nodeA, const RTNode *nodeB) {
			
			return (nodeA->vGaussian.size() >= nodeB->vGaussian.size());
		}		
		
		// detroy the regression tree
		void destroy(RTNode *node);
		
		// remove centroids
		void removeCentroids(RTNode *node);	
		
		// counts the number of leaves in the tree (recursively)
		int countLeaves(RTNode *rtNode);	
		
		// perform k-means to split a tree-node
		// note: it returns the new tree-nodes as children of the given tree-node
		bool kMeansClustering(RTNode *rtNode, int iMinimumComponentsCluster);
		
		// compute P(x|A)	
		float computeClusterLikelihood(float *fMean, float *fCovariance, float *fX);
		
		// perform EM clustering to split a tree-node
		// note: it returns the new tree-nodes as children of the given tree-node
		bool expectationMaximizationClustering(RTNode *cluster, int iMinimumComponentsCluster);	
		
		// compute the weight vector for the computation of the weighted euclidean distance
		float *computeWeight(RTNode *rtNode);
		
		// compute the euclidean distance between two vectors
		// note: weight is optional
		float euclideanDistance(float *fV1, float *fV2, float *fW = NULL) {
		
			float fAcc = 0.0;
			float fAux;
			for(int i=0 ; i<m_iDim ; ++i) {
				fAux = fV1[i]-fV2[i];
				if (fW) {
					fAux *= fW[i];
				}
				fAcc += fAux*fAux;
			}
			fAcc /= ((float)m_iDim);
			
			return fAcc;
		}	
		
		// print the regression tree
		void print();
		
		// compute a transform from the adaptation data (the given vector of gaussians)
		void computeTransform(MLLRTransform *transform, VGaussianStats &vGaussianStats, bool bMeanOnly);		
		
		MLLRTransform *newMLLRTransform(RTNode *node, double dOccupation) {
		
			MLLRTransform *transform = new MLLRTransform;	
			transform->node = node;
			transform->dOccupation = dOccupation;
			transform->matrixK = NULL;
			transform->matrixG = NULL;
			transform->matrixMean = NULL;
			transform->matrixCovariance = NULL;
			
			return transform;
		}
		
		// identify the tranforms that will be computed from the tree (recursively)
		float identifyTransformsToCompute(float fMinimumOccupationTransform, int iMinimumGaussianComponentsObserved, 
			RTNode *node);
		
		// apply a transform to a Gaussian
		void applyTransform(MLLRTransform *transform, GaussianDecoding *gaussian, bool bMean, bool bCovariance) {
		
			// mean transformation
			if (bMean) {
				assert(transform->matrixMean);
		
				// create the extended mean vector
				VectorStatic<float> vMean(gaussian->fMean,m_iDim);
				Vector<float> vMeanEx(vMean);	
				vMeanEx.appendFront(1.0);
				vMean.mul(*transform->matrixMean,vMeanEx);
			}
			
			// covariance transformation
			if (bCovariance) {
				assert(transform->matrixCovariance);
				
				VectorStatic<float> vCovariance(gaussian->fCovariance,m_iDim);
				Vector<double> vCInv(vCovariance);
			#ifdef OPTIMIZED_COMPUTATION
				vCInv.mul(2.0);
				vCInv.invertElements();
			#endif
				vCInv.invertElements();
				vCInv.sqrt();
				vCInv.invertElements();
				
				// apply the transform
				vCovariance.copy(vCInv); 
				vCovariance.mulDiagonal(*transform->matrixCovariance); 
				vCovariance.mulElements(vCInv); 
			}
		}
					
		// count the number of nodes in the tree
		int countNodes(RTNode *node) {
		
			if (node->left == NULL) {
				assert(node->right == NULL);
				return 1;
			} else {
				return countNodes(node->left)+countNodes(node->right)+1;	
			}	
		}
		

	public:

		// constructor
		RegressionTree(HMMManager *hmmManager);

		// destructor
		~RegressionTree();	
		
		// compute MLLR adaptation transforms
		void computeTransforms(float fMinimumOccupationTransform, int iMinimumGaussianComponentsObserved,
			GaussianStats **gaussianStats, bool bMeanOnly);	
			
		// apply the available transforms to the set of HMM-models
		void applyTransforms() {
			
			return applyTransforms(true,isCovarianceTransform());
		}
		
		// apply the transforms to the set of HMM-models
		void applyTransforms(bool bMean, bool bCovariance);
		
		// print the transforms
		void printTransforms();	
		
		// store the transforms 
		void storeTransforms(const char *strFile);
		
		// load the transforms from a file
		void loadTransforms(const char *strFile);
		
		// build the regression tree by clustering the gaussian components and 
		// assigns a baseclass to every gaussian component
		void build(int iBaseClasses, unsigned char iMethod, unsigned int iMinimumComponentsCluster);
		
		// store the tree into disk
		void store(const char *strFile);
		
		// load the tree from disk
		void load(const char *strFile);
		
		// return the number of base classes
		int getBaseClasses() {
		
			return m_iBaseClasses;
		}		
		
		// accumulate statistics into a base-class
		void accumulateStatistics(float *fObservation, double dOccupation, GaussianStats *gaussianStats) {
			
			assert(m_baseClasses[m_iGaussianBaseClass[gaussianStats->gaussian->iId]]->matrixObservationTranspose);
			
			// compute the cholesky factor of the inverse covariance matrix
			for(int i=0 ; i<m_iDim ; ++i) {
				float fCovariance = gaussianStats->gaussian->fCovariance[i];
			#ifdef OPTIMIZED_COMPUTATION
				fCovariance = 1.0f/(2.0f*fCovariance);
			#endif
				m_fCholeskyFactor[i] = sqrt(1.0f/fCovariance);
			}	
			
			// feature vector transpose
			int iBaseClass = m_iGaussianBaseClass[gaussianStats->gaussian->iId];
			Matrix<double> *matrixObservationTranspose = m_baseClasses[iBaseClass]->matrixObservationTranspose;
			for(int i=0 ; i<m_iDim ; ++i) {
				for(int j=0 ; j<m_iDim ; ++j) {	
					(*matrixObservationTranspose)(i,j) += m_fCholeskyFactor[i]*dOccupation*fObservation[j]*fObservation[i]*m_fCholeskyFactor[j];
				}
			}
		}
		
		// return the transform associated to the given Gaussian component
		MLLRTransform *getTransform(GaussianStats *gaussianStats) {
		
			int iBaseClass = m_iGaussianBaseClass[gaussianStats->gaussian->iId];
			return m_transformMap[iBaseClass];
		}
		
		// return the transforms
		VMLLRTransform *getTransforms() {
				
			return &m_vMLLRTransform;
		}
		
		// return whether there is a covariance transform
		bool isCovarianceTransform() {
		
			assert(m_vMLLRTransform.empty() == false);
		
			return (m_vMLLRTransform[0]->matrixCovariance != NULL);
		}

};

};	// end-of-namespace

#endif
