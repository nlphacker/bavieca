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


#include "RegressionTree.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "TimeUtils.h"

#include <iomanip>

#ifdef _WIN32
#include <time.h>
#endif

namespace Bavieca {

// contructor
RegressionTree::RegressionTree(HMMManager *hmmManager)
{
	m_hmmManager = hmmManager;
	m_iDim = hmmManager->getFeatureDimensionality();	
	m_hmmStates = m_hmmManager->getHMMStatesDecoding(&m_iHMMStates);	
	m_nodeRoot = NULL;
	m_baseClasses = NULL;
	m_transformMap = NULL;
	m_iGaussians = -1;
	
	m_fCholeskyFactor = new float[m_iDim];
	m_iGaussianBaseClass = new int [m_hmmManager->getNumberGaussianComponents()];
}

// destructor
RegressionTree::~RegressionTree()
{
	if (m_vMLLRTransform.empty() == false) {
		for(VMLLRTransform::iterator it = m_vMLLRTransform.begin() ; it != m_vMLLRTransform.end() ; ++it) {
			if ((*it)->matrixK) {
				delete (*it)->matrixK;
			}
			if ((*it)->matrixG) {
				for(int i=0 ; i < m_iDim ; ++i) {
					delete (*it)->matrixG[i];
				}
				delete [] (*it)->matrixG;
			}
			if ((*it)->matrixMean) {
				delete (*it)->matrixMean;
			}
			if ((*it)->matrixCovariance) {
				delete (*it)->matrixCovariance;
			}
			delete *it;
		}
	}

	if (m_transformMap != NULL) {
		delete [] m_transformMap;
	}	
	if (m_baseClasses != NULL) {
		delete [] m_baseClasses;
	}
	
	// destroy the actual regression tree
	destroy(m_nodeRoot);
	m_nodeRoot = NULL;
	
	if (m_fCholeskyFactor != NULL) {
		delete [] m_fCholeskyFactor;
	}
}


// build the regression tree by clustering the gaussian components and 
// assigns a baseclass to every gaussian component
void RegressionTree::build(int iBaseClasses, unsigned char iClusteringMethod, unsigned int iMinimumComponentsCluster) {

	m_iBaseClasses = iBaseClasses;
	float fDistortionTotal = 0.0;
	LRTNode lLeafNode;
	m_iId = 0;
	m_iNodes = 0;
	
	assert((iClusteringMethod == CLUSTERING_METHOD_KMEANS) || (iClusteringMethod == CLUSTERING_METHOD_EM));
	
	// starting time
	double dTimeStart = TimeUtils::getTimeMilliseconds();
	
	// initialize the random number generator
	srand((unsigned int)time(NULL));

	// (1) create a sinle root cluster containing all the gaussians
	m_nodeRoot = new RTNode;	
	++m_iNodes;	
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		int iGaussians = -1;
		GaussianDecoding *gaussians = m_hmmStates[i].getGaussians(iGaussians);
		for(int j=0 ; j < iGaussians ; ++j) {
			m_nodeRoot->vGaussian.push_back(&gaussians[j]);
			m_iGaussianBaseClass[gaussians[j].iId] = -1;
		}
	}	
	m_iGaussians = m_nodeRoot->vGaussian.size();
	
	// important: the number of gaussians is an upper bound to the number of base-classes
	if (m_iGaussians < m_iBaseClasses) {
		printf("warning: # base-classes can't be bigger than the # of gaussians, reduced to %d\n",m_iGaussians);
		m_iBaseClasses = m_iGaussians; 
	}
	
	m_nodeRoot->iId = m_iId++;
	m_nodeRoot->iBaseClass = -1;
	m_nodeRoot->fOccupation = 0.0;
	m_nodeRoot->fDistortion = 0.0;
	m_nodeRoot->iDepth = 0;
	m_nodeRoot->bTouched = true;
	m_nodeRoot->left = NULL;
	m_nodeRoot->right = NULL;
	m_nodeRoot->matrixObservationTranspose = NULL;	
	m_nodeRoot->transform = NULL;
	assert(lLeafNode.empty() == true);
	lLeafNode.push_back(m_nodeRoot);
	
	// split a cluster at each iteration until the desired number of clusters is reached
	while((int)lLeafNode.size() < m_iBaseClasses) {
	
		RTNode *clusterToSplit = NULL;
		LRTNode::iterator jt;
		
		// (2) k-means clustering
		if (iClusteringMethod == CLUSTERING_METHOD_KMEANS) {
		
			bool bSplit = false;
			do {
			
				clusterToSplit = NULL;
				
				// select the cluster that contributes more to the overall distortion
				float fDistortionHighest = -FLT_MAX;
				for(LRTNode::iterator it = lLeafNode.begin() ; it != lLeafNode.end() ; ++it) {
				
					// if we already failed to cluster the node, skip it
					if ((*it)->bTouched == false) {
						continue;
					}				
				
					// skip clusters with not enough Gaussian distributions
					if ((*it)->vGaussian.size() < 2*iMinimumComponentsCluster) {
						continue;
					}	
					
					// keep the cluster with the highest distortion
					if ((*it)->fDistortion > fDistortionHighest) {
						fDistortionHighest = (*it)->fDistortion;
						clusterToSplit = (*it);
						jt = it;
					}
				}	
				if (clusterToSplit != NULL) {
					// try to split the cluster using k-means (k=2), maybe not enough Gaussian distributions for the clustering
					// note that this also depends on the clusters initialization
					if (kMeansClustering(clusterToSplit,iMinimumComponentsCluster) == false) {
						clusterToSplit->bTouched = false;
					} else {
						bSplit = true;
						break;
					}
				} else {
					break;
				}
				
			} while(bSplit == false);
			
			// no cluster to split: end!
			if (clusterToSplit == NULL) {
				BVC_ERROR << "clustering failed to produce " << iBaseClasses << "base-classes";
			}
		} 
		// (2') EM clustering
		else {		
			
			// (2.1) sort the nodes in the list by the number of samples
			lLeafNode.sort(RegressionTree::compareSize);
			
			// (2.2) try to split the heaviest node first
			bool bSplitSuccessful = false;
			for(LRTNode::iterator it = lLeafNode.begin() ; it != lLeafNode.end() ; ++it) {
			
				//printf("size: %d\n",(*it)->vGaussian.size());
				
				// check that the number of samples is enough to split the node into two clusters
				if ((*it)->vGaussian.size() < 2*iMinimumComponentsCluster) {
					break;
				}
				
				// if we alread failed to cluster the node, skip it
				if ((*it)->bTouched == false) {
					continue;
				}
			
				// split the cluster using EM-clustering (2 classes) 
				clusterToSplit = *it;
				jt = it;
				if (expectationMaximizationClustering(clusterToSplit,iMinimumComponentsCluster) == true) {
					assert((clusterToSplit->left->vGaussian.size() > iMinimumComponentsCluster) && 
						(clusterToSplit->right->vGaussian.size() > iMinimumComponentsCluster));
					bSplitSuccessful = true;
					break;	
				} 
				// mark the leaf-node as unable to be clustered
				else {
					(*it)->bTouched = false;
				}
			}
			if (bSplitSuccessful == false) {
				BVC_ERROR << "clustering failed to produce " << iBaseClasses << "base-classes";
			}
		}		
		
		assert(clusterToSplit->left->vGaussian.size() >= iMinimumComponentsCluster);
		assert(clusterToSplit->right->vGaussian.size() >= iMinimumComponentsCluster);
		
		// remove the split cluster and insert the cluster resulting form the split
		assert((*jt)->vGaussian.size() == 0);
		lLeafNode.erase(jt);	
		lLeafNode.push_back(clusterToSplit->left);
		lLeafNode.push_back(clusterToSplit->right);	
		
		fDistortionTotal -= clusterToSplit->fDistortion;
		fDistortionTotal += clusterToSplit->left->fDistortion;
		fDistortionTotal += clusterToSplit->right->fDistortion;
		cout << "leaves: " << setw(4) << lLeafNode.size() << " (" << setw(6) 
			<< clusterToSplit->left->vGaussian.size() << " " << setw(6) 
			<< clusterToSplit->right->vGaussian.size() << ") distortion: " << setw(12) << fDistortionTotal 
			<< " Gaussian components: " << setw(7) 
			<< clusterToSplit->left->vGaussian.size()+clusterToSplit->right->vGaussian.size() << endl;
	}
	
	// sanity checks
	int iTreeLeaves = countLeaves(m_nodeRoot);
	assert(iTreeLeaves == (int)lLeafNode.size());
	assert(iTreeLeaves == m_iBaseClasses);
	
	for(LRTNode::iterator it = lLeafNode.begin() ; it != lLeafNode.end() ; ++it) {
		assert((*it)->vGaussian.size() >= iMinimumComponentsCluster);	
	}	
	
	// store the base-class in each of the gaussians and empty the vector of gaussians
	// create a vector (a list does not allow random access) containing the base classes
	m_baseClasses = new RTNode*[m_iBaseClasses];
	int iBaseClass = 0;
	for(LRTNode::iterator it = lLeafNode.begin() ; it != lLeafNode.end() ; ++it, ++iBaseClass) {
		(*it)->iBaseClass = iBaseClass;
		for(VGaussianStats::iterator jt = (*it)->vGaussianStats.begin() ; jt != (*it)->vGaussianStats.end() ; ++jt) {
			m_iGaussianBaseClass[(*jt)->gaussian->iId] = iBaseClass;
		}
		m_baseClasses[iBaseClass] = (*it);
		// allocate memory for the statistics
		(*it)->matrixObservationTranspose = new Matrix<double>(m_iDim);
	}	
	
	// ending time
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dMillisecondsInterval = dTimeEnd - dTimeStart;
	
	// generate graph statistics
	float fDepthAverage = 0.0;
	int iDepthMinimum = INT_MAX;	
	int iDepthMaximum = 0;	
	float fClusterSizeAverage = 0.0;
	unsigned int iClusterSizeMinimum = INT_MAX;
	unsigned int iClusterSizeMaximum = 0;
	for(LRTNode::iterator it = lLeafNode.begin() ; it != lLeafNode.end() ; ++it) {
		fClusterSizeAverage += (*it)->vGaussian.size();
		if ((*it)->vGaussian.size() > iClusterSizeMaximum) {
			iClusterSizeMaximum = (*it)->vGaussian.size();
		}
		if ((*it)->vGaussian.size() < iClusterSizeMinimum) {
			iClusterSizeMinimum = (*it)->vGaussian.size();
		}
		fDepthAverage += (*it)->iDepth;
		if ((*it)->iDepth > iDepthMaximum) {
			iDepthMaximum = (*it)->iDepth;
		}
		if ((*it)->iDepth < iDepthMinimum) {
			iDepthMinimum = (*it)->iDepth;
		}
	}
	fDepthAverage /= (float)m_iBaseClasses;
	fClusterSizeAverage /= (float)m_iBaseClasses;
	
	// show tree building statistics
	printf("-------------------------------------------------------------\n");
	printf("# HMM-states:          %10d\n",m_iHMMStates);
	printf("# Gaussian components: %10d\n",m_hmmManager->getNumberGaussianComponents());
	printf("# tree-nodes:          %10d\n",m_iNodes);
	printf("# base-classes:        %10d\n",m_iBaseClasses);
	printf("average depth:         %10.2f (min= %5d max= %5d)\n",fDepthAverage,iDepthMinimum,iDepthMaximum);
	printf("average cluster size:  %10.2f (min= %5d max= %5d)\n",fClusterSizeAverage,iClusterSizeMinimum,iClusterSizeMaximum); 
   printf("building time:         %10.2f (s)\n",dMillisecondsInterval/1000.0);	
	printf("-------------------------------------------------------------\n");
 
   // clean-up
   // remove all the Gaussian distrib from the leave nodes, once the tree is built the data
   // structure is used to keep the gaussians that contain adaptation data
   for(int i=0 ; i<m_iBaseClasses ; ++i) {
		m_baseClasses[i]->vGaussian.clear();
	}
}

// destroy the regression tree
void RegressionTree::destroy(RTNode *node) {

	if (node->left) {
		destroy(node->left);
		destroy(node->right);
	}
	delete node->matrixObservationTranspose;
	delete node;
}

// counts the number of leaves in the tree (recursively)
int RegressionTree::countLeaves(RTNode *rtNode) {

	if (rtNode->left == NULL) {
		assert(rtNode->right == NULL);
		return 1;
	} else {
		return countLeaves(rtNode->left)+countLeaves(rtNode->right);
	}
}

// perform k-means to split a tree-node
// note: it returns the new tree-nodes as children of the given tree-node
bool RegressionTree::kMeansClustering(RTNode *cluster, int iMinimumComponentsCluster) {

	float *fCentroid1 = new float[m_iDim];
	float *fCentroid2 = new float[m_iDim];

	// (1) get the initial centroids
	// (1.1) choose a random point as the first centroid
	int iGaussians = cluster->vGaussian.size();
	int iIndex1 = getRandomNumber(0,iGaussians-1);
	for(int i=0 ; i<m_iDim ; ++i) {
		fCentroid1[i] = cluster->vGaussian[iIndex1]->fMean[i];
	}
	
	// (1.2) choose a different point as the second centroid
	float *fWeight = computeWeight(cluster);
	float *fCentroidAux = NULL;
	int iIndex2;
	do {
		iIndex2 = getRandomNumber(0,iGaussians-1);
	} while(iIndex1 == iIndex2);
	fCentroidAux = cluster->vGaussian[iIndex2]->fMean;	
	assert(fCentroidAux != NULL);
	for(int i=0 ; i<m_iDim ; ++i) {
		fCentroid2[i] = fCentroidAux[i];
	}
	
	int iIterations = 0;
	float fDistortion1 = 0.0;
	float fDistortion2 = 0.0;
	float fDistortionPrevious = 0.0;
	float fDistortionCurrent = 0.0;
	double *dAccumulator1 = new double[m_iDim];
	double *dAccumulator2 = new double[m_iDim];
	int iClusterElements1 = 0;
	int iClusterElements2 = 0;	
	
	// (2) iterative process assignment/update
	do {
		// initialization
		fDistortion1 = 0.0;
		fDistortion2 = 0.0;
		for(int i=0;i<m_iDim;++i) {
			dAccumulator1[i] = 0.0;
			dAccumulator2[i] = 0.0;
		}
		iClusterElements1 = 0;
		iClusterElements2 = 0;
		fDistortionPrevious = fDistortionCurrent;
		// assignment step: assign each gaussian to the closest cluster
		for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
			float fDistance1 = euclideanDistance(fCentroid1,(*it)->fMean,fWeight);
			float fDistance2 = euclideanDistance(fCentroid2,(*it)->fMean,fWeight);
			if (fDistance1 < fDistance2) {	
				m_iGaussianBaseClass[(*it)->iId] = 0;
				fDistortion1 += fDistance1;
				for(int i=0;i<m_iDim;++i) {
					dAccumulator1[i] += (*it)->fMean[i];
				}	
				iClusterElements1++;
			} else {
				m_iGaussianBaseClass[(*it)->iId] = 0;
				fDistortion2 += fDistance2;
				for(int i=0;i<m_iDim;++i) {
					dAccumulator2[i] += (*it)->fMean[i];
				}	
				iClusterElements2++;
			}	
		}
		// update step: update the centroids
		for(int i=0;i<m_iDim;++i) {
			fCentroid1[i] = (float)(dAccumulator1[i]/((float)iClusterElements1));
			fCentroid2[i] = (float)(dAccumulator2[i]/((float)iClusterElements2));
		}		
		++iIterations;
		//printf("iteration: %d distortion: %f (%d %d)\n",iIterations,fDistortion1+fDistortion2,iClusterElements1,iClusterElements2);
		fDistortionCurrent = fDistortion1+fDistortion2;
	} while(fDistortionPrevious != fDistortionCurrent);
	
	//printf("iteration: %d distortion: %f (%d %d)\n",iIterations,fDistortion1+fDistortion2,iClusterElements1,iClusterElements2);
	
	// check that both clusters have enough number of elements
	if ((iClusterElements1 < iMinimumComponentsCluster) || 
		(iClusterElements2 < iMinimumComponentsCluster)) {
		return false;
	}	
	
	// create the new clusters
	cluster->left = new RTNode;
	++m_iNodes;	
	cluster->left->iId = m_iId++;	
	cluster->left->iBaseClass = -1;	
	cluster->left->fOccupation = 0.0;
	cluster->left->fDistortion = fDistortion1;
	cluster->left->iDepth = cluster->iDepth+1;
	cluster->left->bTouched = true;
	cluster->left->left = NULL; 
	cluster->left->right = NULL; 	
	cluster->left->matrixObservationTranspose = NULL;
	cluster->left->transform = NULL; 	
	cluster->right = new RTNode();
	++m_iNodes;
	cluster->right->iId = m_iId++;
	cluster->right->iBaseClass = -1;	
	cluster->right->fOccupation = 0.0;
	cluster->right->fDistortion = fDistortion2;
	cluster->right->iDepth = cluster->iDepth+1;
	cluster->right->bTouched = true;
	cluster->right->left = NULL; 
	cluster->right->right = NULL; 
	cluster->right->matrixObservationTranspose = NULL;
	cluster->right->transform = NULL;
	// move the gaussians to the corresponding cluster
	for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		if (m_iGaussianBaseClass[(*it)->iId] == 0) {
			cluster->left->vGaussian.push_back(*it);	
		} else {
			cluster->right->vGaussian.push_back(*it);
		}
	}
	cluster->vGaussian.clear();
	
	delete [] fWeight;
	delete [] fCentroid1;
	delete [] fCentroid2;
	delete [] dAccumulator1;
	delete [] dAccumulator2;
		
	return true;
}

// compute P(x|A)	
float RegressionTree::computeClusterLikelihood(float *fMean, float *fCovariance, float *fX) {
	
	double dAcc,dDeterminant;
	float fProbability = 0.0;
	
	double dConstant = pow(PI_NUMBER*2.0,((double)m_iDim)/2.0);
	
	dDeterminant = 1.0;
	double dExponent = 0.0;
	
	for(int j = 0 ; j < m_iDim ; ++j) {
		// compute the determinant of the covariance matrix
		dDeterminant *= fCovariance[j];
		// compute the numerator
		dAcc = fX[j]-fMean[j];
		dExponent += dAcc*dAcc*(1.0/fCovariance[j]);
	}		
	dDeterminant = pow(dDeterminant,0.5);
	fProbability = (float)(exp(-0.5*dExponent)/(dConstant*dDeterminant));
	
	// floor the likelihood
	if (fProbability < EM_MINIMUM_LIKELIHOOD) {
		fProbability = EM_MINIMUM_LIKELIHOOD;
	} else if (fProbability > EM_MAXIMUM_LIKELIHOOD) {
		fProbability = EM_MAXIMUM_LIKELIHOOD;
	}
	
	if (finite(fProbability) == 0.0) {
		fProbability = EM_MINIMUM_LIKELIHOOD;
	}
	
	return fProbability;
}

// perform EM clustering to split a tree-node into two subnodes
// note: it returns the new tree-nodes as children of the given tree-node
bool RegressionTree::expectationMaximizationClustering(RTNode *cluster, int iMinimumComponentsCluster) {

	// (1) k-means clustering (it is used to initialize the clusters in order to speed up the EM convergence)	
	float *fCentroid1 = new float[m_iDim];
	float *fCentroid2 = new float[m_iDim];

	// (1.1) get the initial centroids
	// choose a random point as the first centroid
	int iGaussians = cluster->vGaussian.size();
	int iIndex1 = getRandomNumber(0,iGaussians-1);
	for(int i=0 ; i<m_iDim ; ++i) {
		fCentroid1[i] = cluster->vGaussian[iIndex1]->fMean[i];
	}
	
	// choose a different point as the second centroid
	float *fWeight = computeWeight(cluster);
	float *fCentroidAux = NULL;
	int iIndex2;
	do {
		iIndex2 = getRandomNumber(0,iGaussians-1);
	} while(iIndex1 == iIndex2);	
	fCentroidAux = cluster->vGaussian[iIndex2]->fMean;	
	assert(fCentroidAux != NULL);
	for(int i=0 ; i<m_iDim ; ++i) {
		fCentroid2[i] = fCentroidAux[i];
	}
	
	int iIterations = 0;
	float fDistortion1 = 0.0;
	float fDistortion2 = 0.0;
	float fDistortionPrevious = 0.0;
	float fDistortionCurrent = 0.0;
	double *dAccumulator1 = new double[m_iDim];
	double *dAccumulator2 = new double[m_iDim];
	int iClusterElements1 = 0;
	int iClusterElements2 = 0;	
	
	// (1.2) iterative process assignment/update
	do {
		// initialization
		fDistortion1 = 0.0;
		fDistortion2 = 0.0;
		for(int i=0;i<m_iDim;++i) {
			dAccumulator1[i] = 0.0;
			dAccumulator2[i] = 0.0;
		}
		iClusterElements1 = 0;
		iClusterElements2 = 0;
		fDistortionPrevious = fDistortionCurrent;
		// assignment step: assign each gaussian to the closest cluster
		for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
			float fDistance1 = euclideanDistance(fCentroid1,(*it)->fMean,fWeight);
			float fDistance2 = euclideanDistance(fCentroid2,(*it)->fMean,fWeight);
			if (fDistance1 < fDistance2) {
				m_iGaussianBaseClass[(*it)->iId] = 0;
				fDistortion1 += fDistance1;
				for(int i=0;i<m_iDim;++i) {
					dAccumulator1[i] += (*it)->fMean[i];
				}	
				iClusterElements1++;
			} else {
				m_iGaussianBaseClass[(*it)->iId] = 1;
				fDistortion2 += fDistance2;
				for(int i=0;i<m_iDim;++i) {
					dAccumulator2[i] += (*it)->fMean[i];
				}	
				iClusterElements2++;
			}	
		}
		// update step: update the centroids
		for(int i=0;i<m_iDim;++i) {
			fCentroid1[i] = (float)(dAccumulator1[i]/((float)iClusterElements1));
			fCentroid2[i] = (float)(dAccumulator2[i]/((float)iClusterElements2));
		}		
		++iIterations;
		//printf("iteration: %d distortion: %f (%d %d)\n",iIterations,fDistortion1+fDistortion2,iClusterElements1,iClusterElements2);
		fDistortionCurrent = fDistortion1+fDistortion2;
		
	} while(fDistortionPrevious != fDistortionCurrent);
	
	delete [] fCentroid1;
	delete [] fCentroid2;
	delete [] fWeight;
	
	// (2) initialize the mean and covariance to the cluster's, and compute the priors using the cluster elements count

	// each cluster is modeled by a gaussian 
	float *fMean1 = new float[m_iDim];	
	float *fMean2 = new float[m_iDim];
	float *fCovariance1 = new float[m_iDim];	
	float *fCovariance2 = new float[m_iDim];
	
	// compute the initial priors
	float fPrior1 = ((float)iClusterElements1)/((float)(iClusterElements1+iClusterElements2));
	float fPrior2 = 1.0f-fPrior1;
	
	double *dObservation1 = new double[m_iDim];
	double *dObservation2 = new double[m_iDim];
	double *dObservationSquare1 = new double[m_iDim];
	double *dObservationSquare2 = new double[m_iDim];
	
	// reset
	for(int i=0 ; i<m_iDim ;++i) {
		dObservation1[i] = 0.0;
		dObservationSquare1[i] = 0.0;
		dObservation2[i] = 0.0;
		dObservationSquare2[i] = 0.0;
	}

	// accumulate statistics
	int iElements1 = 0;
	int iElements2 = 0;
	for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		if (m_iGaussianBaseClass[(*it)->iId] == 0) {
			for(int i=0 ; i<m_iDim ;++i) {
				dObservation1[i] += fPrior1*(*it)->fMean[i];
				dObservationSquare1[i] += fPrior1*(*it)->fMean[i]*(*it)->fMean[i];		
			}
			++iElements1;
		} else {
			for(int i=0 ; i<m_iDim ;++i) {
				dObservation2[i] += fPrior2*(*it)->fMean[i];
				dObservationSquare2[i] += fPrior2*(*it)->fMean[i]*(*it)->fMean[i];
			}
			++iElements2;
		}		
	}	
	
	// compute the initial clusters mean and covariance
	for(int i=0 ; i<m_iDim ;++i) {
		fMean1[i] = (float)(dObservation1[i]/((float)iElements1));
		fCovariance1[i] = (float)(dObservationSquare1[i]/((float)iElements1))-(fMean1[i]*fMean1[i]);
		fMean2[i] = (float)(dObservation2[i]/((float)iElements2));
		fCovariance2[i] = (float)(dObservationSquare2[i]/((float)iElements2))-(fMean2[i]*fMean2[i]);
	}	
	
	// (3) actual clustering: Expectation Maximization (iterative process assignment/update)
	
	float fLikelihoodGlobalPrev = -FLT_MAX;
	float fLikelihoodGlobal = -FLT_MAX;
	float fLikelihoodGlobal1 = 0.0;
	float fLikelihoodGlobal2 = 0.0;
	
	iIterations = 0;
	do {
	
		fLikelihoodGlobalPrev = fLikelihoodGlobal;
		fLikelihoodGlobal = 0.0;
		fLikelihoodGlobal1 = 0.0;
		fLikelihoodGlobal2 = 0.0;
	
		// initialization
		float fOccupation1 = 0.0;
		float fOccupation2 = 0.0;
		
		iElements1 = 0;
		iElements2 = 0;
		
		for(int i=0 ; i<m_iDim ;++i) {
			dObservation1[i] = 0.0;
			dObservationSquare1[i] = 0.0;
			dObservation2[i] = 0.0;
			dObservationSquare2[i] = 0.0;
		}	
		
		float fResponsibility1;
		float fResponsibility2;

		// assignment step: there is a soft assignment of each gaussian to every cluster
		for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		
			// compute the log-likelihood of the point
			float fLikelihood1 = computeClusterLikelihood(fMean1,fCovariance1,(*it)->fMean);
			float fLikelihood2 = computeClusterLikelihood(fMean2,fCovariance2,(*it)->fMean);	
		
			// compute the probability that the point (a gaussian mean) belongs to each cluster
			fResponsibility1 = (fPrior1*fLikelihood1)/(fPrior1*fLikelihood1+fPrior2*fLikelihood2);
			fResponsibility2 = 1.0f-fResponsibility1;
			
			assert(fResponsibility1 <= 1.0f);
			assert(fResponsibility1 >= 0.0f);
			
			// accumulate data
			for(int i=0 ; i<m_iDim ; ++i) {
				dObservation1[i] += fResponsibility1*(*it)->fMean[i];
				dObservationSquare1[i] += fResponsibility1*((*it)->fMean[i]*(*it)->fMean[i]);
				dObservation2[i] += fResponsibility2*(*it)->fMean[i];
				dObservationSquare2[i] += fResponsibility2*((*it)->fMean[i]*(*it)->fMean[i]);
			}
			
			// accumulate occupation
			fOccupation1 += fResponsibility1;
			fOccupation2 += fResponsibility2;
			
			// compute the likelihood of the datapoint
			float fLikelihoodPoint = log(fPrior1*fLikelihood1+fPrior2*fLikelihood2);
			
			if (fResponsibility1 > fResponsibility2) {
				fLikelihoodGlobal1 += fLikelihoodPoint;
				m_iGaussianBaseClass[(*it)->iId] = 0;
				++iElements1;
			} else {
				fLikelihoodGlobal2 += fLikelihoodPoint;
				m_iGaussianBaseClass[(*it)->iId] = 1;
				++iElements2;
			}
			fLikelihoodGlobal += fLikelihoodPoint;
		}
		
		//printf("global likelihood: %12.2f\n",fLikelihoodGlobal);
		
		// update step
		for(int i=0 ; i<m_iDim ;++i) {
			fMean1[i] = (float)(dObservation1[i]/fOccupation1);
			fCovariance1[i] = (float)(dObservationSquare1[i]/fOccupation1-(fMean1[i]*fMean1[i]));
			fMean2[i] = (float)(dObservation2[i]/fOccupation2);
			fCovariance2[i] = (float)(dObservationSquare2[i]/fOccupation2-(fMean2[i]*fMean2[i]));
		}	
		fPrior1 = fOccupation1/(fOccupation1+fOccupation2);
		fPrior2 = 1.0f-fPrior1;	
		
		++iIterations;
		
		if (iIterations >= EM_MAXIMUM_ITERATIONS) {
			break;
		}
		
	} while ((iIterations < 5) || (fabs(fLikelihoodGlobal-fLikelihoodGlobalPrev) > fabs(fLikelihoodGlobal*EM_CONVERGENCE_EPSILON)));
	
	// clean-up
	delete [] dObservation1;
	delete [] dObservationSquare1;
	delete [] dObservation2;
	delete [] dObservationSquare2;
	delete [] dAccumulator1;
	delete [] dAccumulator2;
	
	delete [] fMean1;
	delete [] fMean2;
	delete [] fCovariance1;
	delete [] fCovariance2;
	
	// check that both clusters have enough number of elements
	if ((iElements1 < iMinimumComponentsCluster) || 
		(iElements2 < iMinimumComponentsCluster)) {
		return false;
	}
	
	// create the new clusters
	cluster->left = new RTNode;
	++m_iNodes;	
	cluster->left->iId = m_iId++;	
	cluster->left->iBaseClass = -1;	
	cluster->left->fOccupation = 0.0;
	cluster->left->fDistortion = fLikelihoodGlobal1;
	cluster->left->iDepth = cluster->iDepth+1;
	cluster->left->bTouched = true;
	cluster->left->left = NULL; 
	cluster->left->right = NULL; 	
	cluster->left->matrixObservationTranspose = NULL;
	cluster->left->transform = NULL; 	
	cluster->right = new RTNode();
	++m_iNodes;
	cluster->right->iId = m_iId++;
	cluster->right->iBaseClass = -1;	
	cluster->right->fOccupation = 0.0;
	cluster->right->fDistortion = fLikelihoodGlobal2;
	cluster->right->iDepth = cluster->iDepth+1;
	cluster->right->bTouched = true;
	cluster->right->left = NULL; 
	cluster->right->right = NULL; 
	cluster->right->matrixObservationTranspose = NULL;
	cluster->right->transform = NULL;
	// move the gaussians to the corresponding cluster
	for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		if (m_iGaussianBaseClass[(*it)->iId] == 0) {
			cluster->left->vGaussian.push_back(*it);	
		} else {
			cluster->right->vGaussian.push_back(*it);
		}
	}
	cluster->vGaussian.clear();
	
	cout << "global likelihood: " << setw(12) << fLikelihoodGlobal 
		<< " (" << cluster->left->vGaussian.size() << " " 
		<< cluster->right->vGaussian.size() << ") iterations: " << iIterations << endl;
		
	return true;
}


// compute the weight vector for the computation of the weighted euclidean distance
float *RegressionTree::computeWeight(RTNode *cluster) {

	// (1) compute the average of the covariance of all the gaussians
	double *fCovAverage = new double[m_iDim];
	for(int i=0;i<m_iDim;++i) {
		fCovAverage[i] = 0.0;
	}
	for(VGaussianDecoding::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		for(int i=0;i<m_iDim;++i) {
			fCovAverage[i] += (*it)->fCovariance[i]; 
		}
	}
	for(int i=0;i<m_iDim;++i) {
		fCovAverage[i] /= ((float)cluster->vGaussian.size());
	}
	
	// (2) compute the inverse square root
	float *fWeight = new float[m_iDim];
	for(int i=0;i<m_iDim;++i) {
		fWeight[i] = (float)(1.0f/sqrt(fCovAverage[i]));
	}	
	
	delete [] fCovAverage;
	
	return fWeight;
}

// print the regression tree
void RegressionTree::print() {

	QRTNode qNode;
	
	qNode.push(m_nodeRoot);
	
	printf("Regression tree -----------------------------------\n");
	while(qNode.empty() == false) {
		RTNode *node = qNode.front();
		qNode.pop();
		// add the children to the queue (if any)
		if (node->left != NULL) {
			assert(node->right != NULL);
			qNode.push(node->left);
			qNode.push(node->right);
			// print the node
			printf("internal node %4d (%4d)(%4d) depth= %d occ= %8.2f\n",node->iId,node->left->iId,node->right->iId,node->iDepth,node->fOccupation);	
		} else {
			// print the node
			printf("leaf node     %4d depth= %d occ= %8.2f \n",node->iId,node->iDepth,node->fOccupation);
		}
	}
	printf("---------------------------------------------------\n");
}

// store the tree into disk
void RegressionTree::store(const char *strFile) {

	FileOutput file(strFile,true);
	file.open();

	IOBase::write(file.getStream(),m_iHMMStates);
	IOBase::write(file.getStream(),m_iGaussians);
	IOBase::write(file.getStream(),m_iNodes);				// tree nodes
	IOBase::write(file.getStream(),m_iBaseClasses);
	
	// write the nodes one by one	
	QRTNode qNode;	
	qNode.push(m_nodeRoot);	
	while(qNode.empty() == false) {
		RTNode *node = qNode.front();
		qNode.pop();
		int iNodeLeftId = -1;
		int iNodeRightId = -1;
		unsigned char iPosition = NODE_POSITION_LEAF;
		// add the children to the queue (if any)	
		if (node->left != NULL) {
			assert(node->right != NULL);
			qNode.push(node->left);
			qNode.push(node->right);
			iNodeLeftId = node->left->iId;
			iNodeRightId = node->right->iId;
			iPosition = NODE_POSITION_INTERNAL;
		}	
		// write the node information
		IOBase::write(file.getStream(),node->iId);
		IOBase::write(file.getStream(),iPosition);
		// internal node
		if (iPosition == NODE_POSITION_INTERNAL) {
			IOBase::write(file.getStream(),iNodeLeftId);
			IOBase::write(file.getStream(),iNodeRightId);
		}
		// leaf node
		else {
			IOBase::write(file.getStream(),node->iBaseClass);
		}
	}
	
	// write the mapping from Gaussian components to base classes
	int iHMMStates = -1;
	HMMStateDecoding *hmmStates = m_hmmManager->getHMMStatesDecoding(&iHMMStates);
	for(int i=0 ; i < iHMMStates ; ++i) {
		// state-id
		int iId = hmmStates[i].getId();
		IOBase::write(file.getStream(),iId);
		// # Gaussian components
		int iComponents = hmmStates[i].getGaussianComponents();
		IOBase::write(file.getStream(),iComponents);
		// base-class for each Gaussian component
		for(int j=0 ; j<hmmStates[i].getGaussianComponents() ; ++j) {
			int iBaseClass = m_iGaussianBaseClass[hmmStates[i].getGaussian(j)->iId];
			IOBase::write(file.getStream(),iBaseClass);
		}
	}

	file.close();
}

// load the tree from disk
void RegressionTree::load(const char *strFile) {

	FileInput file(strFile,true);
	file.open();	
	
	IOBase::read(file.getStream(),&m_iHMMStates);
	IOBase::read(file.getStream(),&m_iGaussians);
	IOBase::read(file.getStream(),&m_iNodes);				// tree nodes
	assert(m_iNodes > 0);
	IOBase::read(file.getStream(),&m_iBaseClasses);
	
	RTNode **nodes = new RTNode*[m_iNodes];
	// allocate memory for the nodes
	for(int i=0 ; i<m_iNodes ; ++i) {
		nodes[i] = new RTNode;
	}
	m_nodeRoot = nodes[0];
	
	// read the nodes one by one	
	LRTNode lLeafNode;
	int iNodeLeftId = -1;
	int iNodeRightId = -1;
	for(int i=0 ; i<m_iNodes ; ++i) {
		// read the node information
		// node id
		int iId;
		IOBase::read(file.getStream(),&iId);
		nodes[iId]->iId = iId;
		// position
		unsigned char iPosition;
		IOBase::read(file.getStream(),&iPosition);
		// internal node
		if (iPosition == NODE_POSITION_INTERNAL) {
			IOBase::read(file.getStream(),&iNodeLeftId);
			IOBase::read(file.getStream(),&iNodeRightId);
			assert((iNodeLeftId >= 0) && (iNodeLeftId < m_iNodes));
			assert((iNodeRightId >= 0) || (iNodeRightId < m_iNodes));
			nodes[iId]->iBaseClass = -1;
			nodes[iId]->left = nodes[iNodeLeftId];
			nodes[iId]->right = nodes[iNodeRightId];
		}
		// leaf node
		else {
			// base-class
			IOBase::read(file.getStream(),&nodes[iId]->iBaseClass);
			nodes[iId]->left = NULL;
			nodes[iId]->right = NULL;
			lLeafNode.push_back(nodes[iId]);
		}
		// initialize the remaining fields
		nodes[iId]->fOccupation = 0.0;
		nodes[iId]->fDistortion = 0.0;
		nodes[iId]->transform = NULL;	
		nodes[iId]->matrixObservationTranspose = NULL;
	}	
	delete [] nodes;
	
	// load the mapping from Gaussian components to base classes
	int iHMMStates = -1;
	HMMStateDecoding *hmmStates = m_hmmManager->getHMMStatesDecoding(&iHMMStates);
	for(int i=0 ; i < iHMMStates ; ++i) {
		// state-id
		int iId = -1;
		IOBase::read(file.getStream(),&iId);
		assert(iId == hmmStates[i].getId());
		// # Gaussian components
		int iComponents = -1;
		IOBase::read(file.getStream(),&iComponents);
		assert(iComponents == hmmStates[i].getGaussianComponents());
		// base-class for each Gaussian component
		for(int j=0 ; j<iComponents ; ++j) {
			int iBaseClass = -1;
			IOBase::read(file.getStream(),&iBaseClass);
			assert((iBaseClass >= 0) && (iBaseClass < m_iBaseClasses));
			m_iGaussianBaseClass[hmmStates[i].getGaussian(j)->iId] = iBaseClass;
		}
	}
	
	file.close();
	
	// store the base-class in each of the gaussians and empty the vector of gaussians
	// create a vector (a list does not allow random access) containing the base classes
	m_baseClasses = new RTNode*[m_iBaseClasses];
	for(LRTNode::iterator it = lLeafNode.begin() ; it != lLeafNode.end() ; ++it) {
		m_baseClasses[(*it)->iBaseClass] = (*it);
		// allocate memory for the statistics
		(*it)->matrixObservationTranspose = new Matrix<double>(m_iDim);
		(*it)->matrixObservationTranspose->zero();	
	}	
}

// compute a transform from the adaptation data (the given vector of gaussians)
void RegressionTree::computeTransform(MLLRTransform *transform, VGaussianStats &vGaussianStats,
	bool bMeanOnly) {

	// (1) mean transform
	transform->matrixMean = new Matrix<float>(m_iDim,m_iDim+1);
	transform->matrixMean->zero();
	
	// allocate memory for all the k(i)
	if (transform->matrixK == NULL) {	
		transform->matrixK = new Matrix<float>(m_iDim,m_iDim+1);
		transform->matrixK->zero();	
		assert(transform->matrixG == NULL);
		transform->matrixG = new Matrix<float>*[m_iDim];
		for(int i=0 ; i < m_iDim ; ++i) {
			transform->matrixG[i] = new Matrix<float>(m_iDim+1,m_iDim+1);
		}
	}
	
	// compute each of the k(i) (rows in the K matrix)
	Matrix<float> matrixTemp(m_iDim,m_iDim+1);
	for(VGaussianStats::iterator it = vGaussianStats.begin() ; it != vGaussianStats.end() ; ++it) {
		// divide first order stats by covariance
		VectorStatic<float> vCovariance((*it)->gaussian->fCovariance,m_iDim); 
		Vector<float> vColumn(m_iDim);	
		vColumn.copy(*(*it)->vObservation);
		vColumn.divide(vCovariance);
		// copy to each column in the matrix
		matrixTemp.copyCols(vColumn);
		// create the extended mean vector
		VectorStatic<float> vMean((*it)->gaussian->fMean,m_iDim);
		Vector<float> vMeanEx(vMean);	
		vMeanEx.appendFront(1.0);
		// multiply each row by the extended mean
		matrixTemp.mulRows(vMeanEx);
		// accumulate
		transform->matrixK->add(1.0,matrixTemp);
	}	
	
	for(int i=0 ; i < m_iDim ; ++i) {

		// compute G(i)
		Matrix<float> *matrixG = transform->matrixG[i];
		matrixG->zero();
		for(VGaussianStats::iterator it = vGaussianStats.begin() ; it != vGaussianStats.end() ; ++it) {
			
			// create the extended mean vector
			VectorStatic<float> vMean((*it)->gaussian->fMean,m_iDim);
			Vector<float> vMeanEx(vMean);	
			vMeanEx.appendFront(1.0);
				
			matrixG->addVecMul((float)((*it)->dOccupation/(*it)->gaussian->fCovariance[i]),vMeanEx,vMeanEx);
		}
			
		// compute the ith row of the transformation matrix
		Matrix<float> matrixGInverted(*matrixG);
		matrixGInverted.invert();
		transform->matrixMean->getRow(i).mul(transform->matrixK->getRow(i),matrixGInverted);
	}
	
	// (2) covariance transform (optional)
	if (bMeanOnly == false) {
	
		// allocate memory for the covariance transform
		transform->matrixCovariance = new Matrix<double>(m_iDim);
		transform->matrixCovariance->zero();	
		float fDenominatorH = 0.0;
		
		// keep the baseclasses whose data is already accumulated (first term of the numerator)
		bool *bAccumulated = new bool[m_iBaseClasses];
		for(int i=0 ; i<m_iBaseClasses ; ++i) {
			bAccumulated[i] = false;
		}
		
		// compute the first term of the numerator: accumulate statistics for each base-class
		for(VGaussianStats::iterator it = vGaussianStats.begin() ; it != vGaussianStats.end() ; ++it) {
			int iBaseClass = m_iGaussianBaseClass[(*it)->gaussian->iId];
			if (bAccumulated[iBaseClass]) {
				continue;
			} 
			bAccumulated[iBaseClass] = true;
			transform->matrixCovariance->add(1.0,*m_baseClasses[iBaseClass]->matrixObservationTranspose);
		}
		
		delete [] bAccumulated;
		
		Matrix<double> matrixAcc(m_iDim);
		
		for(VGaussianStats::iterator it = vGaussianStats.begin() ; it != vGaussianStats.end() ; ++it) {
		
			matrixAcc.zero();
			
			// compute the Cholesky factor of the inverted covariance matrix 
			// note: covariance is converted to standard format
			Vector<double> vC(VectorStatic<float>((*it)->gaussian->fCovariance,m_iDim));
		#ifdef OPTIMIZED_COMPUTATION
			vC.mul(2.0);
			vC.invertElements();
		#endif
			vC.invertElements();
			vC.sqrt();
				
			// compute the adapted mean
			Vector<double> vMeanAdapted(VectorStatic<float>((*it)->gaussian->fMean,m_iDim));
			Vector<double> vMeanEx(vMeanAdapted);	
			vMeanEx.appendFront(1.0);
			vMeanAdapted.mul(*transform->matrixMean,vMeanEx);
			
			// compute the second, third and fourth terms
			matrixAcc.addVecMul(-1.0,vMeanAdapted,*(*it)->vObservation);
			matrixAcc.addVecMul(-1.0,*(*it)->vObservation,vMeanAdapted);
			matrixAcc.addVecMul((*it)->dOccupation,vMeanAdapted,vMeanAdapted);
			matrixAcc.mulRows(vC);
			matrixAcc.mulCols(vC);
			
			transform->matrixCovariance->add(1.0,matrixAcc);
			
			// step 4: denominator, which is the regression class occupation
			fDenominatorH += (float)(*it)->dOccupation;
		}	
		transform->matrixCovariance->mul(1.0/fDenominatorH);
	}
}

// compute MLLR adaptation transforms using the regression tree
void RegressionTree::computeTransforms(float fMinimumOccupationTransform, int iMinimumGaussianComponentsObserved,
	GaussianStats **gaussianStats, bool bMeanOnly) {

	// attach each Gaussian with adaptation data to its corresponding leaf-node (base-class)
	double dOccupationTotal = 0.0;
	int iComponentsOccupied = 0;
	for(int i=0 ; i < m_iGaussians ; ++i) {	
		if (gaussianStats[i]) {
			m_baseClasses[m_iGaussianBaseClass[gaussianStats[i]->gaussian->iId]]->vGaussianStats.push_back(gaussianStats[i]);
			m_baseClasses[m_iGaussianBaseClass[gaussianStats[i]->gaussian->iId]]->fOccupation += (float)gaussianStats[i]->dOccupation;
			dOccupationTotal += gaussianStats[i]->dOccupation;
			++iComponentsOccupied;
		}
	}
	
	// check if enough data to compute a transform
	if ((dOccupationTotal < fMinimumOccupationTransform) ||
		(iComponentsOccupied < iMinimumGaussianComponentsObserved)) {
		BVC_ERROR << "not enough adaptation data to compute a transform";
	}	
	
	// get the transforms that will be computed
	assert(m_nodeRoot != NULL);
	identifyTransformsToCompute(fMinimumOccupationTransform,iMinimumGaussianComponentsObserved,m_nodeRoot);
	
	// assign an id to each transform
	int iId = 0;
	for(VMLLRTransform::iterator it = m_vMLLRTransform.begin() ; it != m_vMLLRTransform.end() ; ++it, ++iId) {
		(*it)->iId = iId;
	}
	
	// print the transforms
	printTransforms();
	
	// initialize the base-class to transform map
	if (m_transformMap != NULL) {
		delete [] m_transformMap;
	}
	m_transformMap = new MLLRTransform*[m_iBaseClasses];
	for(int i=0 ; i < m_iBaseClasses ; ++i) {
		m_transformMap[i] = NULL;
	}
	
	// compute the transforms one by one
	VGaussianStats vGaussianStats;
	QRTNode qNode;
	for(VMLLRTransform::iterator it = m_vMLLRTransform.begin() ; it != m_vMLLRTransform.end() ; ++it) {
		
		// (1)  get the list of Gaussian components that will be used to compute the transform
		qNode.push((*it)->node);
		while(qNode.empty() == false) {
			RTNode *node = qNode.front();
			qNode.pop();
			// if the node has already a tranform, add the k and G matrices
			if ((node->transform) && (node != (*it)->node)) {
				// (a) the current node does not have a transform yet: allocate memory and initialize
				if ((*it)->matrixK == NULL) {
					assert((*it)->matrixG == NULL);
					(*it)->matrixK = new Matrix<float>(m_iDim,m_iDim+1);
					(*it)->matrixG = new Matrix<float>*[m_iDim];
					for(int i=0 ; i < m_iDim ; ++i) {
						(*it)->matrixG[i] = new Matrix<float>(m_iDim+1,m_iDim+1);
					}
				}
				// (b) the current node has already a transform: add matrices 
				else {
					(*it)->matrixK->add(1.0,*node->transform->matrixK);
					for(int i=0 ; i<m_iDim ; ++i) {
						(*it)->matrixG[i]->add(1.0,*node->transform->matrixG[i]);
					}
				}
				// nothing else to do: if there are children they are ignored
			}
			else {
				// the node does not have a transform: copy the gaussians and move the children to the queue
				for(VGaussianStats::iterator jt = node->vGaussianStats.begin() ; jt != node->vGaussianStats.end() ; ++jt) {
					vGaussianStats.push_back(*jt);
				}	
				if (node->left != NULL) {
					qNode.push(node->left);
					qNode.push(node->right);
				} else {
					m_transformMap[node->iBaseClass] = *it;	
				}
			}
		}
		
		// (2) compute the tranform from the list of Gaussian distrib
		computeTransform(*it,vGaussianStats,bMeanOnly);
		vGaussianStats.clear();
	}
	
	// sanity check
	for(int i=0 ; i < m_iBaseClasses ; ++i) {
		assert(m_transformMap[i] != NULL);
	}	
}

// identify the transforms that will be computed from the tree (recursively)
float RegressionTree::identifyTransformsToCompute(float fMinimumOccupationTransform, 
	int iMinimumGaussianComponentsObserved, RTNode *node) {

	// (a) always compute a transform for leaf nodes if the occupation/observed components is big enough
	if (node->left == NULL) {
		assert(node->right == NULL);
		// check if there is enough data to compute a transform	
		if ((node->fOccupation >= fMinimumOccupationTransform) && 
			(getGaussianComponentsObserved(node) >= iMinimumGaussianComponentsObserved)) {
			// create a new transform
			MLLRTransform *transform = newMLLRTransform(node,node->fOccupation);
			m_vMLLRTransform.push_back(transform);
			node->transform = transform;
		}
		return node->fOccupation;	
	}
	
	// (b) compute a transform for the node if either child has not enough occupation
	// get the occupation at child nodes
	float fOccupationLeft = identifyTransformsToCompute(fMinimumOccupationTransform,
		iMinimumGaussianComponentsObserved,node->left);
	float fOccupationRight = identifyTransformsToCompute(fMinimumOccupationTransform,
		iMinimumGaussianComponentsObserved,node->right);
	int iObservedComponentsLeft = getGaussianComponentsObserved(node->left);
	int iObservedComponentsRight = getGaussianComponentsObserved(node->right);	
	if (((fOccupationLeft < fMinimumOccupationTransform) || 
			(fOccupationRight < fMinimumOccupationTransform) ||
			(iObservedComponentsLeft < iMinimumGaussianComponentsObserved) ||
			(iObservedComponentsRight < iMinimumGaussianComponentsObserved)) 
			&& (fOccupationLeft+fOccupationRight >= fMinimumOccupationTransform)
			&& (iObservedComponentsLeft+iObservedComponentsRight >= iMinimumGaussianComponentsObserved)) {
		// create a new transform
		MLLRTransform *transform = newMLLRTransform(node,fOccupationLeft+fOccupationRight);
		m_vMLLRTransform.push_back(transform);
		node->transform = transform;
	}
	
	return fOccupationLeft+fOccupationRight;
}

// apply the transforms to the set of HMM-models
void RegressionTree::applyTransforms(bool bMean, bool bCovariance) {

	for(int i=0 ; i < m_iHMMStates ; ++i) {
		int iGaussians = -1;
		GaussianDecoding *gaussians = m_hmmStates[i].getGaussians(iGaussians);
		for(int j=0 ; j < iGaussians ; ++j) {
			applyTransform(m_transformMap[m_iGaussianBaseClass[gaussians[j].iId]],&gaussians[j],bMean,bCovariance);	
		}
	}
	
	// if the covariances are updated, the precomputed constants and covariance matrix determinant used to speed-up
	// the emission probability computation have to be updated too
	if (bCovariance) {	
		m_hmmManager->initializeDecoding();
	}
}

// print the transforms
void RegressionTree::printTransforms() {

	cout << "Total transforms: " << m_vMLLRTransform.size() << endl;
	for(VMLLRTransform::iterator it = m_vMLLRTransform.begin() ; it != m_vMLLRTransform.end() ; ++it) {
		cout << "node: " << setw(5) << (*it)->node->iId << " occupation: " << setw(10) 
		<< (*it)->dOccupation << " Gaussian components: " << setw(8) 
		<< getGaussianComponentsObserved((*it)->node) << endl;
	}
}

// store the transforms 
void RegressionTree::storeTransforms(const char *strFile) {

	assert(!m_vMLLRTransform.empty());
	
	FileOutput file(strFile,true);
	file.open();
	
	// (1) store the MLLR transforms
	
	int iTransforms = m_vMLLRTransform.size();
	IOBase::write(file.getStream(),iTransforms);	
	
	// mean/covariance?
	unsigned char iType = UCHAR_MAX;
	(m_vMLLRTransform[0]->matrixCovariance == NULL) ? iType = TRANSFORM_MLLR_MEAN : iType = TRANSFORM_MLLR_MEAN_COVARIANCE;
	IOBase::write(file.getStream(),iType);	
		
	// actual transforms
	for(VMLLRTransform::iterator it = m_vMLLRTransform.begin() ; it != m_vMLLRTransform.end() ; ++it) {	

		// transform id
		IOBase::write(file.getStream(),(*it)->iId);
		// mean transform
		(*it)->matrixMean->write(file.getStream());
		// covariance transform
		if (iType == TRANSFORM_MLLR_MEAN_COVARIANCE) {
			(*it)->matrixCovariance->write(file.getStream());
		}		
	}
	
	// (2) store the number of base-classes
	IOBase::write(file.getStream(),m_iBaseClasses);
	
	// (3) store the base-class to transform map
	for(int i = 0 ; i < m_iBaseClasses ; ++i) {
		IOBase::write(file.getStream(),m_transformMap[i]->iId);
	}	
	
	file.close();
}


// load the transforms from a file
void RegressionTree::loadTransforms(const char *strFile) {

	FileInput file(strFile,true);
	file.open();
	
	// (1) load the MLLR transforms
	
	int iTransforms = -1;
	unsigned char iType = UCHAR_MAX;	
	
	IOBase::read(file.getStream(),&iTransforms);	
	IOBase::read(file.getStream(),&iType);

	// actual transforms
	for(int i=0 ; i < iTransforms ; ++i) {	
	
		MLLRTransform *transform = new MLLRTransform;
		
		// read the transform id
		IOBase::read(file.getStream(),&transform->iId);
		// mean transform
		transform->matrixMean = Matrix<float>::read(file.getStream());
		// covariance transform
		transform->matrixCovariance = NULL;
		if (iType == TRANSFORM_MLLR_MEAN_COVARIANCE) {
			transform->matrixCovariance = Matrix<double>::read(file.getStream());
		}
		transform->node = NULL;
		transform->matrixK = NULL;
		transform->matrixG = NULL;	
		m_vMLLRTransform.push_back(transform);
	}
	
	// (2) number of base-classes
	int iBaseClasses = -1;
	IOBase::read(file.getStream(),&iBaseClasses);
	assert(iBaseClasses == m_iBaseClasses);
	
	// (3) base-class to transform map
	m_transformMap = new MLLRTransform*[m_iBaseClasses];
	int iId = -1;
	for(int i = 0 ; i < m_iBaseClasses ; ++i) {
		IOBase::read(file.getStream(),&iId);
		assert((iId >= 0) || (iId < m_iBaseClasses));
		m_transformMap[i] = m_vMLLRTransform[iId];	
	}	
	
	file.close();
}

};	// end-of-namespace

