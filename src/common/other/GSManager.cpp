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


#include "GSManager.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructur
GSManager::GSManager(HMMManager *hmmManager)
{
	m_hmmStatesDecoding = hmmManager->getHMMStatesDecoding(&m_iHMMStates);
	m_clusterRoot = NULL;
}

// destructor
GSManager::~GSManager()
{
}

// perform k-means on a cluster 
// note: it returns the new clusters as children of the given cluster
void GSManager::kMeans(GaussianCluster *cluster) {

	float fCentroid1[DIMENSIONALITY];
	float fCentroid2[DIMENSIONALITY];
	int iAttempts = 10;
	
//start:

	// (1) get the initial couple of centroids
	// (1.1) choose a random point as centroid
	int iGaussians = cluster->vGaussian.size();
	int iIndex1 = getRandomNumber(0,iGaussians-1);
	assert((iIndex1 >= 0) && (iIndex1 < iGaussians));
	for(int i=0 ; i<DIMENSIONALITY ; ++i) {
		fCentroid1[i] = cluster->vGaussian[iIndex1]->gaussian->fMean[i];
	}
	
	// (1.2) choose a different point as the second centroid	
	float *fWeight = computeWeight(cluster);
	float *fCentroidAux = NULL;
	int iIndex2;
	do {
		iIndex2 = getRandomNumber(0,iGaussians-1);
	} while(iIndex1 == iIndex2);
	fCentroidAux = cluster->vGaussian[iIndex2]->gaussian->fMean;
	assert(fCentroidAux != NULL);
	for(int i=0 ; i<DIMENSIONALITY ; ++i) {
		fCentroid2[i] = fCentroidAux[i];
	}
	
	int iIterations = 0;
	float fDistortion1 = 0.0;
	float fDistortion2 = 0.0;
	float fDistortionPrevious = 0.0;
	float fDistortionCurrent = 0.0;
	double dAccumulator1[DIMENSIONALITY];
	double dAccumulator2[DIMENSIONALITY];
	int iClusterElements1 = 0;
	int iClusterElements2 = 0;	
	
	// (2) iterative process assignment/update
	do {
		// initialization
		fDistortion1 = 0.0;
		fDistortion2 = 0.0;
		for(int i=0;i<DIMENSIONALITY;++i) {
			dAccumulator1[i] = 0.0;
			dAccumulator2[i] = 0.0;
		}
		iClusterElements1 = 0;
		iClusterElements2 = 0;
		fDistortionPrevious = fDistortionCurrent;
		// assignment step: assign each gaussian to the closest cluster
		for(VGaussianState::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
			float fDistance1 = weightedEuclideanDistance(fCentroid1,(*it)->gaussian->fMean,fWeight);
			float fDistance2 = weightedEuclideanDistance(fCentroid2,(*it)->gaussian->fMean,fWeight);
			if (fDistance1 < fDistance2) {
				(*it)->iCluster = 0;
				fDistortion1 += fDistance1;
				for(int i=0;i<DIMENSIONALITY;++i) {
					dAccumulator1[i] += (*it)->gaussian->fMean[i];
				}	
				iClusterElements1++;
			} else {
				(*it)->iCluster = 1;
				fDistortion2 += fDistance2;
				for(int i=0;i<DIMENSIONALITY;++i) {
					dAccumulator2[i] += (*it)->gaussian->fMean[i];
				}	
				iClusterElements2++;
			}	
		}
		// update step: update the centroids
		for(int i=0;i<DIMENSIONALITY;++i) {
			fCentroid1[i] = (float)(dAccumulator1[i]/((float)iClusterElements1));
			fCentroid2[i] = (float)(dAccumulator2[i]/((float)iClusterElements2));
		}		
		++iIterations;
		//printf("iteration: %d distortion: %f (%d %d)\n",iIterations,fDistortion1+fDistortion2,iClusterElements1,iClusterElements2);
		fDistortionCurrent = fDistortion1+fDistortion2;
	} while(fDistortionPrevious != fDistortionCurrent);
	
	//it seems like there is an implementation bug regarding the minimum cluster size, under k-means clusters are always going to be of at
	//least one element (the centroid), however that is not enough, fix that here and also in the REgressionTree class. Also by imposing
	
	//another bug: regresion tree should have a parameter which is the minimum number of observed gaussians to compute a transform
	
	//printf("iteration: %d distortion: %f (%d %d)\n",iIterations,fDistortion1+fDistortion2,iClusterElements1,iClusterElements2);
	
	if ((((float)iClusterElements1)/((float)iClusterElements2) > 5.0) || (((float)iClusterElements1)/((float)iClusterElements2) < 0.2)) {
		if (iAttempts > 0) {
			--iAttempts;
			//printf("-> %d\n",iAttempts);
			//goto start;
		}
	}
	
	// create the new clusters
	cluster->left = new GaussianCluster;
	cluster->left->fCentroid = new float[DIMENSIONALITY];
	cluster->left->fDistortion = fDistortion1;
	cluster->left->iDepth = cluster->iDepth+1;
	cluster->left->left = NULL; 
	cluster->left->right = NULL; 
	cluster->right = new GaussianCluster();
	cluster->right->fCentroid = new float[DIMENSIONALITY];
	cluster->right->fDistortion = fDistortion2;
	cluster->right->iDepth = cluster->iDepth+1;
	cluster->right->left = NULL; 
	cluster->right->right = NULL; 
	// store the centroids
	for(int i=0;i<DIMENSIONALITY;++i) {
		cluster->left->fCentroid[i] = fCentroid1[i];
		cluster->right->fCentroid[i] = fCentroid2[i];
	}
	// move the gaussians to the corresponding cluster
	for(VGaussianState::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		if ((*it)->iCluster == 0) {
			cluster->left->vGaussian.push_back(*it);	
		} else {
			cluster->right->vGaussian.push_back(*it);
		}
	}
	cluster->vGaussian.clear();
	
	delete [] fWeight;
}


// build the Vector Quantization tree
// at each step of the tree building process a node is selected for splitting, the split is done using k-means (k=2)
// the node selected for splitting is the one that contributes the most to the overall distortion
bool GSManager::buildVQTree(int iPrototypes) {

	//LGaussianCluster lGaussianCluster;	// keeps the tree-leaves (initially the root)
	float fDistortionTotal = 0.0;
	
	// starting time
	double dTimeStart = TimeUtils::getTimeMilliseconds();
	
	// initialize the random number generator
#ifdef __linux__
	srand(time(NULL));
#elif _WIN32
	srand(GetTickCount());
#endif

	// (1) create a sinle root cluster containing all the gaussians
	m_clusterRoot = new GaussianCluster;	
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		int iGaussians = -1;
		GaussianDecoding *gaussians = m_hmmStatesDecoding[i].getGaussians(iGaussians);
		for(int j=0 ; j < iGaussians ; ++j) {
			GaussianState *gaussianState = new GaussianState;
			gaussianState->gaussian = &gaussians[j];
			gaussianState->hmmStateDecoding = &m_hmmStatesDecoding[i];
			m_clusterRoot->vGaussian.push_back(gaussianState);	
		}
	}	
	m_clusterRoot->fCentroid = NULL;
	m_clusterRoot->fDistortion = 0.0;
	m_clusterRoot->iDepth = 0;
	m_clusterRoot->left = NULL;
	m_clusterRoot->right = NULL;
	m_lGaussianCluster.push_back(m_clusterRoot);
	
	// split a cluster at each iteration until the desired number of clusters is reached
	while(1) {
	
		// (2) select the cluster that contributes more to the overall distortion
		GaussianCluster *clusterToSplit = NULL;
		float fDistortionHighest = -FLT_MAX;
		LGaussianCluster::iterator jt;
		for(LGaussianCluster::iterator it = m_lGaussianCluster.begin() ; it != m_lGaussianCluster.end() ; ++it) {
			if ((*it)->fDistortion > fDistortionHighest) {
				fDistortionHighest = (*it)->fDistortion;
				clusterToSplit = (*it);
				jt = it;
			}
		}	
		assert(clusterToSplit != NULL);
		
		// (3) split the cluster using k-means (k=2) 
		kMeans(clusterToSplit);
		
		// remove the split cluster and insert the cluster resulting from the split
		m_lGaussianCluster.erase(jt);
		m_lGaussianCluster.push_back(clusterToSplit->left);
		m_lGaussianCluster.push_back(clusterToSplit->right);
		fDistortionTotal -= clusterToSplit->fDistortion;
		fDistortionTotal += clusterToSplit->left->fDistortion;
		fDistortionTotal += clusterToSplit->right->fDistortion;
		cout << "leaves: " << m_lGaussianCluster.size() << " distortion: " << fDistortionTotal << endl;
		if ((int)m_lGaussianCluster.size() == iPrototypes) {
			break;
		}
	}
	
	// ending time
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dMillisecondsInterval = dTimeEnd - dTimeStart;	
	
	// sanity check
	unsigned int iTreeLeaves = countLeaves(m_clusterRoot);
	assert(iTreeLeaves == m_lGaussianCluster.size());
	
	// show tree statistics
	float fDepthAverage = 0.0;
	int iDepthMinimum = INT_MAX;	
	int iDepthMaximum = 0;
	float fClusterSizeAverage = 0.0;
	unsigned int iClusterSizeMinimum = INT_MAX;
	unsigned int iClusterSizeMaximum = 0;
	for(LGaussianCluster::iterator it = m_lGaussianCluster.begin() ; it != m_lGaussianCluster.end() ; ++it) {
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
	fDepthAverage /= (float)iPrototypes;
	fClusterSizeAverage /= (float)iPrototypes;
	printf("Average depth: %.2f (min= %d max= %d)\n",fDepthAverage,iDepthMinimum,iDepthMaximum);
	printf("Average cluster size:  %.2f (min= %d max= %d)\n",fClusterSizeAverage,iClusterSizeMinimum,iClusterSizeMaximum); 
   printf("Building time: %.2f\n",dMillisecondsInterval/1000.0);	
   
   printClusterOccupationHistogram();

	return true;
}

// destroy the Vector Quantization tree (recursively)
void GSManager::destroyVQTree(GaussianCluster *cluster) {

	delete [] cluster->fCentroid;
	for(VGaussianState::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		delete *it;
	}

	if (cluster->left == NULL) {
		assert(cluster->right == NULL);	
	} else {
		destroyVQTree(cluster->left);
		destroyVQTree(cluster->right);
	}	
	delete cluster;
}

// counts the number of leaves in the tree (recursively)
int GSManager::countLeaves(GaussianCluster *cluster) {

	if (cluster->left == NULL) {
		assert(cluster->right == NULL);
		return 1;
	} else {
		return countLeaves(cluster->left)+countLeaves(cluster->right);
	}
}

// compute the weight vector for the computation of the weighted euclidean distance
float *GSManager::computeWeight(GaussianCluster *cluster) {

	// (1) compute the average of the covariance of all the gaussians
	double *dCovAverage = new double[DIMENSIONALITY];
	for(int i=0;i<DIMENSIONALITY;++i) {
		dCovAverage[i] = 0.0;
	}
	for(VGaussianState::iterator it = cluster->vGaussian.begin() ; it != cluster->vGaussian.end() ; ++it) {
		for(int i=0;i<DIMENSIONALITY;++i) {
			dCovAverage[i] += (*it)->gaussian->fCovariance[i]; 
		}
	}
	for(int i=0;i<DIMENSIONALITY;++i) {
		dCovAverage[i] /= ((float)cluster->vGaussian.size());
	}
	
	// (2) compute the inverse square root
	float *fWeight = new float[DIMENSIONALITY];
	for(int i=0;i<DIMENSIONALITY;++i) {
		fWeight[i] = (float)(1.0f/sqrt(dCovAverage[i]));
	}	
	
	delete [] dCovAverage;
	
	return fWeight;
}

// quantize a feature vector by finding the closest prototype (it performs a tree traversal)
int GSManager::quantize(float *fVector) {

	assert(m_clusterRoot != NULL);

	return 0;
}

// print the cluster occupation histogram
void GSManager::printClusterOccupationHistogram() {

	map<unsigned int,unsigned int> mClusterSizeCount;
	unsigned int iGaussiansClustered = 0;
	
	// create the histogram
	for(LGaussianCluster::iterator it = m_lGaussianCluster.begin() ; it != m_lGaussianCluster.end() ; ++it) {
		map<unsigned int,unsigned int>::iterator jt = mClusterSizeCount.find((*it)->vGaussian.size());
		if (jt != mClusterSizeCount.end()) {
			jt->second++;
		} else {
			mClusterSizeCount.insert(map<unsigned int,unsigned int>::value_type((*it)->vGaussian.size(),1));
		}
	}	
	
	// print the histogram
	for(map<unsigned int,unsigned int>::iterator it = mClusterSizeCount.begin() ; it != mClusterSizeCount.end() ; ++it) {
		printf("%6u %6u\n",it->first,it->second);	
		iGaussiansClustered += it->first*it->second;
	}
	printf("# gaussians clustered: %u\n",iGaussiansClustered);
}

};	// end-of-namespace



