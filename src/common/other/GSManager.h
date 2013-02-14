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


#ifndef GSMANAGER_H
#define GSMANAGER_H

#include "HMMManager.h"

using namespace std;

#include <deque>

namespace Bavieca {

// keeps a gaussian and a pointer to the state it comes from
typedef struct {
	GaussianDecoding *gaussian;				// Gaussian pointer
	HMMStateDecoding *hmmStateDecoding;		// HMM-state owner of the Gaussian	
	int iCluster;									// temporal cluster to wich the gaussian is assigned
} GaussianState;

typedef vector<GaussianState*> VGaussianState;

typedef struct _GaussianCluster {
	VGaussianState vGaussian;		// cluster gaussians 
	float *fCentroid;					// cluster centroid
	float fDistortion;				// cluster distortion contribution
	int iDepth;							// cluster depth in the tree
	_GaussianCluster *left;			// left children cluster
	_GaussianCluster *right;		// right children cluster
} GaussianCluster;

typedef list<GaussianCluster*> LGaussianCluster;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
// Gaussian Selection Manager
class GSManager {

	private:
	
		HMMStateDecoding *m_hmmStatesDecoding;			// HMM-states
		int m_iHMMStates;										// number of HMM-states	
		GaussianCluster *m_clusterRoot;					// root cluster
		LGaussianCluster m_lGaussianCluster;				// gaussian clusters
		
		// compute the weighted euclidean distance between two vectors
		inline float weightedEuclideanDistance(float *fV1, float *fV2, float *fW) {
		
			float fAcc = 0.0;
			float fAux = 0.0;
			for(int i=0 ; i<DIMENSIONALITY ; ++i) {
				fAux = (fW[i]*(fV1[i]-fV2[i]));
				fAcc += fAux*fAux;
			}
			fAcc /= ((float)DIMENSIONALITY);
		
			return fAcc;
		}
			
		// compute the euclidean distance between two vectors
		inline float euclideanDistance(float *fV1, float *fV2) {
		
			float fAux = 0.0;
			for(int i=0 ; i<DIMENSIONALITY ; ++i) {
				fAux += (fV1[i]-fV2[i])*(fV1[i]-fV2[i]);
			}
			
			return fAux;
		}			

		// destroy the Vector Quantization tree (recursively)
		void destroyVQTree(GaussianCluster *cluster);
	
	public:
		
		// constructor
		GSManager(HMMManager *hmmManager);

		// destructor
		~GSManager();
		
		// perform k-means on a cluster 
		// note: it returns the new clusters as children of the given cluster
		void kMeans(GaussianCluster *cluster);
		
		// build the Vector Quantization tree
		bool buildVQTree(int iPrototypes);	
		
		// destroy the Vector Quantization tree (recursively)
		inline void destroyVQTree() {
		
			destroyVQTree(m_clusterRoot);
		}
		
		// counts the number of leaves in the tree (recursively)
		int countLeaves(GaussianCluster *cluster);	
		
		// return a random number from a range
		inline int getRandomNumber(int iStart, int iEnd) {
		
			#ifdef __linux__	
			// on linux 32 bits RAND_MAX = 2147483647 (this is the maximum signed integer)
			assert(RAND_MAX > iEnd-iStart);			
			// this is a very weak random number
			return iStart + rand()%(iEnd-iStart+1);
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
		
		// compute the weight vector for the computation of the weighted euclidean distance
		float *computeWeight(GaussianCluster *cluster);
		
		// quantize a feature vector by finding the closest prototype (it performs a tree traversal)
		int quantize(float *fVector);
		
		// print the cluster occupation histogram
		void printClusterOccupationHistogram();

};

};	// end-of-namespace

#endif
