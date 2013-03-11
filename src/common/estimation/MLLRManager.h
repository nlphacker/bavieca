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


#ifndef MLLRMANAGER_H
#define MLLRMANAGER_H

using namespace std;

#include <map>

#include <iostream>
#include <fstream>

#include "HMMManager.h"
#include "Matrix.h"
#include "RegressionTree.h"

namespace Bavieca {

class Alignment;
class BestPath;
class Gaussian;
class LexiconManager;
class PhoneSet;

typedef struct {
	string strSpeakerId;
	float *fMeanTransform;
	float *fCovarianceTransform;	
} SpeakerMLLRData;

typedef map<string,SpeakerMLLRData*> MSpeakerMLLRData;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class MLLRManager {

	private:
		
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManager;
		string m_strFileRegressionTree;
		RegressionTree *m_regressionTree;
		VGaussianDecoding m_vGaussianWithOccupation;
		int m_iAdaptationFrames;
		int m_iBaseClasses;
		unsigned char m_iClusteringMethod;
		bool m_bMeanOnly;												// whether only the mean transformation will be applied
		int m_iDim;
		
		float m_fMinimumOccupationTransform;
		int m_iMinimumGaussianComponentsObserved;
		bool m_bBestComponentOnly;
		GaussianStats **m_gaussianStats;

	public:
		
		// constructor
		MLLRManager(PhoneSet *phoneSet, HMMManager *hmmManager, 
			const char *strFileRegressionTree, float fMinimumOccupationTransform, 
			int iMinimumGaussianComponentsObserved, bool bBestComponentOnly, bool bMeanOnly);

		// destructor
		~MLLRManager();
		
		// initialization
		void initialize();	
		
		// compute transforms from the adaptation data (the adaptation data is stored in the gaussian accumulator)
		// there is an array containing all the gaussians with associated adaptation data
		void computeTransforms();
		
		// update the HMM-state parameters using the computed tranforms
		void applyTransforms();	
		
		// feed adaptation data from an alignment
		void feedAdaptationData(float *fFeatures, unsigned int iFeatures, Alignment *alignment, double *dLikelihood);	
		
		// feed adaptation data from an alignment into the adaptation process
		void feedAdaptationData(float *fFeatures, unsigned int iFeatures, VPhoneAlignment *vPhoneAlignment, double *dLikelihood);
		
		// feed adaptation data from a batch file containing entries (rawFile alignmentFile)
		void feedAdaptationData(const char *strBatchFile, const char *strAlignmentFormat, double *dLikelihood, bool bVerbose);
			
		// store the transforms to the given file
		void storeTransforms(const char *strFile);
};

};	// end-of-namespace

#endif
