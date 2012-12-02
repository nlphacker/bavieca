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


#include "MLLRManager.h"
#include "AlignmentFile.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "BatchFile.h"

#include "Alignment.h"
#include "BestPath.h"
#include "HMMManager.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "RegressionTree.h"

namespace Bavieca {

// constructor
MLLRManager::MLLRManager(PhoneSet *phoneSet, HMMManager *hmmManager, 
	const char *strFileRegressionTree, float fMinimumOccupationTransform, 
	int iMinimumGaussianComponentsObserved, bool bBestComponentOnly, bool bMeanOnly)
{
	m_phoneSet = phoneSet;
	m_hmmManager = hmmManager;	
	m_iAdaptationFrames = 0;
	m_fMinimumOccupationTransform = fMinimumOccupationTransform;
	m_iMinimumGaussianComponentsObserved = iMinimumGaussianComponentsObserved;
	m_bBestComponentOnly = bBestComponentOnly;
	m_bMeanOnly = bMeanOnly;
	m_iDim = m_hmmManager->getFeatureDimensionality();
	
	m_regressionTree = NULL;
	m_strFileRegressionTree = strFileRegressionTree;
}

// destructor
MLLRManager::~MLLRManager()
{
	// destroy the regression tree
	if (m_regressionTree != NULL) {
		delete m_regressionTree;
	}
	if (m_gaussianStats) {
		for(int i=0 ; i < m_hmmManager->getNumberGaussianComponents() ; ++i) {
			if (m_gaussianStats[i]) {	
				delete m_gaussianStats[i]->vObservation;
				delete m_gaussianStats[i];
			}
		}
		delete [] m_gaussianStats;
	}
}

// initialization
void MLLRManager::initialize() {

	m_regressionTree = new RegressionTree(m_hmmManager);	

	// load the regression tree from disk
	m_regressionTree->load(m_strFileRegressionTree.c_str());
	
	// get the number of base-classes
	m_iBaseClasses = m_regressionTree->getBaseClasses();
	
	m_gaussianStats = new GaussianStats*[m_hmmManager->getNumberGaussianComponents()];
	for(int i=0 ; i < m_hmmManager->getNumberGaussianComponents() ; ++i) {
		m_gaussianStats[i] = NULL;
	}
}


// compute transforms from the adaptation data (the adaptation data is stored in the gaussian accumulator)
// there is an array containing all the gaussians with associated adaptation data
void MLLRManager::computeTransforms() {

	float fAdpatationSeconds = m_iAdaptationFrames/100.0f;

	cout << "# seconds of adaptation data: " << fAdpatationSeconds << endl;
	cout << "# Gaussian distributions containing adaptation data: " << m_vGaussianWithOccupation.size() << endl;
	
	// check that there is enough data to compute at least a single tranform
	if (m_fMinimumOccupationTransform > m_iAdaptationFrames) {
		BVC_ERROR << "not enough adaptation data to robustly estimate a transform (required " << m_fMinimumOccupationTransform << ", available " << m_iAdaptationFrames << ")";
	}

	// compute the transform for each base-class
	m_regressionTree->computeTransforms(m_fMinimumOccupationTransform,m_iMinimumGaussianComponentsObserved
		,m_gaussianStats,m_bMeanOnly);
}

// update the HMM-state parameters using the computed tranforms
void MLLRManager::applyTransforms() {

	m_regressionTree->applyTransforms(true,!m_bMeanOnly);
}

// feed adaptation data from an alignment
void MLLRManager::feedAdaptationData(float *fFeatures, int iFeatures, Alignment *alignment, double *dLikelihood) {

	// sanity check
	assert(iFeatures == alignment->getFrames());

	m_iAdaptationFrames += iFeatures;

	*dLikelihood = 0.0;
	for(int t=0 ; t<iFeatures ; ++t) {
		FrameAlignment *frameAlignment = alignment->getFrameAlignment(t);
		float *fFeatureVector = fFeatures+(t*m_iDim);
		for(FrameAlignment::iterator it = frameAlignment->begin() ; it != frameAlignment->end() ; ++it) {
			HMMStateDecoding *hmmStateDecoding = m_hmmManager->getHMMStateDecoding((*it)->iHMMState);
			// compute the contribution of each Gaussian 
			// (case 1) all the frame-level adaptation data goes to the best scoring Gaussian component (faster)
			if (m_bBestComponentOnly) {
				float fLikelihood = -FLT_MAX;
				GaussianDecoding *gaussian = hmmStateDecoding->getBestScoringGaussian(fFeatureVector,&fLikelihood);
				*dLikelihood += std::max<float>(fLikelihood,LOG_LIKELIHOOD_FLOOR);
				GaussianStats *gaussianStats = NULL;
				if (m_gaussianStats[gaussian->iId] == NULL) {
					gaussianStats = new GaussianStats;
					gaussianStats->gaussian = gaussian;
					gaussianStats->dOccupation = 0.0;
					gaussianStats->vObservation = new Vector<double>(m_iDim);
					gaussianStats->vObservation->zero();
					m_gaussianStats[gaussian->iId] = gaussianStats;
				} else {
					gaussianStats = m_gaussianStats[gaussian->iId];
				}
				gaussianStats->dOccupation += 1.0;
				gaussianStats->vObservation->add(1.0,VectorStatic<float>(fFeatureVector,m_iDim));
				m_regressionTree->accumulateStatistics(fFeatureVector,1.0,gaussianStats);	
			} 
			// (case 2) adaptation data is shared across all components (slightly more accurate)
			else {
				double *dProbGaussian = new double[hmmStateDecoding->getGaussianComponents()];
				double dProbTotal = 0.0;
				for(int g=0 ; g < hmmStateDecoding->getGaussianComponents() ; ++g) {
					dProbGaussian[g] = exp(hmmStateDecoding->computeGaussianProbability(g,fFeatureVector));
					dProbTotal += dProbGaussian[g];
					assert((finite(dProbGaussian[g])) && (finite(dProbTotal)));
				}
				*dLikelihood += max(log(dProbTotal),LOG_LIKELIHOOD_FLOOR);
				for(int g=0 ; g < hmmStateDecoding->getGaussianComponents() ; ++g) {
					dProbGaussian[g] /= dProbTotal;
					double dOccupation = dProbGaussian[g]*1.0;
					GaussianDecoding *gaussian = hmmStateDecoding->getGaussian(g);
					GaussianStats *gaussianStats = NULL;
					if (m_gaussianStats[gaussian->iId] == NULL) {
						gaussianStats = new GaussianStats;
						gaussianStats->gaussian = gaussian;
						gaussianStats->dOccupation = 0.0;
						gaussianStats->vObservation = new Vector<double>(m_iDim);
						gaussianStats->vObservation->zero();
						m_gaussianStats[gaussian->iId] = gaussianStats;
					} else {
						gaussianStats = m_gaussianStats[gaussian->iId];
					}
					gaussianStats->dOccupation += dOccupation;
					gaussianStats->vObservation->add(dOccupation,VectorStatic<float>(fFeatureVector,m_iDim));
					m_regressionTree->accumulateStatistics(fFeatureVector,dOccupation,gaussianStats);
				}
				delete [] dProbGaussian;
			}	
		}
	}
}

// feed adaptation data from a batch file containing entries (rawFile alignmentFile)
void MLLRManager::feedAdaptationData(const char *strBatchFile, const char *strAlignmentFormat, 
	double *dLikelihood, bool bVerbose) {

	BatchFile *batchFile = new BatchFile(strBatchFile,"features|alignment");
	batchFile->load();
	
	*dLikelihood = 0.0;
	
	//for(int i=0 ; i < batchFile->size() ; ++i) {	
	for(int i=0 ; i < 5 ; ++i) {	
		
		// load the alignment
		Alignment *alignment = NULL;
		// text format
		if (strcmp(strAlignmentFormat,"text") == 0) {	
			AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,NULL);
			VPhoneAlignment *vPhoneAlignment = alignmentFile->load(batchFile->getField(i,"alignment"));
			assert(vPhoneAlignment);
			alignment = AlignmentFile::toAlignment(m_phoneSet,m_hmmManager,vPhoneAlignment);
			AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
			delete alignmentFile;
		} 
		// binary format
		else {
			alignment = Alignment::load(batchFile->getField(i,"alignment"),NULL);
			assert(alignment);	
		}
		
		// load the feature vectors
		FeatureFile *featureFile = new FeatureFile(batchFile->getField(i,"features"),MODE_READ);
		featureFile->load();
		int iFeatureVectors = -1;
		float *fFeatures = featureFile->getFeatureVectors(&iFeatureVectors);
		
		// check consistency	
		if (iFeatureVectors != alignment->getFrames()) {
			BVC_ERROR << "inconsistent number of feature vectors / alignment file";
		}
		
		// accumulate adaptation data
		double dLikelihoodAlignment = 0.0;
		feedAdaptationData(fFeatures,iFeatureVectors,alignment,&dLikelihoodAlignment);
		if (bVerbose) {
			printf("loaded file: %s likelihood: %12.2f\n",batchFile->getField(i,"alignment"),dLikelihoodAlignment);
		}
		*dLikelihood += dLikelihoodAlignment;
		
		// clean-up
		delete alignment;
		delete featureFile;
		delete [] fFeatures;
	}
	
	if (bVerbose == true) {
		printf("total likelihood: %14.4f\n",*dLikelihood);
	}

	delete batchFile;
}

// store the transforms to the given file
void MLLRManager::storeTransforms(const char *strFile) {

	m_regressionTree->storeTransforms(strFile);
}

};	// end-of-namespace


