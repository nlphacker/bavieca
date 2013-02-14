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


#include "BatchFile.h"
#include "HMMInitializer.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
HMMInitializer::HMMInitializer(int iDim, int iCovarianceModelling, PhoneSet *phoneSet)
{
	m_iDim = iDim;
	m_iCovarianceModelling = iCovarianceModelling;
	m_phoneSet = phoneSet;
	m_iState = STATE_MODELS_UNINITIALIZED;
}

// destructor
HMMInitializer::~HMMInitializer()
{
	
}

// compute the global mean and covariance
void HMMInitializer::computeGlobalDistribution(VMLFUtterance *vMLFUtterance, const char *strFolderFeatures) {

	// (1) load the MLF and compute the global distribution
	Vector<double> vObservation(m_iDim);
	SMatrix<double> mObservationSquare(m_iDim);	
	long lVectors = 0;
	
	cout << "loading MLF: ";
	fflush(stdout);
	
	for(VMLFUtterance::const_iterator it = vMLFUtterance->begin() ; it != vMLFUtterance->end() ; ++it) {
	
		// read the feature vectors from the file
		ostringstream ossFileFeatures;
		ossFileFeatures << strFolderFeatures << PATH_SEPARATOR << (*it)->strFilePattern;
		FeatureFile featureFile(ossFileFeatures.str().c_str(),MODE_READ);
		featureFile.load();
		int iFeatures;
		float *fFeatures = featureFile.getFeatureVectors(&iFeatures);
		
		// update counters
		for(int i = 0 ; i < iFeatures ; ++i) {
			
			VectorStatic<float> vFeatures(fFeatures+(i*m_iDim),m_iDim);
			vObservation.add(1.0,vFeatures);
			mObservationSquare.addSquare(1.0,vFeatures);
		}	
		lVectors += iFeatures;
		
		delete [] fFeatures;
	}
	if (lVectors == 0) {
		BVC_ERROR << "no feature vectors found to compute the global distribution";
	}
	
	// show the final count
	int iHours = lVectors/(100*60*60);
	int iMinutes = (lVectors%(100*60*60))/(100*60);
	printf("%dh:%02d' of speech loaded\n",iHours,iMinutes);
	
	// (2) compute the global mean and variance
	m_vMeanGlobal = new Vector<float>(m_iDim);
	m_vMeanGlobal->mul(1.0f/((float)lVectors),vObservation);
	
	m_mCovarianceGlobal = new SMatrix<float>(m_iDim);
	m_mCovarianceGlobal->add(1.0f/((float)lVectors),mObservationSquare);
	m_mCovarianceGlobal->addSquare(-1.0f,*m_vMeanGlobal);	
}

// initialize the HMMs to have the mean and covariance of the global set of features
HMMManager *HMMInitializer::initializeModelsFlatStart(VMLFUtterance *vMLFUtterance, const char *strFolderFeatures) {

	// check state
	assert(m_iState == STATE_MODELS_UNINITIALIZED);
	
	// (1) load the MLF and compute the global distribution
	computeGlobalDistribution(vMLFUtterance,strFolderFeatures);
	
	// (3) initialize the acoustic models 
	HMMManager *hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
	
	// create the HMM prototype for monophones
	int iComponents = 1;
	if (hmmManager->createSingleGaussianMonophoneModelsPrototype(m_iDim,
		m_iCovarianceModelling,iComponents) == false) {
		return NULL;
	}

	hmmManager->initializeModels(*m_vMeanGlobal,*m_mCovarianceGlobal);	
	hmmManager->setInitialized(true);
	
	m_iState = STATE_MODELS_INITIALIZED;

	return hmmManager;
}

// initialize the HMMs from baseline HMMs
HMMManager *HMMInitializer::initializeModelsBootstrap(VMLFUtterance *vMLFUtterance, const char *strFolderFeatures, const char *strBootstrapModelsFile, bool bFloorCovariance, float fCovarianceFloor, float *fMeanGlobal, float *fCovarianceGlobal) {

	// check state
	if (m_iState != STATE_MODELS_UNINITIALIZED) {
		return NULL;
	}
	
	// (1) load the MLF and compute the global distribution
	computeGlobalDistribution(vMLFUtterance,strFolderFeatures);
	
	HMMManager *hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);	
	
	// load the HMMs from a file
	hmmManager->load(strBootstrapModelsFile);
	
	m_iState = STATE_MODELS_INITIALIZED;

	return hmmManager;
}

};	// end-of-namespace

