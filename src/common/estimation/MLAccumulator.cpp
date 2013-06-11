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

#include <stdexcept>

#include "ConfigurationFeatures.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "ForwardBackwardX.h"
#include "HMMManager.h"
#include "LexUnitsFile.h"
#include "MLAccumulator.h"
#include "MLFFile.h"
#include "PhoneSet.h"
#include "PhoneticRulesManager.h"
#include "TimeUtils.h"

namespace Bavieca {

// constructor
MLAccumulator::MLAccumulator(const char *strFilePhoneSet, const char *strFileConfigurationFeaturesAlignment, 
			const char *strFolderFeaturesAlignment, const char *strFileModelsAlignment, 
			unsigned char iAccumulatorType, const char *strFileOptionalSymbols, bool bMultiplePronunciations,
			unsigned char iContextModelingOrderAccumulatorsWW, unsigned char iContextModelingOrderAccumulatorsCW, 
			const char *strFileConfigurationFeaturesAcc, const char *strFolderFeaturesAcc, int iCovarianceModellingAcc, const char *strFileLexicon, const char *strFileMLF, const char *strFileAccumulators, float fForwardPruningBeam, float fBackwardPruningBeam, int iTrellisMaxSize, bool bTrellisCache, int iTrellisCacheMaxSize)
{
	m_strFilePhoneSet = strFilePhoneSet;
	m_strFileConfigurationFeaturesAlignment = strFileConfigurationFeaturesAlignment; 
	m_strFolderFeaturesAlignment = strFolderFeaturesAlignment;
	m_strFileModelsAlignment = strFileModelsAlignment;
	m_iAccumulatorType = iAccumulatorType;
	m_iContextModelingOrderAccumulatorsWW = iContextModelingOrderAccumulatorsWW;
	m_iContextModelingOrderAccumulatorsCW = iContextModelingOrderAccumulatorsCW;
	m_strFileOptionalSymbols = strFileOptionalSymbols;
	m_bMultiplePronunciations = bMultiplePronunciations;
	m_strFileConfigurationFeaturesAcc = strFileConfigurationFeaturesAcc;
	m_strFolderFeaturesAcc = strFolderFeaturesAcc;
	m_iCovarianceModellingAcc = iCovarianceModellingAcc;
	m_strFileLexicon = strFileLexicon;
	m_strFileMLF = strFileMLF;
	m_strFileAccumulators = strFileAccumulators;
	m_fForwardPruningBeam = fForwardPruningBeam;
	m_fBackwardPruningBeam = fBackwardPruningBeam; 
	m_iTrellisMaxSize = iTrellisMaxSize;
	m_bTrellisCache = bTrellisCache;
	m_iTrellisCacheMaxSize = iTrellisCacheMaxSize;
	
	// single feature stream?
	if ((m_strFolderFeaturesAcc == NULL) || (strcmp(m_strFolderFeaturesAlignment,m_strFolderFeaturesAcc) == 0)) {
		m_bSingleFeatureStream = true;
	} else {
		m_bSingleFeatureStream = false;
	}	
	
	m_hmmManagerAlignment = NULL;
	m_hmmManagerAccumulation = NULL;
}

// destructor
MLAccumulator::~MLAccumulator()
{
	delete m_phoneSet;
	delete m_lexiconManager;
	delete m_mlfFile;
	delete m_hmmManagerAlignment;
	if (m_hmmManagerAccumulation != m_hmmManagerAlignment) {
		delete m_hmmManagerAccumulation;
	}
	delete m_forwardBackwardX;
}


// initialize the accumulation
void MLAccumulator::initialize() {

   // load the phone set
   m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();

   // load the feature configuration (alignment)
   ConfigurationFeatures *configurationFeaturesAlignment = new ConfigurationFeatures(m_strFileConfigurationFeaturesAlignment);
   configurationFeaturesAlignment->load();
   m_iFeatureDimensionalityAlignment = configurationFeaturesAlignment->getDimensionality();
   delete configurationFeaturesAlignment;
   
	// load the feature configuration (accumulation)
	if (m_bSingleFeatureStream == false) {
		assert(m_strFileConfigurationFeaturesAcc != NULL);
		ConfigurationFeatures *configurationFeaturesAcc = new ConfigurationFeatures(m_strFileConfigurationFeaturesAcc);
		configurationFeaturesAcc->load();
		m_iFeatureDimensionalityAcc = configurationFeaturesAcc->getDimensionality();
		delete configurationFeaturesAcc;
	}
	else {
		m_iFeatureDimensionalityAcc = m_iFeatureDimensionalityAlignment;
	}

   // load the lexicon
   m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
   m_lexiconManager->load();
	
	// load the Master Label File
	m_mlfFile = new MLFFile(m_lexiconManager,m_strFileMLF,MODE_READ);	
	m_mlfFile->load();
	
	// load the HMMs used for the alignment
	m_hmmManagerAlignment = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
	m_hmmManagerAlignment->load(m_strFileModelsAlignment);
	if (m_iCovarianceModellingAcc == COVARIANCE_MODELLING_TYPE_DEFAULT) {
		m_iCovarianceModellingAcc = m_hmmManagerAlignment->getCovarianceModelling();
	}
	
	m_hmmManagerAlignment->initializeEstimation(m_iAccumulatorType,m_iContextModelingOrderAccumulatorsWW,
		m_iContextModelingOrderAccumulatorsCW);
	
	// load the HMMs used for the accumulation, one of these conditions should be met:
	// - double feature stream
	// - different convariance modeling
	if ((m_bSingleFeatureStream == false) || 
		(m_hmmManagerAlignment->getCovarianceModelling() != m_iCovarianceModellingAcc)) {
		m_hmmManagerAccumulation = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
		m_hmmManagerAccumulation->initializeModels(m_hmmManagerAlignment,
			m_iFeatureDimensionalityAcc,m_iCovarianceModellingAcc);
		m_hmmManagerAccumulation->initializeEstimation(m_iAccumulatorType,
			m_iContextModelingOrderAccumulatorsWW,m_iContextModelingOrderAccumulatorsCW );
	}
	else {
		m_hmmManagerAccumulation = m_hmmManagerAlignment;
	}	
	
	// create the Forward-Backward object
	m_forwardBackwardX = new ForwardBackwardX(m_phoneSet,m_lexiconManager,m_hmmManagerAlignment,m_hmmManagerAccumulation,
		m_fForwardPruningBeam,m_fBackwardPruningBeam,m_iTrellisMaxSize,m_bTrellisCache,m_iTrellisCacheMaxSize);
	
	// transcription properties
	if (m_strFileOptionalSymbols) {
		// get the optional symbols
		LexUnitsFile lexUnitsFile(m_lexiconManager,m_strFileOptionalSymbols);
		lexUnitsFile.load();
		lexUnitsFile.getLexUnits(m_vLexUnitOptional);
	} else {
		m_vLexUnitOptional.push_back(m_lexiconManager->getLexUnitSilence());
	}
}
		
// accumulate statistics
void MLAccumulator::accumulate() {

	// make sure the HMMs are already initialized
	assert(m_hmmManagerAlignment->areInitialized());
	assert(m_hmmManagerAccumulation->areInitialized());
	
	double dLikelihoodTotal = 0.0;
	long iFeatureVectorsTotal = 0;
	long iFeatureVectorsUsedTotal = 0;
		
	double dBegin = TimeUtils::getTimeMilliseconds();
		
	// empty the accumulators
	m_hmmManagerAccumulation->resetAccumulators();
	
	// precompute constants used to speed-up emission probability computation
	m_hmmManagerAlignment->precomputeConstants();

	// (2) process each utterance in the MLF file 
	int iUtterance = 0;
	VMLFUtterance *vMLFUtterance = m_mlfFile->getUtterances();
	// at this point we might not know the total amount of audio but we do know the total number of utterances
	unsigned int iUtterancesTotal = (unsigned int)vMLFUtterance->size();
	float fPercentageDisplayed = 0.0;
	for(VMLFUtterance::iterator it = vMLFUtterance->begin() ; it != vMLFUtterance->end() ; ++it, ++iUtterance) {
	
		// load the features for the estimation
		ostringstream strFileFeatures;
		strFileFeatures << m_strFolderFeaturesAlignment << PATH_SEPARATOR << (*it)->strFilePattern;
		FeatureFile featureFileAlignment(strFileFeatures.str().c_str(),MODE_READ,FORMAT_FEATURES_FILE_DEFAULT,
			m_iFeatureDimensionalityAlignment);
		try {
			featureFileAlignment.load();
		} catch (std::runtime_error &e) {
			std::cerr << e.what() << std::endl;
			BVC_WARNING << "unable to load the features file: " << strFileFeatures.str();
			continue;
		}	
			
		Matrix<float> *mFeaturesAlignment = featureFileAlignment.getFeatureVectors();
		
		// load the features for the accumulation (if necessary)
		Matrix<float> *mFeaturesAcc = mFeaturesAlignment;
		if (m_bSingleFeatureStream == false) {
			ostringstream strFileFeatures;
			strFileFeatures << m_strFolderFeaturesAcc << PATH_SEPARATOR << (*it)->strFilePattern;
			FeatureFile featureFileAcc(strFileFeatures.str().c_str(),MODE_READ,FORMAT_FEATURES_FILE_DEFAULT,
				m_iFeatureDimensionalityAcc);
			try {
				featureFileAcc.load();
			} catch (std::runtime_error &e) {
				std::cerr << e.what() << std::endl;
				BVC_WARNING << "unable to load the features file: " << strFileFeatures.str();
				delete mFeaturesAlignment;
				continue;
			}	
			mFeaturesAcc = featureFileAcc.getFeatureVectors();
		} 
		
		iFeatureVectorsTotal += mFeaturesAlignment->getRows();	
		
		// process the utterance using Forward-Backward (get the occupation counts)
		double dUtteranceLikelihood;
		const char *strReturnCode = NULL;
		Alignment *alignment = m_forwardBackwardX->processUtterance((*it)->vLexUnit,m_bMultiplePronunciations,
			m_vLexUnitOptional,*mFeaturesAlignment,*mFeaturesAcc,&dUtteranceLikelihood,&strReturnCode);
		if (strcmp(strReturnCode,FB_RETURN_CODE_SUCCESS) ==0) {
			// count the audio actually used for training
			dLikelihoodTotal += dUtteranceLikelihood;
			iFeatureVectorsUsedTotal += mFeaturesAlignment->getRows();
		}
		// utterance discarded: show a message
		else {
			BVC_WARNING << "unable to process utterance: \"" << strFileFeatures.str() << "\", reason: " << strReturnCode;
		}
		delete alignment;
		
		// clean-up
		if (mFeaturesAcc != mFeaturesAlignment) {
			delete mFeaturesAcc;
		}
		delete mFeaturesAlignment;
		
		// update the progress bar if necessary
		float fPercentage = (((float)iUtterance)*100)/((float)iUtterancesTotal);
		if (fPercentage >= fPercentageDisplayed + 10.0) {
			fPercentageDisplayed += 10.0;
			printf("*");
			fflush(stdout);
		}
	}
	// update the progress bar if necessary
	while (fPercentageDisplayed < 100.0) {
		printf("*");
		fPercentageDisplayed += 10.0;
	}
	
	// get the iteration end time
	double dEnd = TimeUtils::getTimeMilliseconds();
	double dMillisecondsInterval = dEnd - dBegin;
	
	// compute the Real Time Factor of the reestimation process
	float fRTF = ((float)dMillisecondsInterval/10.0f)/((float)iFeatureVectorsTotal);
	
	int iGaussians = m_hmmManagerAlignment->getNumberGaussianComponents();
	float fLikelihoodFrame = ((float)dLikelihoodTotal)/((float)iFeatureVectorsUsedTotal);
	// compute audio available for training
	int iHours,iMinutes,iSeconds;
	TimeUtils::convertHundredths((double)iFeatureVectorsTotal,iHours,iMinutes,iSeconds);
	// compute audio used
	int iHoursUsed,iMinutesUsed,iSecondsUsed;
	TimeUtils::convertHundredths((double)iFeatureVectorsUsedTotal,iHoursUsed,iMinutesUsed,iSecondsUsed);
	
	// show the accumulation information
	printf(" likelihood= %.4f (%.2f) [%8d Gauss][RTF=%.4f][%d:%02d'%02d''][%d:%02d'%02d'']\n",dLikelihoodTotal,fLikelihoodFrame,iGaussians,
		fRTF,iHours,iMinutes,iSeconds,iHoursUsed,iMinutesUsed,iSecondsUsed);	
	
	// (8) dump the accumulators
	m_hmmManagerAccumulation->dumpAccumulators(m_strFileAccumulators);
}

};	// end-of-namespace


