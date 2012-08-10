/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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

#include "MLAccumulator.h"

// constructor
MLAccumulator::MLAccumulator(const char *strFilePhoneSet, const char *strFileConfigurationFeaturesAlignment, 
			const char *strFolderFeaturesAlignment, const char *strFileModelsAlignment, 
			unsigned char iAccumulatorType, const char *strFileOptionalSymbols, bool bMultiplePronunciations,
			unsigned char iContextModelingOrderAccumulatorsWW, unsigned char iContextModelingOrderAccumulatorsCW, 
			const char *strFileConfigurationFeaturesAcc, const char *strFolderFeaturesAcc, int iCovarianceModellingAcc, const char *strFileLexicon, const char *strFileMLF, const char *strFileAccumulators, const char *strFileGlobalDistribution, float fForwardPruningBeam, float fBackwardPruningBeam, int iTrellisMaxSize, bool bTrellisCache, int iTrellisCacheMaxSize, const char *strFileLog, int iLogVerbosity)
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
	m_strFileGlobalDistribution = strFileGlobalDistribution;
	m_fForwardPruningBeam = fForwardPruningBeam;
	m_fBackwardPruningBeam = fBackwardPruningBeam; 
	m_iTrellisMaxSize = iTrellisMaxSize;
	m_bTrellisCache = bTrellisCache;
	m_iTrellisCacheMaxSize = iTrellisCacheMaxSize;
	m_strFileLog = strFileLog;
	m_iLogVerbosity = iLogVerbosity;
	
	// single feature stream?
	if ((m_strFolderFeaturesAcc == NULL) || (strcmp(m_strFolderFeaturesAlignment,m_strFolderFeaturesAcc) == 0)) {
		m_bSingleFeatureStream = true;
	} else {
		m_bSingleFeatureStream = false;
	}	
	
	m_hmmManagerAlignment = NULL;
	m_hmmManagerAccumulation = NULL;
		
	return;
}

// destructor
MLAccumulator::~MLAccumulator()
{
	delete m_log;
	delete m_phoneSet;
	delete m_lexiconManager;
	/*if (m_hmmManagerAlignment != NULL) {
		delete m_hmmManagerAlignment;
		m_hmmManagerAlignment = NULL;
	}
	if (m_hmmManagerAccumulation != NULL) {
		delete m_hmmManagerAccumulation;
		m_hmmManagerAccumulation = NULL;
	}*/
	delete m_mlfFile;
	delete m_forwardBackwardX;
}


// initialize the accumulation
bool MLAccumulator::initialize() {

   // create the log
   if (m_strFileLog != NULL) {
   	m_log = new Log(m_strFileLog,m_iLogVerbosity);
   } else {
   	m_log = new Log();
   }

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
   m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet,m_log); 
   m_lexiconManager->load();
	
	// load the Master Label File
	m_mlfFile = new MLFFile(m_lexiconManager,m_strFileMLF,m_log,MODE_READ);	
	m_mlfFile->load();
	
	// load the HMMs used for the alignment
	m_hmmManagerAlignment = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
	m_hmmManagerAlignment->load(m_strFileModelsAlignment);
	if (m_iCovarianceModellingAcc == COVARIANCE_MODELLING_TYPE_DEFAULT) {
		m_iCovarianceModellingAcc = m_hmmManagerAlignment->getCovarianceModelling();
	}
	
	if (m_hmmManagerAlignment->initializeEstimation(m_iAccumulatorType,m_iContextModelingOrderAccumulatorsWW,
		m_iContextModelingOrderAccumulatorsCW) == false) {
		return false;
	}
	
	// load the HMMs used for the accumulation, one of these conditions should be met:
	// - double feature stream
	// - different convariance modeling
	if ((m_bSingleFeatureStream == false) || 
		(m_hmmManagerAlignment->getCovarianceModelling() != m_iCovarianceModellingAcc)) {
		m_hmmManagerAccumulation = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
		if (m_hmmManagerAccumulation->initializeModels(m_hmmManagerAlignment,
			m_iFeatureDimensionalityAcc,m_iCovarianceModellingAcc) == false) {
			return false;
		}
		if (m_hmmManagerAccumulation->initializeEstimation(m_iAccumulatorType,
			m_iContextModelingOrderAccumulatorsWW,m_iContextModelingOrderAccumulatorsCW ) == false) {
			return false;
		}
	}
	else {
		m_hmmManagerAccumulation = m_hmmManagerAlignment;
	}	
	
	// create the Forward-Backward object
	m_forwardBackwardX = new ForwardBackwardX(m_phoneSet,m_lexiconManager,m_hmmManagerAlignment,m_hmmManagerAccumulation,
		m_fForwardPruningBeam,m_fBackwardPruningBeam,m_iTrellisMaxSize,m_bTrellisCache,m_iTrellisCacheMaxSize,m_log);
	
	// transcription properties
	if (m_strFileOptionalSymbols != NULL) {
		m_vLexUnitOptional.push_back(m_lexiconManager->getLexUnitSilence());
	}
		
	// global distribution of the data
	//initializeGlobalDistribution();	

	return true;
}
		
// accumulate statistics
bool MLAccumulator::accumulate() {

	// make sure the HMMs are already initialized
	if (m_hmmManagerAlignment->areInitialized() == false) {
		return false;
	}
	if (m_hmmManagerAccumulation->areInitialized() == false) {
		return false;
	}
	
	char strFileFeatures[MAX_PATH_LENGTH+1];
	char strFileAccumulators[MAX_PATH_LENGTH];
	char strAux[MAX_PATH_LENGTH];
	char strMessage[MAX_LOG_MESSAGE_LENGTH+1];	
	
	double dLikelihoodTotal = 0.0;
	long iFeatureVectorsTotal = 0;
	long iFeatureVectorsUsedTotal = 0;
		
   double dBegin = Time::getTimeMilliseconds();
		
	// empty the accumulators
	m_hmmManagerAccumulation->resetAccumulators();
	
	// reset the global distribution of the data	
	//resetGlobalDistribution();
	
	// precompute constants used to speed-up emission probability computation
	if (m_hmmManagerAlignment->precomputeConstants() == false) {
		m_log->logError("unable to precompute the Gaussian estimation constants");	
		return false;
	}

	// (2) process each utterance in the MLF file 
	int iUtterance = 0;
	VMLFUtterance *vMLFUtterance = m_mlfFile->getUtterances();
	// at this point we might not know the total amount of audio but we do know the total number of utterances
	int iUtterancesTotal = vMLFUtterance->size();
	float fPercentageDisplayed = 0.0;
	for(VMLFUtterance::iterator it = vMLFUtterance->begin() ; it != vMLFUtterance->end() ; ++it, ++iUtterance) {
	
		// load the features for the estimation
		sprintf(strFileFeatures,"%s%c%s",m_strFolderFeaturesAlignment,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
		FeatureFile *featureFileAlignment = new FeatureFile(strFileFeatures,MODE_READ,FORMAT_FEATURES_FILE_DEFAULT,m_iFeatureDimensionalityAlignment);
		featureFileAlignment->load();
		int iFeatureVectorsAlignment = 0;
		float *fFeaturesAlignment = (float*)featureFileAlignment->getFeatureVectors(&iFeatureVectorsAlignment);
		delete featureFileAlignment;
		
		// load the features for the accumulation (if necessary)
		int iFeatureVectorsAcc = -1;
		float *fFeaturesAcc = NULL;
		if (m_bSingleFeatureStream == false) {
			sprintf(strFileFeatures,"%s%c%s",m_strFolderFeaturesAcc,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
			FeatureFile *featureFileAcc = new FeatureFile(strFileFeatures,MODE_READ,FORMAT_FEATURES_FILE_DEFAULT,m_iFeatureDimensionalityAcc);
			featureFileAcc->load();
			fFeaturesAcc = (float*)featureFileAcc->getFeatureVectors(&iFeatureVectorsAcc);
			delete featureFileAcc;
		} 
		// same feature stream for estimation and parameter update (default behaviour)
		else {
			iFeatureVectorsAcc = iFeatureVectorsAlignment;
			fFeaturesAcc = fFeaturesAlignment;
		}
		
		iFeatureVectorsTotal += iFeatureVectorsAlignment;	
		
		// process the utterance using Forward-Backward (get the occupation counts)
		double dUtteranceLikelihood;
		int iReturnValue = -1;
		Alignment *alignment = m_forwardBackwardX->processUtterance((*it)->vLexUnit,m_bMultiplePronunciations,
			m_vLexUnitOptional,fFeaturesAlignment,fFeaturesAcc,iFeatureVectorsAlignment,
			&dUtteranceLikelihood,iReturnValue);
		if (iReturnValue == UTTERANCE_PROCESSED_SUCCESSFULLY) {
			// count the audio actually used for training
			dLikelihoodTotal += dUtteranceLikelihood;
			iFeatureVectorsUsedTotal += iFeatureVectorsAlignment;
			// accumulate the statistics to compute the global distribution of the data
			//accumulateGlobalDistribution(fFeaturesAcc,iFeatureVectorsAcc);
		}
		// utterance discarded: show a message
		else {
			sprintf(strMessage,"unable to process utterance: \"%s\", reason: %s",strFileFeatures,
			m_forwardBackwardX->getErrorMessage(iReturnValue));
			m_log->logInformation(strMessage);
		}
		delete alignment;
		
		// clean-up
		if (fFeaturesAcc != fFeaturesAlignment) {
			delete [] fFeaturesAcc;
		}
		delete [] fFeaturesAlignment;
	
		if (iFeatureVectorsTotal >= 100*60*30) {
			//break;
		}
		
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
	double dEnd = Time::getTimeMilliseconds();
	double dMillisecondsInterval = dEnd - dBegin;
	
	// compute the Real Time Factor of the reestimation process
	float fRTF = (dMillisecondsInterval/10.0)/((float)iFeatureVectorsTotal);
	
	int iGaussians = m_hmmManagerAlignment->getNumberGaussianComponents();
	float fLikelihoodFrame = dLikelihoodTotal/((float)iFeatureVectorsUsedTotal);
	// compute audio available for training
	int iHours,iMinutes,iSeconds;
	Time::convertHundredths((double)iFeatureVectorsTotal,iHours,iMinutes,iSeconds);
	// compute audio used
	int iHoursUsed,iMinutesUsed,iSecondsUsed;
	Time::convertHundredths((double)iFeatureVectorsUsedTotal,iHoursUsed,iMinutesUsed,iSecondsUsed);
	
	// show the accumulation information
	printf(" likelihood= %.4f (%.2f) [%8d Gauss][RTF=%.4f][%d:%02d'%02d''][%d:%02d'%02d'']\n",dLikelihoodTotal,fLikelihoodFrame,iGaussians,
		fRTF,iHours,iMinutes,iSeconds,iHoursUsed,iMinutesUsed,iSecondsUsed);	
	
	// (8) dump the accumulators
	m_hmmManagerAccumulation->dumpAccumulators(m_strFileAccumulators);

	return true;
}



