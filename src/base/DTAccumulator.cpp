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

#include "DTAccumulator.h"

// constructor
DTAccumulator::DTAccumulator(const char *strFilePhoneSet, const char *strFileConfigurationFeatures, 
			const char *strFolderFeatures, const char *strFileModels, const char *strFileOptionalSymbols, 
			bool bMultiplePronunciations, const char *strFileLexicon, const char *strFileMLF, 
			const char *strFolderLattices, float fScaleAM, float fScaleLM, const char *strFileAccumulatorsNum, 
			const char *strFileAccumulatorsDen, const char *strObjectiveFunction, float fBoostingFactor,
			bool bCanceledStatistics, float fForwardPruningBeam, float fBackwardPruningBeam, 
			int iTrellisMaxSize, bool bTrellisCache, int iTrellisCacheMaxSize, const char *strFileLog, 
			int iLogVerbosity)
{
	m_strFilePhoneSet = strFilePhoneSet;
	m_strFileConfigurationFeatures = strFileConfigurationFeatures; 
	m_strFolderFeatures = strFolderFeatures;
	m_strFileModels = strFileModels;
	m_strFileOptionalSymbols = strFileOptionalSymbols;
	m_bMultiplePronunciations = bMultiplePronunciations;
	m_strFileLexicon = strFileLexicon;
	m_strFileMLF = strFileMLF;
	m_strFolderLattices = strFolderLattices;
	m_fScaleAM = fScaleAM;
	m_fScaleLM = fScaleLM;
	m_strFileAccumulatorsNum = strFileAccumulatorsNum;
	m_strFileAccumulatorsDen = strFileAccumulatorsDen;
	m_strObjectiveFunction = strObjectiveFunction;
	m_fBoostingFactor = fBoostingFactor;
	m_bCanceledStatistics = bCanceledStatistics;
	m_fForwardPruningBeam = fForwardPruningBeam;
	m_fBackwardPruningBeam = fBackwardPruningBeam; 
	m_iTrellisMaxSize = iTrellisMaxSize;
	m_bTrellisCache = bTrellisCache;
	m_iTrellisCacheMaxSize = iTrellisCacheMaxSize;
	m_strFileLog = strFileLog;
	m_iLogVerbosity = iLogVerbosity;
	
	m_hmmManager = NULL;
		
	return;	
}

// destructor
DTAccumulator::~DTAccumulator()
{
	delete m_log;
	delete m_phoneSet;
	delete m_lexiconManager;
	delete m_mlfFile;
	delete m_forwardBackwardX;
	delete m_forwardBackward;
	delete m_hmmManager;
}


// initialize the accumulation
bool DTAccumulator::initialize() {

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
   ConfigurationFeatures *configurationFeatures = new ConfigurationFeatures(m_strFileConfigurationFeatures);
   configurationFeatures->load();
   m_iFeatureDimensionality = configurationFeatures->getDimensionality();
   delete configurationFeatures;

   // load the lexicon
   m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet,m_log); 
   m_lexiconManager->load();
	
	// load the Master Label File
	m_mlfFile = new MLFFile(m_lexiconManager,m_strFileMLF,m_log,MODE_READ);	
	m_mlfFile->load();
	
	// load the HMMs used for the alignment
	m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
	m_hmmManager->load(m_strFileModels);
	
	if (m_hmmManager->initializeEstimation(ACCUMULATOR_TYPE_PHYSICAL,
		m_hmmManager->getContextModelingOrderHMM(),m_hmmManager->getContextModelingOrderHMMCW()) == false) {
		return false;
	}
	
	// create the Forward-Backward objects
	// (a) numerator statistics
	m_forwardBackwardX = new ForwardBackwardX(m_phoneSet,m_lexiconManager,m_hmmManager,m_hmmManager,
		m_fForwardPruningBeam,m_fBackwardPruningBeam,m_iTrellisMaxSize,m_bTrellisCache,m_iTrellisCacheMaxSize,m_log);
	// (b) denominator statistics
	m_forwardBackward = new ForwardBackward(m_phoneSet,m_hmmManager,m_hmmManager,
		m_fForwardPruningBeam,m_fBackwardPruningBeam,m_iTrellisMaxSize,m_bTrellisCache,m_iTrellisCacheMaxSize,m_log);
	
	// transcription properties
	if (m_strFileOptionalSymbols != NULL) {
		m_vLexUnitOptional.push_back(m_lexiconManager->getLexUnitSilence());
	}
	
	m_iCovarianceModeling = m_hmmManager->getCovarianceModelling();
	m_iGaussianComponents = m_hmmManager->getNumberGaussianComponents();
	m_iHMMStates = m_hmmManager->getNumberHMMStatesPhysical();	

	return true;
}
		
// accumulate statistics
bool DTAccumulator::accumulate() {

	// make sure the HMMs are already initialized
	if (m_hmmManager->areInitialized() == false) {
		return false;
	}
	if (m_hmmManager->areInitialized() == false) {
		return false;
	}
	
	char strFileFeatures[MAX_PATH_LENGTH+1];
	char strFileLattice[MAX_PATH_LENGTH+1];
	char strFileAccumulators[MAX_PATH_LENGTH+1];
	char strFileAux[MAX_PATH_LENGTH+1];
	char strMessage[MAX_LOG_MESSAGE_LENGTH+1];	
	
	int iErrorCode;
	double dLikelihoodTotalNum = 0.0;
	double dLikelihoodTotalNumAcoustic = 0.0;
	double dLikelihoodTotalDen = 0.0;
	long iFeatureVectorsTotal = 0;
	long iFeatureVectorsUsedTotal = 0;
		
   double dBegin = Time::getTimeMilliseconds();
		
	// empty the accumulators
	m_hmmManager->resetAccumulators();
	
	m_bMMI = false;
	if (strcmp(m_strObjectiveFunction,DISCRIMINATIVE_TRAINING_OBJECTIVE_FUNCTION_BMMI) == 0) {
		m_bMMI = true;
	}
	
	// create accumulators for each HMM-state and Gaussian component	
	for(int i=0 ; i < m_hmmManager->getNumberHMMStatesPhysical() ; ++i) {
		for(int g=0 ; g < m_hmmManager->getHMMState(i)->getGaussianComponents() ; ++g) {
			unsigned int iKey = Accumulator::getPhysicalAccumulatorKey(i,g);
			// numerator
			Accumulator *accumulatorNum = new Accumulator(m_iFeatureDimensionality,
				m_hmmManager->getCovarianceModelling(),i,g);
			m_mAccumulatorNum.insert(MAccumulatorPhysical::value_type(iKey,accumulatorNum));
			// denominator
			Accumulator *accumulatorDen = new Accumulator(m_iFeatureDimensionality,
				m_hmmManager->getCovarianceModelling(),i,g);
			m_mAccumulatorDen.insert(MAccumulatorPhysical::value_type(iKey,accumulatorDen));
		}
	}
	
	// precompute constants used to speed-up emission probability computation
	if (m_hmmManager->precomputeConstants() == false) {
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
	
		// (2.1) load the features for the estimation
		sprintf(strFileFeatures,"%s%c%s",m_strFolderFeatures,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
		FeatureFile *featureFile = new FeatureFile(strFileFeatures,MODE_READ,FORMAT_FEATURES_FILE_DEFAULT,m_iFeatureDimensionality);
		featureFile->load();
		int iFeatureVectors = 0;
		float *fFeatures = (float*)featureFile->getFeatureVectors(&iFeatureVectors);
		delete featureFile;
		
		iFeatureVectorsTotal += iFeatureVectors;	
		
		// process the utterance using Forward-Backward (get the occupation counts)
		/*double dLikelihoodNum = -DBL_MAX;
		Alignment *alignmentNum = m_forwardBackwardX->processUtterance((*it)->vLexUnit,m_bMultiplePronunciations,
			m_vLexUnitOptional,fFeatures,fFeatures,iFeatureVectors,&dLikelihoodNum,iErrorCode);
		if (iErrorCode != UTTERANCE_PROCESSED_SUCCESSFULLY) {
			iFeatureVectorsTotal -= iFeatureVectors;
			delete [] fFeatures;	
			sprintf(strMessage,"unable to process utterance: \"%s\", reason: %s",strFileFeatures,
			m_forwardBackwardX->getErrorMessage(iErrorCode));
			m_log->logInformation(strMessage);
			continue;
		}*/
		
		// (2.2) load the hypothesis lattice
		sprintf(strFileAux,"%s%c%s",m_strFolderLattices,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
		FileUtils::replaceExtension(strFileLattice,strFileAux,"bin");
		printf("lattice: %s\n",strFileLattice);
		HypothesisLattice *lattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_log);
		try {
			lattice->load(strFileLattice,FILE_FORMAT_BINARY);
		} catch (ExceptionBase &e) {
			std::cerr << e.what() << std::endl;
			WARNING << "unable to load the lattice: " << strFileLattice;
			iFeatureVectorsTotal -= iFeatureVectors;
			delete [] fFeatures;
			continue;
		}
		
		//lattice->printProperties();
		// check lattice properties
		if ((lattice->isProperty(LATTICE_PROPERTY_AM_PROB) == false) ||
			(lattice->isProperty(LATTICE_PROPERTY_LM_PROB) == false) ||
			(lattice->isProperty(LATTICE_PROPERTY_INSERTION_PENALTY) == false) ||
			(lattice->isProperty(LATTICE_PROPERTY_HMMS) == false) ||
			(lattice->isProperty(LATTICE_PROPERTY_PHONE_ALIGN) == false)) {
			return false;
		}
		
		//LatticeDepth *depth = lattice->computeDepth();
		//printf("depth: %12.4f\n",depth->fDepth);
		
		//lattice->store("./lattice.txt",FILE_FORMAT_TEXT);
		
		// mark best path
		m_lexiconManager->removeNonStandardLexUnits((*it)->vLexUnit);
		LatticeWER *latticeWER = lattice->computeWER((*it)->vLexUnit);
		//HypothesisLattice::print(latticeWER);
		
		// get the best-path from the lattice (transcription) and compute numerator stats from it
		//BestPath *bestPath = lattice->getBestPath();
		//bestPath->print();
		//m_lexiconManager->print((*it)->vLexUnit);
		VLPhoneAlignment *vLPhoneAlignment = lattice->getBestPathAlignment();
		if (vLPhoneAlignment == NULL) {
			return false;
		}
		
		// compute numerator statistics
		iErrorCode = -1;
		double dLikelihoodNum = -DBL_MAX;		
		Alignment *alignmentNum = m_forwardBackward->processPhoneAlignment(fFeatures,iFeatureVectors,vLPhoneAlignment,dLikelihoodNum,iErrorCode);
		
		// get denominator statistics from the lattice
		iErrorCode = -1;
		double dLikelihoodDen = -DBL_MAX;		
		MOccupation *mOccupationDen = m_forwardBackward->processLattice(lattice,
			fFeatures,iFeatureVectors,m_fScaleAM,m_fScaleLM,dLikelihoodDen,m_bMMI,m_fBoostingFactor,iErrorCode);	
		
		// perform statistics cancellation between numerator and denominator
		if (m_bCanceledStatistics) {
			statisticsCancellation(alignmentNum,mOccupationDen);
		}
		
		Alignment *alignmentDen = ForwardBackward::getAlignment(mOccupationDen,iFeatureVectors);
		
		// accumulate statistics for both numerator and denominator
		accumulate(alignmentNum,fFeatures,iFeatureVectors,true);
		accumulate(alignmentDen,fFeatures,iFeatureVectors,false);	
		
		// get the best path with updated am-scores, lm-prob and insertion penalties
		BestPath *bestPath = lattice->getBestPath();
		assert(bestPath != NULL);
		LBestPathElement *lBestPathElement = bestPath->getBestPathElements();
		double dLMIP = 0.0;
		for(LBestPathElement::iterator it = lBestPathElement->begin() ; it != lBestPathElement->end() ; ++it) {
			dLMIP += (*it)->fScoreLanguageModel+m_fScaleAM*(*it)->fInsertionPenalty;
		}
		delete bestPath;
		dLikelihoodTotalNumAcoustic += dLikelihoodNum;
		// apply acoustic scaling and add the lm and insertion penalty scores
		dLikelihoodNum *= m_fScaleAM;
		dLikelihoodNum += dLMIP;
		
		dLikelihoodTotalNum += dLikelihoodNum;
		dLikelihoodTotalDen += dLikelihoodDen;
		
		//printf("%12.4f %12.4f\n",dLikelihoodNum,dLikelihoodDen);
		
		// clean-up
		delete latticeWER;
		delete vLPhoneAlignment;
		delete lattice;
		delete alignmentNum;
		delete alignmentDen;
		delete mOccupationDen;
		delete [] fFeatures;
	
		if (iFeatureVectorsTotal >= 100*60*30) {
			//break;
		}
		
		if (iUtterance == 5) {
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
	
	iFeatureVectorsUsedTotal = iFeatureVectorsTotal;
	
	// get the iteration end time
	double dEnd = Time::getTimeMilliseconds();
	double dMillisecondsInterval = dEnd - dBegin;
	
	// compute the Real Time Factor of the reestimation process
	float fRTF = (dMillisecondsInterval/10.0)/((float)iFeatureVectorsTotal);
	
	int iGaussians = m_hmmManager->getNumberGaussianComponents();
	float fLikelihoodFrameNum = dLikelihoodTotalNum/((float)iFeatureVectorsUsedTotal);
	float fLikelihoodFrameDen = dLikelihoodTotalDen/((float)iFeatureVectorsUsedTotal);
	// compute audio available for training
	int iHours,iMinutes,iSeconds;
	Time::convertHundredths((double)iFeatureVectorsTotal,iHours,iMinutes,iSeconds);
	// compute audio used
	int iHoursUsed,iMinutesUsed,iSecondsUsed;
	Time::convertHundredths((double)iFeatureVectorsUsedTotal,iHoursUsed,iMinutesUsed,iSecondsUsed);
	
	// show the accumulation information
	printf(" likelihood= (%.4f) %.4f %.4f %.4f [%8d Gauss][RTF=%.4f][%d:%02d'%02d''][%d:%02d'%02d'']\n",
		dLikelihoodTotalNumAcoustic,dLikelihoodTotalNum,dLikelihoodTotalDen,dLikelihoodTotalNum-dLikelihoodTotalDen,
		iGaussians,fRTF,iHours,iMinutes,iSeconds,iHoursUsed,iMinutesUsed,iSecondsUsed);	
	
	// dump the accumulators
	Accumulator::storeAccumulators(m_strFileAccumulatorsNum,m_iFeatureDimensionality,m_iCovarianceModeling,
		m_iHMMStates,m_iGaussianComponents,m_mAccumulatorNum);
	Accumulator::storeAccumulators(m_strFileAccumulatorsDen,m_iFeatureDimensionality,m_iCovarianceModeling,
		m_iHMMStates,m_iGaussianComponents,m_mAccumulatorDen);
	
	// destroy the accumulators
	Accumulator::destroy(m_mAccumulatorNum);
	Accumulator::destroy(m_mAccumulatorDen);

	return true;
}

// statistics cancellation (between numerator and denominator)
void DTAccumulator::statisticsCancellation(Alignment *alignmentNum, MOccupation *mOccupationDen) {

	for(int t=0 ; t < alignmentNum->getFrames() ; ++t) {
		FrameAlignment *frameAlignment = alignmentNum->getFrameAlignment(t);
		for(FrameAlignment::iterator it = frameAlignment->begin() ; it != frameAlignment->end() ; ++it) {
			double dOccupationNum = (*it)->dOccupation;
			MOccupation::iterator jt = mOccupationDen->find(pair<int,int>(t,(*it)->iHMMState));
			if (jt != mOccupationDen->end()) {
				double dOccupationDen = jt->second;
				double dOccupationShared = min(dOccupationNum,dOccupationDen);
				jt->second -= dOccupationShared;
				(*it)->dOccupation -= dOccupationShared;
			}
		}
	}
}

// accumulate statistics
void DTAccumulator::accumulate(Alignment *alignment, float *fFeatures, int iFeatures, bool bNumerator) {

	Accumulator *accumulator = NULL;
	double dOccupationTotal = 0.0;
	for(int t=0 ; t < alignment->getFrames() ; ++t) {
		FrameAlignment *frameAlignment = alignment->getFrameAlignment(t);
		for(FrameAlignment::iterator it = frameAlignment->begin() ; it != frameAlignment->end() ; ++it) {
			double dOccupationNum = (*it)->dOccupation;
			HMMState *hmmState = m_hmmManager->getHMMState((*it)->iHMMState);
			// get Gaussian occupation from the mixture occupation
			
			// (1) compute the mixture likelihood (all Gaussian components)
			double dLikelihoodTotal = -DBL_MAX;
			int iGaussianComponents = hmmState->getGaussianComponents();
			double *dLikelihoodGaussian = new double[iGaussianComponents];
			for(int iGaussian = 0 ; iGaussian < iGaussianComponents ; ++iGaussian) {
				dLikelihoodGaussian[iGaussian] = hmmState->computeEmissionProbabilityGaussian(iGaussian,
					fFeatures+(t*m_iFeatureDimensionality),-1);
				dLikelihoodGaussian[iGaussian] += log(hmmState->getGaussian(iGaussian)->fWeight);
				dLikelihoodTotal = NumericalFunctions::logAddition(dLikelihoodTotal,dLikelihoodGaussian[iGaussian]);
			}
			// (2) accumulate statistics for each mixture component
			for(int iGaussian = 0 ; iGaussian < iGaussianComponents ; ++iGaussian) {
				
				double dProbGaussian = exp(dLikelihoodGaussian[iGaussian]-dLikelihoodTotal);
				assert(dProbGaussian >= 0.0);
				double dOccupationGaussian = dOccupationNum*dProbGaussian;
				unsigned int iKey = Accumulator::getPhysicalAccumulatorKey(hmmState->getId(),iGaussian);
				if (bNumerator) {
					accumulator = m_mAccumulatorNum[iKey];
				} else {
					accumulator = m_mAccumulatorDen[iKey];	
				}
				accumulator->accumulateObservation(fFeatures+(t*m_iFeatureDimensionality),dOccupationGaussian);
				dOccupationTotal += dOccupationGaussian;
			}
			delete [] dLikelihoodGaussian;
		}
	}	
	printf("accumulated: %12.6f\n",dOccupationTotal);
	
	return;
}




