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


#include "CommandLineManager.h"
#include "HMMManager.h"
#include "LexiconManager.h"
#include "MLLRManager.h"
#include "TimeUtils.h"

using namespace Bavieca;

// main for the tool "mllrestimator"
int main(int argc, char *argv[]) {
	
	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("mllrestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","bootstrap acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-bat","batch file containing pairs (featureFile alignmentFile)",
		PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");
	m_commandLineManager->defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-tra","folder to store the transforms",PARAMETER_TYPE_FOLDER,true);	
	m_commandLineManager->defineParameter("-rgt","file containing the regression tree",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-occ","minimum number of frames to compute a transform",
		PARAMETER_TYPE_INTEGER,true,"[1000|6000]","3500");
	m_commandLineManager->defineParameter("-gau","minimum number of observed Gaussian distributions to compute a transform",PARAMETER_TYPE_INTEGER,true,NULL,"1");
	m_commandLineManager->defineParameter("-bst","whether to assign all occupation to best scoring Gaussian component",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");
	m_commandLineManager->defineParameter("-cov","whether to compute the covariance transformation",PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");
	
	// (2) parse command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *m_strFileBatch = m_commandLineManager->getParameterValue("-bat");
	const char *m_strAlignmentFormat = m_commandLineManager->getParameterValue("-for");
	const char *m_strFileModelsOutput = m_commandLineManager->getParameterValue("-out");	
	const char *m_strFolderTransforms = m_commandLineManager->getParameterValue("-tra");	
	const char *m_strFileRegressionTree = NULL;
	if (m_commandLineManager->isParameterSet("-rgt") == true) {
		m_strFileRegressionTree = m_commandLineManager->getParameterValue("-rgt");
	}
	float m_fMinimumOccupationTransform = -1.0;
	if (m_commandLineManager->isParameterSet("-occ")) {
		m_fMinimumOccupationTransform = atof(m_commandLineManager->getParameterValue("-occ"));
	}
	int m_iMinimumGaussianComponentsObserved = -1;
	if (m_commandLineManager->isParameterSet("-gau")) {
		m_iMinimumGaussianComponentsObserved = atoi(m_commandLineManager->getParameterValue("-gau"));
	}
	bool m_bBestComponentOnly = m_commandLineManager->getBoolParameterValue("-bst");
	bool m_bMeanOnly = m_commandLineManager->getBoolParameterValue("-cov");
	
	// load the phone set
	PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
	m_phoneSet->load();
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeDecoding();
		
	// get starting time
	double dBegin = TimeUtils::getTimeMilliseconds();
	
	// create the MLLR object
	MLLRManager *m_mllrManager = new MLLRManager(m_phoneSet,m_hmmManager,m_strFileRegressionTree,
		m_fMinimumOccupationTransform,m_iMinimumGaussianComponentsObserved,m_bBestComponentOnly,m_bMeanOnly);
	
	// initialize MLLR
	m_mllrManager->initialize();
	
	// likelihood
	double dLikelihoodUnadapted = 0.0;
	double dLikelihoodAdapted = 0.0;
	
	// process the batch file
	m_mllrManager->feedAdaptationData(m_strFileBatch,m_strAlignmentFormat,&dLikelihoodUnadapted,true);
	
	double dEndFeed = TimeUtils::getTimeMilliseconds();	
	
	// compute the transforms from the adaptation data
	m_mllrManager->computeTransforms();
	
	double dEndComputation = TimeUtils::getTimeMilliseconds();	
	
	// update the HMM-state parameters using the computed tranforms
	m_mllrManager->applyTransforms();
	
	double dEndFinal = TimeUtils::getTimeMilliseconds();	
   
	// show the processing time
	printf("collecting data:      %6.2f seconds\n",(dEndFeed-dBegin)/1000.0);
	printf("computing transforms: %6.2f seconds\n",(dEndComputation-dEndFeed)/1000.0);
	printf("applying transforms:  %6.2f seconds\n",(dEndFinal-dEndComputation)/1000.0);
	
	// process the batch file
	m_mllrManager->feedAdaptationData(m_strFileBatch,m_strAlignmentFormat,&dLikelihoodAdapted,true);
	
	// write the transforms to disk
	if (m_strFolderTransforms) {
		m_mllrManager->storeTransforms(m_strFolderTransforms);
	}
	
	// write the adapted acoustic models to disk
	if (m_strFileModelsOutput) {
		m_hmmManager->store(m_strFileModelsOutput);
	}	
	
	// clean-up
	delete m_mllrManager;
	delete m_hmmManager;
	delete m_phoneSet;
	delete m_commandLineManager;
		
	return 0;
}
