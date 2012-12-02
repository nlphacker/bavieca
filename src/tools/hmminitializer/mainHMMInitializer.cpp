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


#include <iostream>
#include <cstdlib>

#include "CommandLineManager.h"
#include "ConfigurationFeatures.h"
#include "FeatureExtractor.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "HMMInitializer.h"
#include "HMMManager.h"
#include "MLEstimator.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the HMM estimation tool: "hmminitializer"
int main(int argc, char *argv[]) {

	// define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("hmminitializer",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-fea","feature folder",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-mlf","master label file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-met","HMM initialization method",
		PARAMETER_TYPE_STRING,true,"flatStart|alignment","flatStart");
	m_commandLineManager->defineParameter("-cov","covariance type",PARAMETER_TYPE_STRING,true,"diagonal|full","diagonal");	
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);	
	
	// parse command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get the parameters
	const char *m_strFileFeatureConfiguration = m_commandLineManager->getParameterValue("-cfg");
	const char *m_strFolderFeatures = m_commandLineManager->getParameterValue("-fea");
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileLexicon = m_commandLineManager->getParameterValue("-lex");
	const char *m_strFileMLF = m_commandLineManager->getParameterValue("-mlf");
	const char *m_strHMMInitializationMethod = m_commandLineManager->getParameterValue("-met");
	const char *m_strCovarianceType = m_commandLineManager->getParameterValue("-cov");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	
	// load the feature configuration file
	ConfigurationFeatures *m_configurationFeatures = new ConfigurationFeatures(m_strFileFeatureConfiguration);
	m_configurationFeatures->load();
	// get the feature dimensionality
	int m_iFeatureDimensionality = m_configurationFeatures->getDimensionality();
	delete m_configurationFeatures;
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
	
	// load the lexicon
   LexiconManager *m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
   m_lexiconManager->load();
   //m_lexiconManager->print(true);
	
	// load the Master Label File
	MLFFile *m_mlfFile = new MLFFile(m_lexiconManager,m_strFileMLF,MODE_READ);	
	m_mlfFile->load();
	
	// initialize the acoustic models
	int iCovarianceType = Gaussian::getCovarianceModellingType(m_strCovarianceType);	
	HMMInitializer *m_hmmInitializer = new HMMInitializer(m_iFeatureDimensionality,iCovarianceType,m_phoneSet);
	
	HMMManager *m_hmmManager = NULL;
	
	// initialization of HMM-parameters
	// (flat-start)
	if (strcmp(m_strHMMInitializationMethod,AM_INITIALIZATION_METHOD_FLATSTART) == 0) {
	
		// initialize the acoustic models to the global distribution of the data
		m_hmmManager = m_hmmInitializer->initializeModelsFlatStart(m_mlfFile->getUtterances(),m_strFolderFeatures);
		assert(m_hmmManager);
		
		// store the models to disk
		m_hmmManager->store(m_strFileModels);
		
		delete m_hmmManager;
	}
	// (alignment)
	else {
		assert(strcmp(m_strHMMInitializationMethod,AM_INITIALIZATION_METHOD_ALIGNMENT) == 0);	
		BVC_ERROR << "initialization method not supported yet";
	}
		
	// clean-up
	delete m_hmmInitializer;
	delete m_mlfFile;
	delete m_lexiconManager;
	delete m_phoneSet;
	delete m_commandLineManager;	
	
	return 0;
}

