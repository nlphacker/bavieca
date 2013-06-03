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

	try {

		// define command line parameters
		CommandLineManager commandLineManager("hmminitializer",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-fea","feature folder",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-mlf","master label file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-met","HMM initialization method",
			PARAMETER_TYPE_STRING,true,"flatStart|alignment","flatStart");
		commandLineManager.defineParameter("-cov","covariance type",PARAMETER_TYPE_STRING,true,"diagonal|full","diagonal");	
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);	
		
		// parse command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the parameters
		const char *strFileFeatureConfiguration = commandLineManager.getParameterValue("-cfg");
		const char *strFolderFeatures = commandLineManager.getParameterValue("-fea");
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileLexicon = commandLineManager.getParameterValue("-lex");
		const char *strFileMLF = commandLineManager.getParameterValue("-mlf");
		const char *strHMMInitializationMethod = commandLineManager.getParameterValue("-met");
		const char *strCovarianceType = commandLineManager.getParameterValue("-cov");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		
		// load the feature configuration file
		ConfigurationFeatures configurationFeatures(strFileFeatureConfiguration);
		configurationFeatures.load();
		// get the feature dimensionality
		int iFeatureDimensionality = configurationFeatures.getDimensionality();
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet); 
		lexiconManager.load();
		
		// load the Master Label File
		MLFFile mlfFile(&lexiconManager,strFileMLF,MODE_READ);	
		mlfFile.load();
		
		// initialize the acoustic models
		int iCovarianceType = Gaussian::getCovarianceModellingType(strCovarianceType);	
		HMMInitializer hmmInitializer(iFeatureDimensionality,iCovarianceType,&phoneSet);
		
		// initialization of HMM-parameters
		// (flat-start)
		if (strcmp(strHMMInitializationMethod,AM_INITIALIZATION_METHOD_FLATSTART) == 0) {
		
			// initialize the acoustic models to the global distribution of the data
			HMMManager *hmmManager = hmmInitializer.initializeModelsFlatStart(mlfFile.getUtterances(),strFolderFeatures);
			assert(hmmManager);
			hmmManager->store(strFileModels);	
			delete hmmManager;
		}
		// (alignment)
		else {
			assert(strcmp(strHMMInitializationMethod,AM_INITIALIZATION_METHOD_ALIGNMENT) == 0);	
			BVC_ERROR << "initialization method not supported yet";
		}
	
	} catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
	
	return 0;
}

