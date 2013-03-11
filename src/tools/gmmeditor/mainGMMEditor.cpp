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
#include "GMMEditor.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "HMMInitializer.h"
#include "HMMManager.h"
#include "MLEstimator.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the HMM edition tool: "hmmeditor"
int main(int argc, char *argv[]) {

	try {

		// define command line parameters
		CommandLineManager commandLineManager("hmmeditor",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-acc","accumulators",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-dbl","mixture doubling",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
		commandLineManager.defineParameter("-inc","mixture increment",PARAMETER_TYPE_INTEGER,true);
		commandLineManager.defineParameter("-crt","splitting criterion",PARAMETER_TYPE_STRING,true,"covariance|weight","covariance");
		commandLineManager.defineParameter("-occ","minimum component occupation",PARAMETER_TYPE_FLOAT,true,NULL,"100.0");
		commandLineManager.defineParameter("-wgh","minimum component weight",PARAMETER_TYPE_FLOAT,true,NULL,"0.00001");
		commandLineManager.defineParameter("-eps","epsilon (mean perturbation)",PARAMETER_TYPE_FLOAT,true,NULL,"0.05");
		commandLineManager.defineParameter("-mrg","merge components below minimum occ/wgh",PARAMETER_TYPE_BOOLEAN,true,NULL,"yes");	
		commandLineManager.defineParameter("-cov","covariance flooring ratio",PARAMETER_TYPE_FLOAT,true,NULL,"0.05");	
		commandLineManager.defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-vrb","verbose",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");	
		
		// parse command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
			
		// get the parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		const char *strFileAccList = commandLineManager.getParameterValue("-acc");	
		bool bMixtureDoubling = CommandLineManager::str2bool(commandLineManager.getParameterValue("-dbl"));
		int iMixtureIncrement = 0;
		if (commandLineManager.isParameterSet("-inc") == true) {	
			iMixtureIncrement = atoi(commandLineManager.getParameterValue("-inc"));
		}
		const char *strSplittingCriterion = commandLineManager.getParameterValue("-crt");
		int iSplittingCriterion = GMMEditor::getSplittingCriterion(strSplittingCriterion);
		float fMinimumComponentOccupation = atof(commandLineManager.getParameterValue("-occ"));
		float fMinimumComponentWeight = atof(commandLineManager.getParameterValue("-wgh"));
		float fEpsilon = atof(commandLineManager.getParameterValue("-eps"));
		bool bMixtureMerging = CommandLineManager::str2bool(commandLineManager.getParameterValue("-mrg"));
		float fCovarianceFlooringRatio = atof(commandLineManager.getStrParameterValue("-cov"));
		const char *strFileModelsOutput = commandLineManager.getParameterValue("-out");
		//bool bVerbose = CommandLineManager::str2bool(commandLineManager.getParameterValue("-vrb"));
		
		// load the phonetic symbol set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the HMMs
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_ESTIMATION);	
		hmmManager.load(strFileModels);
		
		// load the accumulators
		AccMetadata metadata;
		MAccumulatorPhysical mAccumulatorPhysical;
		Accumulator::loadAccumulatorList(strFileAccList,mAccumulatorPhysical,metadata);
		
		// create the GMM editor
		GMMEditor gmmEditor(&hmmManager,&mAccumulatorPhysical);
		gmmEditor.initialize();	
		
		// mixture merging
		if (bMixtureMerging) {
			gmmEditor.mixtureMerge(fMinimumComponentOccupation,fMinimumComponentWeight,fCovarianceFlooringRatio);
		}	
		// mixture doubling
		if (bMixtureDoubling) {
			gmmEditor.mixtureDouble(fMinimumComponentOccupation,fEpsilon);
		}
		// mixture increment
		else {
			gmmEditor.mixtureIncrement(iMixtureIncrement,iSplittingCriterion,false,
				fMinimumComponentOccupation,fEpsilon);
		}
			
		// store the models after the Gaussian merging process is done	
		hmmManager.store(strFileModelsOutput);
	
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;	
}
