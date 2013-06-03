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

#include "CommandLineManager.h"
#include "ConfigurationFeatures.h"
#include "HLDAEstimator.h"
#include "HMMManager.h"
#include "PhoneSet.h"

using namespace Bavieca;
 
// main for the tool "hldaestimator"
int main(int argc, char *argv[]) {

	try {

		// define the parameters
		CommandLineManager commandLineManager("hldaestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-acc","input accumulator filelist",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-itt","number of iterations for transform update",
			PARAMETER_TYPE_INTEGER,true,"[1|100]","10");	
		commandLineManager.defineParameter("-itp","number of iterations for parameter update",
			PARAMETER_TYPE_INTEGER,true,"[1|100]","10");	
		commandLineManager.defineParameter("-red","dimensionality reduction",PARAMETER_TYPE_INTEGER,true,NULL,"13");	
		commandLineManager.defineParameter("-out","output folder",PARAMETER_TYPE_FOLDER,false);	
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// retrieve the parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		const char *strFileAccList = commandLineManager.getParameterValue("-acc");
		int iIterationsTransformUpdate = atoi(commandLineManager.getParameterValue("-itt"));
		int iIterationsParameterUpdate = atoi(commandLineManager.getParameterValue("-itp"));
		int iDimensionalityReduction = atoi(commandLineManager.getParameterValue("-red"));
		const char *strFolderOutput = commandLineManager.getParameterValue("-out");
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
	
		// load the HMM-models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_ESTIMATION);	
		hmmManager.load(strFileModels);
		hmmManager.initializeEstimation(ACCUMULATOR_TYPE_PHYSICAL,
			hmmManager.getContextModelingOrderHMM(),hmmManager.getContextModelingOrderHMMCW());	
	
		// initialize the HLDA estimator
		HLDAEstimator hldaEstimator(&hmmManager,strFileAccList,iDimensionalityReduction,
			iIterationsTransformUpdate,iIterationsParameterUpdate,strFolderOutput);
		
		// do the actual estimation
		hldaEstimator.estimate();
	
	} catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	

	return 0;
}

