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

#include "Accumulator.h"
#include "HMMManager.h"
#include "PhoneSet.h"
#include "CommandLineManager.h"
#include "LDAEstimator.h"

using namespace Bavieca;
 
// main for the tool "ldaestimator"
int main(int argc, char *argv[]) {

	try {

		// define the parameters
		CommandLineManager commandLineManager("hldaestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-acc","input accumulator filelist",PARAMETER_TYPE_FILE,true);	
		commandLineManager.defineParameter("-dim","target dimensionality",PARAMETER_TYPE_INTEGER,false);	
		commandLineManager.defineParameter("-out","output transform",PARAMETER_TYPE_FILE,false);	
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// retrieve the parameters
		/*const char *strFileFeatureConfiguration = commandLineManager.getParameterValue("-cfg");
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
			hmmManager.getContextModelingOrderHMM(),hmmManager.getContextModelingOrderHMMCW());	*/
		
		// create the logging object
		/*Log log(VERBOSITY_LEVEL_INFORMATION);
	
		// initialize the HLDA estimator
		HLDAEstimator *hldaEstimator = new HLDAEstimator(hmmManager,strFileAccList,iDimensionalityReduction,iIterationsTransformUpdate,iIterationsParameterUpdate,strFolderOutput,log);
		
		// do the actual estimation
		if (hldaEstimator->estimate() == false) {
			printf("Error: a problem was encountered during the HLDA estimation\n");
			return -1;
		}*/
		
	} catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	

	return 0;
}

