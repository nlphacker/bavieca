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
#include "DTEstimator.h"
#include "MLEstimator.h"
#include "MLFFile.h"
#include "PhoneSet.h"

using namespace Bavieca;
 
// main for the tool "DTEstimator"
int main(int argc, char *argv[]) {

	try {
	
		// (1) define command line parameters
		CommandLineManager commandLineManager("dtestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-accNum","input accumulator filelist (numerator)",
			PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-accDen","input accumulator filelist (denominator)",
			PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-cov","covariance flooring ratio",PARAMETER_TYPE_FLOAT,true,NULL,"0.05");
		commandLineManager.defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-E","learning rate constant",PARAMETER_TYPE_FILE,true,NULL,"2.0");
		commandLineManager.defineParameter("-I","I-smoothing",PARAMETER_TYPE_STRING,true,"none|prev","none");
		commandLineManager.defineParameter("-tau","I-smoothing constant",PARAMETER_TYPE_FLOAT,true,NULL,"100.0");
		
		// parse the command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		const char *strFileAccListNum = commandLineManager.getParameterValue("-accNum");
		const char *strFileAccListDen = commandLineManager.getParameterValue("-accDen");
		float fCovarianceFlooringRatio = atof(commandLineManager.getParameterValue("-cov"));
		const char *strFileModelsOutput = commandLineManager.getParameterValue("-out");
		float fE = atof(commandLineManager.getParameterValue("-E"));
		const char *strISmoothing = commandLineManager.getParameterValue("-I");
		float fTau = -1.0;
		if (commandLineManager.isParameterSet("-tau")) {
			fTau = atof(commandLineManager.getParameterValue("-tau"));
		}
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
			
		// load the HMMs
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_ESTIMATION);	
		hmmManager.load(strFileModels);
	
		// initialize the HMMs
		hmmManager.initializeEstimation(ACCUMULATOR_TYPE_PHYSICAL,UCHAR_MAX,UCHAR_MAX);
		
		// create the discriminative training estimator
		DTEstimator dtEstimator(&hmmManager);
		
		// estimate parameters
		dtEstimator.estimateParameters(strFileAccListNum,strFileAccListDen,fE,strISmoothing,fTau,true);
		
		// floor covariances
		dtEstimator.floorCovariances(fCovarianceFlooringRatio);
		
		// create the output HMMs
		hmmManager.store(strFileModelsOutput);
	
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	

	return 0;
}
