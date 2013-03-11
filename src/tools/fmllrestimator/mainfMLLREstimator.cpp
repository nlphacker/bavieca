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
#include "FMLLREstimator.h"
#include "PhoneSet.h"
#include "TimeUtils.h"
#include "Transform.h"

using namespace Bavieca;

// main for the tool "mapEstimator"
int main(int argc, char *argv[]) {

	try {
	
		// (1) define command line parameters
		CommandLineManager commandLineManager("fmllrestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","bootstrap acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-bat","batch file containing pairs (featureFile alignmentFile)",
			PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");
		commandLineManager.defineParameter("-tra","feature transform",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-bst","whether to assign all occupation to best scoring Gaussian component",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// load the parameter values
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		const char *strFileBatch = commandLineManager.getParameterValue("-bat");
		const char *strFormat = commandLineManager.getParameterValue("-for");
		const char *strFileTransform = commandLineManager.getParameterValue("-tra");
		bool bBestComponentOnly = commandLineManager.getBoolParameterValue("-bst");
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		hmmManager.initializeDecoding();
		
		// get starting time
		double dBegin = TimeUtils::getTimeMilliseconds();	
			
		// create the fMLLR estimator
		FMLLREstimator fMLLREstimator(&phoneSet,&hmmManager,20,bBestComponentOnly);
		fMLLREstimator.initializeEstimation();
		double dLikelihood = 0.0;
		fMLLREstimator.feedAdaptationData(strFileBatch,strFormat,&dLikelihood,true);	
		
		double dEndFeed = TimeUtils::getTimeMilliseconds();
		
		Transform *transform = fMLLREstimator.estimateTransform(NULL);	
		assert(transform);
		
		// store the transform
		transform->store(strFileTransform);
		delete transform;
		
		fMLLREstimator.uninitializeEstimation();
		
		double dEndComputation = TimeUtils::getTimeMilliseconds();
		
		// show the processing time
		printf("collecting data:      %6.2f seconds\n",(dEndFeed-dBegin)/1000.0);
		printf("computing transforms: %6.2f seconds\n",(dEndComputation-dEndFeed)/1000.0);
   
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
	
	return 0;
}
