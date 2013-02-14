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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "CommandLineManager.h"
#include "HMMManager.h"
#include "HLDAEstimator.h"
#include "PhoneSet.h"
#include "RegressionTree.h"
#include "Transform.h"

using namespace std;

#include <string>
#include <map>

using namespace Bavieca;

// main for the tool "hmmx"
int main(int argc, char *argv[])
{	
	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("hmmx",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-tra","model transform",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-rgt","regression-tree",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-in","inout acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
		
		// (2) process command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// load parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileTransform = commandLineManager.getParameterValue("-tra");
		const char *strFileRegressionTree = commandLineManager.getParameterValue("-rgt");
		const char *strFileModelsInput = commandLineManager.getParameterValue("-in");
		const char *strFileModelsOutput = commandLineManager.getParameterValue("-out");
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// regression tree based transform?	
		if (strFileRegressionTree != NULL) {
		
			// load the acoustic models
			HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
			hmmManager.load(strFileModelsInput);
			hmmManager.initializeDecoding();	
				
			// load the regression tree
			RegressionTree regressionTree(&hmmManager);
			regressionTree.load(strFileRegressionTree);
			
			// load the transforms
			regressionTree.loadTransforms(strFileTransform);
			
			// apply the transforms
			regressionTree.applyTransforms();
			
			// store the output acoustic models to disk
			hmmManager.store(strFileModelsOutput);
		}
		// regular model transform
		else {
		
			// load the input acoustic models
			HMMManager hmmManager(&phoneSet,HMM_PURPOSE_ESTIMATION);
			hmmManager.load(strFileModelsInput);
			
			// load the transformation
			Transform transform;
			transform.load(strFileTransform);
			
			// apply the actual transformation
			HLDAEstimator::applyTransform(transform.getTransform(),&hmmManager);
		
			// store the output acoustic models to disk
			hmmManager.store(strFileModelsOutput);
		}		
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
	
	return 0;
}



