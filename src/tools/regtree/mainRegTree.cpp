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
#include "RegressionTree.h"
#include "TimeUtils.h"

using namespace Bavieca;

// main for the tool "regtree"
int main(int argc, char *argv[]) {

	try {
	
		// (1) define command line parameters
		CommandLineManager commandLineManager("regtree",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-met","clustering method",PARAMETER_TYPE_STRING,true,"kMeans|EM","EM");
		commandLineManager.defineParameter("-rgc","number of regression classes (base-classes)",PARAMETER_TYPE_INTEGER,true,"[1|1000]","50");
		commandLineManager.defineParameter("-gau","minimum number of Gaussian components per base-class",PARAMETER_TYPE_INTEGER,true,NULL,"50");
		commandLineManager.defineParameter("-out","file to store the regression tree",PARAMETER_TYPE_FILE,false);
		
		// (2) parse command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		// clustering method
		unsigned char iClusteringMethod;
		if (strcmp(commandLineManager.getParameterValue("-met"),"kMeans") == 0) {
			iClusteringMethod = CLUSTERING_METHOD_KMEANS;
		} else {
			iClusteringMethod = CLUSTERING_METHOD_EM;
		}	
		int iRegressionClasses = atoi(commandLineManager.getParameterValue("-rgc"));
		int iMinimumGaussianComponentsBaseClass = atoi(commandLineManager.getParameterValue("-gau"));
		const char *strFileRegressionTree = commandLineManager.getParameterValue("-out");
	
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		hmmManager.initializeDecoding();
		
		// build the regression tree
		RegressionTree regressionTree(&hmmManager);
		regressionTree.build(iRegressionClasses,iClusteringMethod,
			iMinimumGaussianComponentsBaseClass);
			
		// store the regression tree to disk
		regressionTree.store(strFileRegressionTree);
		
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
		
	return 0;
}
