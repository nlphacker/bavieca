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

	try {
	
		// (1) define command line parameters
		CommandLineManager commandLineManager("mllrestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","bootstrap acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-bat","batch file containing pairs (featureFile alignmentFile)",
			PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");
		commandLineManager.defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-tra","folder to store the transforms",PARAMETER_TYPE_FOLDER,true);	
		commandLineManager.defineParameter("-rgt","file containing the regression tree",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-occ","minimum number of frames to compute a transform",
			PARAMETER_TYPE_INTEGER,true,"[1000|6000]","3500");
		commandLineManager.defineParameter("-gau","minimum number of observed Gaussian distributions to compute a transform",PARAMETER_TYPE_INTEGER,true,NULL,"1");
		commandLineManager.defineParameter("-bst","whether to assign all occupation to best scoring Gaussian component",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");
		commandLineManager.defineParameter("-cov","whether to compute the covariance transformation",PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");
		
		// (2) parse command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		const char *strFileBatch = commandLineManager.getParameterValue("-bat");
		const char *strAlignmentFormat = commandLineManager.getParameterValue("-for");
		const char *strFileModelsOutput = commandLineManager.getParameterValue("-out");	
		const char *strFolderTransforms = commandLineManager.getParameterValue("-tra");	
		const char *strFileRegressionTree = NULL;
		if (commandLineManager.isParameterSet("-rgt")) {
			strFileRegressionTree = commandLineManager.getParameterValue("-rgt");
		}
		float fMinimumOccupationTransform = -1.0;
		if (commandLineManager.isParameterSet("-occ")) {
			fMinimumOccupationTransform = atof(commandLineManager.getParameterValue("-occ"));
		}
		int iMinimumGaussianComponentsObserved = -1;
		if (commandLineManager.isParameterSet("-gau")) {
			iMinimumGaussianComponentsObserved = atoi(commandLineManager.getParameterValue("-gau"));
		}
		bool bBestComponentOnly = commandLineManager.getBoolParameterValue("-bst");
		bool bMeanOnly = commandLineManager.getBoolParameterValue("-cov");
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		hmmManager.initializeDecoding();
			
		// get starting time
		double dBegin = TimeUtils::getTimeMilliseconds();
		
		// create the MLLR object
		MLLRManager mllrManager(&phoneSet,&hmmManager,strFileRegressionTree,fMinimumOccupationTransform,
			iMinimumGaussianComponentsObserved,bBestComponentOnly,bMeanOnly);
		
		// initialize MLLR
		mllrManager.initialize();
		
		// likelihood
		double dLikelihoodUnadapted = 0.0;
		double dLikelihoodAdapted = 0.0;
		
		// process the batch file
		mllrManager.feedAdaptationData(strFileBatch,strAlignmentFormat,&dLikelihoodUnadapted,true);
		
		double dEndFeed = TimeUtils::getTimeMilliseconds();	
		
		// compute the transforms from the adaptation data
		mllrManager.computeTransforms();
		
		double dEndComputation = TimeUtils::getTimeMilliseconds();	
		
		// update the HMM-state parameters using the computed tranforms
		mllrManager.applyTransforms();
		
		double dEndFinal = TimeUtils::getTimeMilliseconds();	
		
		// show the processing time
		printf("collecting data:      %6.2f seconds\n",(dEndFeed-dBegin)/1000.0);
		printf("computing transforms: %6.2f seconds\n",(dEndComputation-dEndFeed)/1000.0);
		printf("applying transforms:  %6.2f seconds\n",(dEndFinal-dEndComputation)/1000.0);
		
		// process the batch file
		mllrManager.feedAdaptationData(strFileBatch,strAlignmentFormat,&dLikelihoodAdapted,true);
		
		// write the transforms to disk
		if (strFolderTransforms) {
			mllrManager.storeTransforms(strFolderTransforms);
		}
		
		// write the adapted acoustic models to disk
		if (strFileModelsOutput) {
			hmmManager.store(strFileModelsOutput);
		}	

	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
		
	return 0;
}
