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
#include "FeatureExtractor.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "HMMManager.h"
#include "MLAccumulator.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the HMM estimation tool: "mlaccumulator"
int main(int argc, char *argv[]) {

	try {

		// define command line parameters
		CommandLineManager commandLineManager("mlaccumulator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-fea","feature folder",PARAMETER_TYPE_FOLDER,false);
		commandLineManager.defineParameter("-cfg","feature configuration file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-feaA","feature folder (stats accumulation)",PARAMETER_TYPE_FOLDER,true);	
		commandLineManager.defineParameter("-cfgA","feature configuration (stats accumulation)",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-covA","covariance modeling type (stats accumulation)",
			PARAMETER_TYPE_STRING,true,"diagonal|full","diagonal");
		commandLineManager.defineParameter("-mlf","master label file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-ww","within-word context modeling order for the accumulators",PARAMETER_TYPE_STRING,true);
		commandLineManager.defineParameter("-cw","cross-word context modeling order for the accumulators",PARAMETER_TYPE_STRING,true);
		commandLineManager.defineParameter("-opt","file containing optional symbols that can be inserted at word boundaries",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-pro","allow multiple pronunciations",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
		commandLineManager.defineParameter("-fwd","forward pruning",PARAMETER_TYPE_FLOAT,true,"[-100.0|-10.0]","-20");
		commandLineManager.defineParameter("-bwd","backward pruning",PARAMETER_TYPE_FLOAT,true,"[100.0|10000.0]","800");
		commandLineManager.defineParameter("-tre","maximum trellis size (MB)",PARAMETER_TYPE_INTEGER,true,NULL,"500");
		commandLineManager.defineParameter("-dAcc","file to dump accumulators",PARAMETER_TYPE_FILE,false);	
		
		// parse the command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
			
		// load the parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModelsAlignment = commandLineManager.getParameterValue("-mod");
		const char *strFileLexicon = commandLineManager.getParameterValue("-lex");
		const char *strFolderFeaturesAlignment = commandLineManager.getParameterValue("-fea"); 
		const char *strFileFeatureConfigurationAlignment = commandLineManager.getParameterValue("-cfg"); 
		const char *strFolderFeaturesAcc = commandLineManager.getParameterValue("-feaA"); 
		const char *strFileFeatureConfigurationAcc = commandLineManager.getParameterValue("-cfgA"); 
		int iCovarianceModellingAcc;
		if (commandLineManager.isParameterSet("-covA") == false) {
			iCovarianceModellingAcc = COVARIANCE_MODELLING_TYPE_DEFAULT;
		} else if (strcmp(commandLineManager.getParameterValue("-covA"),COVARIANCE_MODELLING_TYPE_FULL_STR) == 0) {
			iCovarianceModellingAcc = COVARIANCE_MODELLING_TYPE_FULL;
		} else {
			assert(strcmp(commandLineManager.getParameterValue("-covA"),COVARIANCE_MODELLING_TYPE_DIAGONAL_STR) == 0);
			iCovarianceModellingAcc = COVARIANCE_MODELLING_TYPE_DIAGONAL;
		}
		const char *strFileMLF = commandLineManager.getParameterValue("-mlf"); 
		unsigned char iAccumulatorType = ACCUMULATOR_TYPE_PHYSICAL;
		int iContextModelingOrderAccumulatorsWW = UCHAR_MAX;
		int iContextModelingOrderAccumulatorsCW = UCHAR_MAX;
		if (commandLineManager.isParameterSet("-ww") == true) {
			iAccumulatorType = ACCUMULATOR_TYPE_LOGICAL;
			iContextModelingOrderAccumulatorsWW = Accumulator::getContextModelingOrder( commandLineManager.getParameterValue("-ww"));
			if (commandLineManager.isParameterSet("-cw") == true) {
				iContextModelingOrderAccumulatorsCW = Accumulator::getContextModelingOrder( commandLineManager.getParameterValue("-cw"));
			} else {
				iContextModelingOrderAccumulatorsCW = iContextModelingOrderAccumulatorsWW;
			}
		}
		const char *strFileOptionalSymbols = commandLineManager.getParameterValue("-opt");
		bool bMultiplePronunciations = CommandLineManager::str2bool(commandLineManager.getParameterValue("-pro"));
		float fForwardPruningBeam = atof(commandLineManager.getParameterValue("-fwd")); 
		float fBackwardPruningBeam = atof(commandLineManager.getParameterValue("-bwd")); 
		int iTrellisMaxSize = atoi(commandLineManager.getParameterValue("-tre")); 
		bool bTrellisCache = true;
		int iTrellisCacheMaxSize = atoi(commandLineManager.getParameterValue("-tre"));	
		const char *strFileAccumulators = commandLineManager.getParameterValue("-dAcc");
			
		// create the accumulator object
		MLAccumulator mlAccumulator(strFilePhoneSet,strFileFeatureConfigurationAlignment,
			strFolderFeaturesAlignment,strFileModelsAlignment,iAccumulatorType,strFileOptionalSymbols,
			bMultiplePronunciations,iContextModelingOrderAccumulatorsWW,iContextModelingOrderAccumulatorsCW,
			strFileFeatureConfigurationAcc,strFolderFeaturesAcc,iCovarianceModellingAcc,strFileLexicon,strFileMLF,
			strFileAccumulators,fForwardPruningBeam,fBackwardPruningBeam,iTrellisMaxSize,
			bTrellisCache,iTrellisCacheMaxSize);
		
		mlAccumulator.initialize();
		mlAccumulator.accumulate();
		
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
	
	return 0;
}

