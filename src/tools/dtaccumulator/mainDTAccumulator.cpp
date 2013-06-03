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
#include "HMMManager.h"
#include "LexiconManager.h"
#include "DTAccumulator.h"
#include "PhoneSet.h"

using namespace Bavieca;
 
// main for the tool "dtaccumulator"
int main(int argc, char *argv[]) {
	
	try {
		
		// define command line parameters
		CommandLineManager commandLineManager("mlaccumulator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation lexicon",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-fea","feature folder",PARAMETER_TYPE_FOLDER,false);
		commandLineManager.defineParameter("-cfg","feature configuration file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mlf","master label file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lat","directory containing hmm-marked lattices",PARAMETER_TYPE_FOLDER,false);
		commandLineManager.defineParameter("-ams","acoustic model scale factor",PARAMETER_TYPE_FLOAT,false);
		commandLineManager.defineParameter("-lms","language model scale factor",PARAMETER_TYPE_FLOAT,true,NULL,"1.0");	
		commandLineManager.defineParameter("-opt","file containing optional symbols that can be inserted at word boundaries",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-pro","allow multiple pronunciations",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
		commandLineManager.defineParameter("-fwd","forward pruning",PARAMETER_TYPE_FLOAT,true,"[-100.0|-10.0]","-20");
		commandLineManager.defineParameter("-bwd","backward pruning",PARAMETER_TYPE_FLOAT,true,"[100.0|10000.0]","800");
		commandLineManager.defineParameter("-tre","maximum trellis size (MB)",PARAMETER_TYPE_INTEGER,true,NULL,"500");
		commandLineManager.defineParameter("-dAccNum","file to dump accumulators (numerator)",
			PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-dAccDen","file to dump accumulators (denominator)",
			PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-obj","objective function",PARAMETER_TYPE_STRING,true,"MMI|bMMI","MMI");	
		commandLineManager.defineParameter("-bst","boosting factor for bMMI",PARAMETER_TYPE_FLOAT,true,NULL,"0.5");
		commandLineManager.defineParameter("-can","statistics cancelation",PARAMETER_TYPE_BOOLEAN,true,NULL,"yes");	
		
		// parse the command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// load the parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");
		const char *strFileLexicon = commandLineManager.getParameterValue("-lex");
		const char *strFolderFeatures = commandLineManager.getParameterValue("-fea"); 
		const char *strFileFeatureConfiguration = commandLineManager.getParameterValue("-cfg"); 
		const char *strFileMLF = commandLineManager.getParameterValue("-mlf"); 
		const char *strFolderLattices = commandLineManager.getParameterValue("-lat");
		float fScaleAM = atof(commandLineManager.getParameterValue("-ams"));
		float fScaleLM = atof(commandLineManager.getParameterValue("-lms"));	 
		const char *strFileOptionalSymbols = commandLineManager.getParameterValue("-opt");
		bool bMultiplePronunciations = CommandLineManager::str2bool(commandLineManager.getParameterValue("-pro"));
		float fForwardPruningBeam = atof(commandLineManager.getParameterValue("-fwd")); 
		float fBackwardPruningBeam = atof(commandLineManager.getParameterValue("-bwd")); 
		int iTrellisMaxSize = atoi(commandLineManager.getParameterValue("-tre")); 
		bool bTrellisCache = true;
		int iTrellisCacheMaxSize = atoi(commandLineManager.getParameterValue("-tre"));	
		const char *strFileAccumulatorsNum = commandLineManager.getParameterValue("-dAccNum");
		const char *strFileAccumulatorsDen = commandLineManager.getParameterValue("-dAccDen");
		const char *strObjectiveFunction = commandLineManager.getParameterValue("-obj");
		float fBoostingFactor = atof(commandLineManager.getParameterValue("-bst"));
		bool bStatisticsCancelation = CommandLineManager::str2bool(commandLineManager.getParameterValue("-can"));
			
		// create the accumulator object
		DTAccumulator dtAccumulator(strFilePhoneSet,strFileFeatureConfiguration,
			strFolderFeatures,strFileModels,strFileOptionalSymbols,bMultiplePronunciations,strFileLexicon,strFileMLF,
			strFolderLattices,fScaleAM,fScaleLM,strFileAccumulatorsNum,strFileAccumulatorsDen,strObjectiveFunction,
			fBoostingFactor,bStatisticsCancelation,fForwardPruningBeam,fBackwardPruningBeam,iTrellisMaxSize,bTrellisCache,
			iTrellisCacheMaxSize);
				
		dtAccumulator.initialize();
		dtAccumulator.accumulate();
	
	} catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;
}

