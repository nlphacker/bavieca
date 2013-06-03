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
#include "FeatureExtractor.h"
#include "HMMManager.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "VTLEstimator.h"

using namespace Bavieca;
 
// main function
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("vtlestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);		
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","lexicon",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-bat","batch file containing pairs (rawFile stateAlignmentFile)",
			PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");	
		commandLineManager.defineParameter("-out","file to store pairs (rawFile warpFactor)",PARAMETER_TYPE_FILE);
		commandLineManager.defineParameter("-fil","list of filler phones that will be ignored",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-floor","warp factor floor",
			PARAMETER_TYPE_FLOAT,true,"[0.75|0.99]","0.80");
		commandLineManager.defineParameter("-ceiling","warp factor ceiling",
			PARAMETER_TYPE_FLOAT,true,"[1.01|1.25]","1.20");
		commandLineManager.defineParameter("-step","warp factor increment between tests",
			PARAMETER_TYPE_FLOAT,true,"[0.01|0.05]","0.02");
		commandLineManager.defineParameter("-ali","whether to realign data for each warp factor",
			PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");	
		commandLineManager.defineParameter("-nrm","cepstral normalization mode",
			PARAMETER_TYPE_STRING,true,"none|utterance|session","utterance");	
		commandLineManager.defineParameter("-met","cepstral normalization method",
			PARAMETER_TYPE_STRING,true,"none|CMN|CMVN","CMN");	
		
		// parse the command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// load the configuration file
		const char *strFileConfiguration = commandLineManager.getParameterValue("-cfg");
		ConfigurationFeatures configurationFeatures(strFileConfiguration);
		configurationFeatures.load();
		
		// load the phone set
		PhoneSet phoneSet(commandLineManager.getParameterValue("-pho"));
		phoneSet.load();
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(commandLineManager.getParameterValue("-mod"));
		hmmManager.initializeDecoding();	
		
		// load the lexicon file
		LexiconManager lexiconManager(commandLineManager.getParameterValue("-lex"),&phoneSet);
		lexiconManager.load();
		
		// get command line parameters
		const char *strAlignmentFormat = commandLineManager.getParameterValue("-for");
		unsigned char iAlignmentFormat;
		if (strcmp(strAlignmentFormat,"binary") == 0) {
			iAlignmentFormat = FILE_FORMAT_BINARY;
		} else {
			assert(strcmp(strAlignmentFormat,"text") == 0);
			iAlignmentFormat = FILE_FORMAT_TEXT;
		}
		float fWarpFactorFloor = commandLineManager.getFloatParameterValue("-floor");
		float fWarpFactorCeiling = commandLineManager.getFloatParameterValue("-ceiling");
		float fWarpFactorStep = commandLineManager.getFloatParameterValue("-step");
		const char *strFileFillerPhones = NULL;
		if (commandLineManager.isParameterSet("-fil")) {
			strFileFillerPhones = commandLineManager.getParameterValue("-fil");
		}
		bool bRealignData = commandLineManager.getBoolParameterValue("-ali");
		int iCepstralNormalizationMode = 
			FeatureExtractor::getNormalizationMode(commandLineManager.getParameterValue("-nrm"));	
		int iCepstralNormalizationMethod = 
			FeatureExtractor::getNormalizationMethod(commandLineManager.getParameterValue("-met"));
		
		// get the output file (it is optional)
		const char *strOutputFile = NULL;
		if (commandLineManager.isParameterSet("-out")) {
			strOutputFile = commandLineManager.getParameterValue("-out");
		}
		
		// estimate the warp factor
		VTLEstimator vtlnManager(&configurationFeatures,&hmmManager,&phoneSet,
			&lexiconManager,strFileFillerPhones,fWarpFactorFloor,fWarpFactorCeiling,fWarpFactorStep,
			bRealignData,iCepstralNormalizationMode,iCepstralNormalizationMethod);
		
		float fLikelihoodGain = -FLT_MAX;
		float fWarpFactor;
		vtlnManager.estimateWarpFactor(commandLineManager.getParameterValue("-bat"),
			iAlignmentFormat,strOutputFile,fLikelihoodGain,fWarpFactor);
		
	} catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
 
	return 0; 
}
