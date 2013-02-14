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
#include "ConfigurationFeatures.h"
#include "FeatureExtractor.h"
#include "Global.h"
#include "TimeUtils.h"

using namespace std;

#include <string>
#include <map>

using namespace Bavieca;

// main for the feature extraction tool: "param"
int main(int argc, char *argv[]) {

	try {

		// define command line parameters
		CommandLineManager commandLineManager("param",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-raw","audio file to parameterize (raw format)",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-fea","file to store the feature vectors extracted",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-bat","batch file with pairs [rawFile featureFile]",
			PARAMETER_TYPE_FILE,true);	
		commandLineManager.defineParameter("-wrp","frequency warp factor",PARAMETER_TYPE_FLOAT,true,"[0.80|1.20]","1.0");	
		commandLineManager.defineParameter("-nrm","cepstral normalization mode",
			PARAMETER_TYPE_STRING,true,"none|utterance|session","utterance");	
		commandLineManager.defineParameter("-met","cepstral normalization method",PARAMETER_TYPE_STRING,true,"none|CMN|CMVN","CMN");	
		commandLineManager.defineParameter("-hlt","whether to halt the batch processing if an error is found",
			PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");	
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the parameters
		
		// load the configuration file
		const char *strFileConfiguration = commandLineManager.getParameterValue("-cfg");
		ConfigurationFeatures configurationFeatures(strFileConfiguration);
		configurationFeatures.load();
		
		// perform the parameterization
			
		// load the parameters
		float fWarpFactor = commandLineManager.getFloatParameterValue("-wrp"); 
		int iCepstralBufferSize = -1;
		int iCepstralNormalizationMode = 
			FeatureExtractor::getNormalizationMode(commandLineManager.getParameterValue("-nrm"));	
		int iCepstralNormalizationMethod = 
			FeatureExtractor::getNormalizationMethod(commandLineManager.getParameterValue("-met"));
				
		FeatureExtractor featureExtractor(&configurationFeatures,fWarpFactor,iCepstralBufferSize,
			iCepstralNormalizationMode,iCepstralNormalizationMethod);
		
		featureExtractor.initialize();
		if (commandLineManager.isParameterSet("-bat") == false) {
			featureExtractor.extractFeatures(commandLineManager.getParameterValue("-raw"),
				commandLineManager.getParameterValue("-fea"));
		} else {
			if (featureExtractor.extractFeaturesBatch(commandLineManager.getParameterValue("-bat"),
				CommandLineManager::str2bool(commandLineManager.getParameterValue("-hlt"))) == false) {
				BVC_ERROR << "problem processing the batch file";	
			}
		}
		
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}
	
	return 0;
}
