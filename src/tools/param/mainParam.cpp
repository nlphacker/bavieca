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
//int mainParam(int argc, char *argv[]) {
int main(int argc, char *argv[]) {

	// define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("param",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-raw","audio file to parameterize (raw format)",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-fea","file to store the feature vectors extracted",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-bat","batch file with pairs [rawFile featureFile]",
		PARAMETER_TYPE_FILE,true);	
	m_commandLineManager->defineParameter("-wrp","frequency warp factor",PARAMETER_TYPE_FLOAT,true,"[0.80|1.20]","1.0");	
	m_commandLineManager->defineParameter("-nrm","cepstral normalization mode",
		PARAMETER_TYPE_STRING,true,"none|utterance|session","utterance");	
	m_commandLineManager->defineParameter("-met","cepstral normalization method",PARAMETER_TYPE_STRING,true,"none|CMN|CMVN","CMN");	
	m_commandLineManager->defineParameter("-hlt","whether to halt the batch processing if an error is found",
		PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");	
	
	// parse the parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get the parameters
	
	// load the configuration file
	const char *m_strFileConfiguration = m_commandLineManager->getParameterValue("-cfg");
	ConfigurationFeatures *m_configurationFeatures = new ConfigurationFeatures(m_strFileConfiguration);
	m_configurationFeatures->load();
	
	// perform the parameterization
		
	// load the parameters
	float m_fWarpFactor = m_commandLineManager->getFloatParameterValue("-wrp"); 
	int m_iCepstralBufferSize = -1;
	int m_iCepstralNormalizationMode = 
		FeatureExtractor::getNormalizationMode(m_commandLineManager->getParameterValue("-nrm"));	
	int m_iCepstralNormalizationMethod = 
		FeatureExtractor::getNormalizationMethod(m_commandLineManager->getParameterValue("-met"));
			
	FeatureExtractor *m_featureExtractor = new FeatureExtractor(m_configurationFeatures,m_fWarpFactor,m_iCepstralBufferSize,m_iCepstralNormalizationMode,m_iCepstralNormalizationMethod);
	
	m_featureExtractor->initialize();
	if (m_commandLineManager->isParameterSet("-bat") == false) {
		m_featureExtractor->extractFeatures(m_commandLineManager->getParameterValue("-raw"),
			m_commandLineManager->getParameterValue("-fea"));
	} else {
		if (m_featureExtractor->extractFeaturesBatch(m_commandLineManager->getParameterValue("-bat"),
			CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-hlt"))) == false) {
			printf("Error processing the batch file\n");	
		}
	}
	delete m_featureExtractor;
	
	// clean-up
	delete m_commandLineManager;
	delete m_configurationFeatures;
	
	return 0;
}
