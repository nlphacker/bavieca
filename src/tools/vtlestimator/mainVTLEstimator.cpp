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
#include "HMMManager.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "VTLEstimator.h"

using namespace Bavieca;
 
// main function
int main(int argc, char *argv[]) {

	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("vtlestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);		
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lex","lexicon",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-bat","batch file containing pairs (rawFile stateAlignmentFile)",
		PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");	
	m_commandLineManager->defineParameter("-out","file to store pairs (rawFile warpFactor)",PARAMETER_TYPE_FILE);
	m_commandLineManager->defineParameter("-fil","list of filler phones that will be ignored",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-floor","warp factor floor",
		PARAMETER_TYPE_FLOAT,true,"[0.75|0.99]","0.80");
	m_commandLineManager->defineParameter("-ceiling","warp factor ceiling",
		PARAMETER_TYPE_FLOAT,true,"[1.01|1.25]","1.20");
	m_commandLineManager->defineParameter("-step","warp factor increment between tests",
		PARAMETER_TYPE_FLOAT,true,"[0.01|0.05]","0.02");
	m_commandLineManager->defineParameter("-ali","whether to realign data for each warp factor",
		PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");	
	m_commandLineManager->defineParameter("-nrm","cepstral normalization mode",
		PARAMETER_TYPE_STRING,true,"none|utterance|session","utterance");	
	m_commandLineManager->defineParameter("-met","cepstral normalization method",
		PARAMETER_TYPE_STRING,true,"none|CMN|CMVN","CMN");	
	
	// parse the command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// load the configuration file
	const char *m_strFileConfiguration = m_commandLineManager->getParameterValue("-cfg");
	ConfigurationFeatures *m_configurationFeatures = new ConfigurationFeatures(m_strFileConfiguration);
	m_configurationFeatures->load();
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_commandLineManager->getParameterValue("-pho"));
   m_phoneSet->load();
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_commandLineManager->getParameterValue("-mod"));
	m_hmmManager->initializeDecoding();	
	
	// load the lexicon file
	LexiconManager *m_lexiconManager = new LexiconManager(m_commandLineManager->getParameterValue("-lex"),m_phoneSet);
	m_lexiconManager->load();
	
	// get command line parameters
	const char *m_strAlignmentFormat = m_commandLineManager->getParameterValue("-for");
	unsigned char m_iAlignmentFormat;
	if (strcmp(m_strAlignmentFormat,"binary") == 0) {
		m_iAlignmentFormat = FILE_FORMAT_BINARY;
	} else {
		assert(strcmp(m_strAlignmentFormat,"text") == 0);
		m_iAlignmentFormat = FILE_FORMAT_TEXT;
	}
	float m_fWarpFactorFloor = m_commandLineManager->getFloatParameterValue("-floor");
	float m_fWarpFactorCeiling = m_commandLineManager->getFloatParameterValue("-ceiling");
	float m_fWarpFactorStep = m_commandLineManager->getFloatParameterValue("-step");
	const char *m_strFileFillerPhones = NULL;
	if (m_commandLineManager->isParameterSet("-fil")) {
		m_strFileFillerPhones = m_commandLineManager->getParameterValue("-fil");
	}
	bool m_bRealignData = m_commandLineManager->getBoolParameterValue("-ali");
	int m_iCepstralNormalizationMode = 
		FeatureExtractor::getNormalizationMode(m_commandLineManager->getParameterValue("-nrm"));	
	int m_iCepstralNormalizationMethod = 
		FeatureExtractor::getNormalizationMethod(m_commandLineManager->getParameterValue("-met"));
	
	// get the output file (it is optional)
	const char *m_strOutputFile = NULL;
	if (m_commandLineManager->isParameterSet("-out")) {
		m_strOutputFile = m_commandLineManager->getParameterValue("-out");
	}
	
	// estimate the warp factor
	VTLEstimator *m_vtlnManager = new VTLEstimator(m_configurationFeatures,m_hmmManager,m_phoneSet,
		m_lexiconManager,m_strFileFillerPhones,m_fWarpFactorFloor,m_fWarpFactorCeiling,m_fWarpFactorStep,
		m_bRealignData,m_iCepstralNormalizationMode,m_iCepstralNormalizationMethod);
	
	float m_fLikelihoodGain = -FLT_MAX;
	float m_fWarpFactor;
	m_vtlnManager->estimateWarpFactor(m_commandLineManager->getParameterValue("-bat"),
		m_iAlignmentFormat,m_strOutputFile,m_fLikelihoodGain,m_fWarpFactor);
	
	// clean-up
	delete m_vtlnManager;
	delete m_lexiconManager;
	delete m_hmmManager;
	delete m_configurationFeatures;
	delete m_phoneSet;
	delete m_commandLineManager;
 
	return 0; 
}
