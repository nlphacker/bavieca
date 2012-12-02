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

	// define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("mlaccumulator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-fea","feature folder",PARAMETER_TYPE_FOLDER,false);
	m_commandLineManager->defineParameter("-cfg","feature configuration file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-feaA","feature folder (stats accumulation)",PARAMETER_TYPE_FOLDER,true);	
	m_commandLineManager->defineParameter("-cfgA","feature configuration (stats accumulation)",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-covA","covariance modeling type (stats accumulation)",
		PARAMETER_TYPE_STRING,true,"diagonal|full","diagonal");
	m_commandLineManager->defineParameter("-mlf","master label file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-ww","within-word context modeling order for the accumulators",PARAMETER_TYPE_STRING,true);
	m_commandLineManager->defineParameter("-cw","cross-word context modeling order for the accumulators",PARAMETER_TYPE_STRING,true);
	m_commandLineManager->defineParameter("-opt","file containing optional symbols that can be inserted at word boundaries",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-pro","allow multiple pronunciations",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
	m_commandLineManager->defineParameter("-fwd","forward pruning",PARAMETER_TYPE_FLOAT,true,"[-100.0|-10.0]","-20");
	m_commandLineManager->defineParameter("-bwd","backward pruning",PARAMETER_TYPE_FLOAT,true,"[100.0|10000.0]","800");
	m_commandLineManager->defineParameter("-tre","maximum trellis size (MB)",PARAMETER_TYPE_INTEGER,true,NULL,"500");
	m_commandLineManager->defineParameter("-dAcc","file to dump accumulators",PARAMETER_TYPE_FILE,false);	
	
	// parse the command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
		
	// load the parameters
	const char *strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *strFileModelsAlignment = m_commandLineManager->getParameterValue("-mod");
	const char *strFileLexicon = m_commandLineManager->getParameterValue("-lex");
	const char *strFolderFeaturesAlignment = m_commandLineManager->getParameterValue("-fea"); 
	const char *strFileFeatureConfigurationAlignment = m_commandLineManager->getParameterValue("-cfg"); 
	const char *strFolderFeaturesAcc = m_commandLineManager->getParameterValue("-feaA"); 
	const char *strFileFeatureConfigurationAcc = m_commandLineManager->getParameterValue("-cfgA"); 
	int iCovarianceModellingAcc;
	if (m_commandLineManager->isParameterSet("-covA") == false) {
		iCovarianceModellingAcc = COVARIANCE_MODELLING_TYPE_DEFAULT;
	} else if (strcmp(m_commandLineManager->getParameterValue("-covA"),COVARIANCE_MODELLING_TYPE_FULL_STR) == 0) {
		iCovarianceModellingAcc = COVARIANCE_MODELLING_TYPE_FULL;
	} else {
		assert(strcmp(m_commandLineManager->getParameterValue("-covA"),COVARIANCE_MODELLING_TYPE_DIAGONAL_STR) == 0);
		iCovarianceModellingAcc = COVARIANCE_MODELLING_TYPE_DIAGONAL;
	}
	const char *strFileMLF = m_commandLineManager->getParameterValue("-mlf"); 
	unsigned char iAccumulatorType = ACCUMULATOR_TYPE_PHYSICAL;
	int iContextModelingOrderAccumulatorsWW = UCHAR_MAX;
	int iContextModelingOrderAccumulatorsCW = UCHAR_MAX;
	if (m_commandLineManager->isParameterSet("-ww") == true) {
		iAccumulatorType = ACCUMULATOR_TYPE_LOGICAL;
		iContextModelingOrderAccumulatorsWW = Accumulator::getContextModelingOrder( m_commandLineManager->getParameterValue("-ww"));
		if (m_commandLineManager->isParameterSet("-cw") == true) {
			iContextModelingOrderAccumulatorsCW = Accumulator::getContextModelingOrder( m_commandLineManager->getParameterValue("-cw"));
		} else {
			iContextModelingOrderAccumulatorsCW = iContextModelingOrderAccumulatorsWW;
		}
	}
	const char *strFileOptionalSymbols = m_commandLineManager->getParameterValue("-opt");
	bool bMultiplePronunciations = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-pro"));
	float fForwardPruningBeam = atof(m_commandLineManager->getParameterValue("-fwd")); 
	float fBackwardPruningBeam = atof(m_commandLineManager->getParameterValue("-bwd")); 
	int iTrellisMaxSize = atoi(m_commandLineManager->getParameterValue("-tre")); 
	bool bTrellisCache = true;
	int iTrellisCacheMaxSize = atoi(m_commandLineManager->getParameterValue("-tre"));	
	const char *strFileAccumulators = m_commandLineManager->getParameterValue("-dAcc");
		
	// create the accumulator object
	MLAccumulator *m_mlAccumulator = new MLAccumulator(strFilePhoneSet,strFileFeatureConfigurationAlignment,
		strFolderFeaturesAlignment,strFileModelsAlignment,iAccumulatorType,strFileOptionalSymbols,
		bMultiplePronunciations,iContextModelingOrderAccumulatorsWW,iContextModelingOrderAccumulatorsCW,
		strFileFeatureConfigurationAcc,strFolderFeaturesAcc,iCovarianceModellingAcc,strFileLexicon,strFileMLF,
		strFileAccumulators,fForwardPruningBeam,fBackwardPruningBeam,iTrellisMaxSize,
		bTrellisCache,iTrellisCacheMaxSize);
	
	m_mlAccumulator->initialize();
	m_mlAccumulator->accumulate();
		
	// clean-up	
	delete m_mlAccumulator;
	delete m_commandLineManager;
	
	return 0;
}

