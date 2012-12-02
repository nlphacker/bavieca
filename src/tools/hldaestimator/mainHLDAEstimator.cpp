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
#include "HLDAEstimator.h"
#include "HMMManager.h"
#include "PhoneSet.h"

using namespace Bavieca;
 
// main for the tool "hldaestimator"
int main(int argc, char *argv[]) {

	// define the parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("hldaestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-acc","input accumulator filelist",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-itt","number of iterations for transform update",
		PARAMETER_TYPE_INTEGER,true,"[1|100]","10");	
	m_commandLineManager->defineParameter("-itp","number of iterations for parameter update",
		PARAMETER_TYPE_INTEGER,true,"[1|100]","10");	
	m_commandLineManager->defineParameter("-red","dimensionality reduction",PARAMETER_TYPE_INTEGER,true,NULL,"13");	
	m_commandLineManager->defineParameter("-out","output folder",PARAMETER_TYPE_FOLDER,false);	
	
	// parse the parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// retrieve the parameters
	const char *m_strFileFeatureConfiguration = m_commandLineManager->getParameterValue("-cfg");
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *m_strFileAccList = m_commandLineManager->getParameterValue("-acc");
	int m_iIterationsTransformUpdate = atoi(m_commandLineManager->getParameterValue("-itt"));
	int m_iIterationsParameterUpdate = atoi(m_commandLineManager->getParameterValue("-itp"));
	int m_iDimensionalityReduction = atoi(m_commandLineManager->getParameterValue("-red"));
	const char *m_strFolderOutput = m_commandLineManager->getParameterValue("-out");
	
   // load the feature configuration
   ConfigurationFeatures *configurationFeatures = new ConfigurationFeatures(m_strFileFeatureConfiguration);
   configurationFeatures->load();
   
   delete configurationFeatures;	
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
 
	// load the HMM-models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);	
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeEstimation(ACCUMULATOR_TYPE_PHYSICAL,
		m_hmmManager->getContextModelingOrderHMM(),m_hmmManager->getContextModelingOrderHMMCW());	

	// initialize the HLDA estimator
	HLDAEstimator *m_hldaEstimator = new HLDAEstimator(m_hmmManager,m_strFileAccList,m_iDimensionalityReduction,m_iIterationsTransformUpdate,m_iIterationsParameterUpdate,m_strFolderOutput);
	
	// do the actual estimation
	m_hldaEstimator->estimate();
	
	// clean-up
	delete m_hldaEstimator;
	delete m_hmmManager;
	delete m_phoneSet;

	return 0;
}

