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
#include "MAPEstimator.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the HMM estimation tool: "mapestimator"
int main(int argc, char *argv[]) {

	// (1) define command line parameters  
	CommandLineManager *m_commandLineManager = new CommandLineManager("mapestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-acc","input accumulator filelist",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-pkw","prior knowledge weight",PARAMETER_TYPE_FLOAT,true,"[2|20]","2");	
	m_commandLineManager->defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
	
	// parse the command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get the parameters
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *m_strFileAccList = m_commandLineManager->getParameterValue("-acc");
	float m_fPriorKnowledgeWeight = atof(m_commandLineManager->getParameterValue("-pkw"));
	const char *m_strFileModelsOutput = m_commandLineManager->getParameterValue("-out");
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
   	
	// load the HMMs
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);	
	m_hmmManager->load(m_strFileModels);
	
	// load the accumulators
	AccMetadata metadata;
	MAccumulatorPhysical mAccumulatorPhysical;
	Accumulator::loadAccumulatorList(m_strFileAccList,mAccumulatorPhysical,metadata);
	 
	// initialize the HMMs 
	m_hmmManager->initializeEstimation(ACCUMULATOR_TYPE_PHYSICAL,UCHAR_MAX,UCHAR_MAX);
	 
	HMMManager *m_hmmManagerUpdate = m_hmmManager;
	
	// estimate the HMM parameters
	MAPEstimator *m_mapEstimator = new MAPEstimator(m_hmmManagerUpdate);
	m_mapEstimator->estimateParameters(mAccumulatorPhysical,m_fPriorKnowledgeWeight);
	
	Accumulator::destroy(mAccumulatorPhysical);
	
	// create the output HMMs
	m_hmmManagerUpdate->store(m_strFileModelsOutput);
	
	// clean-up
	delete m_mapEstimator;
	delete m_hmmManager; 
	delete m_phoneSet;
	delete m_commandLineManager;
	
	return 0;
}
