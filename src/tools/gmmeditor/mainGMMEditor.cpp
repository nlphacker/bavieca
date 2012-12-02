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
#include "ConfigurationFeatures.h"
#include "FeatureExtractor.h"
#include "GMMEditor.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "HMMInitializer.h"
#include "HMMManager.h"
#include "MLEstimator.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the HMM edition tool: "hmmeditor"
int main(int argc, char *argv[]) {

	// define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("hmmeditor",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-acc","accumulators",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-dbl","mixture doubling",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
	m_commandLineManager->defineParameter("-inc","mixture increment",PARAMETER_TYPE_INTEGER,true);
	m_commandLineManager->defineParameter("-crt","splitting criterion",PARAMETER_TYPE_STRING,true,"covariance|weight","covariance");
	m_commandLineManager->defineParameter("-occ","minimum component occupation",PARAMETER_TYPE_FLOAT,true,NULL,"100.0");
	m_commandLineManager->defineParameter("-wgh","minimum component weight",PARAMETER_TYPE_FLOAT,true,NULL,"0.00001");
	m_commandLineManager->defineParameter("-eps","epsilon (mean perturbation)",PARAMETER_TYPE_FLOAT,true,NULL,"0.05");
	m_commandLineManager->defineParameter("-mrg","merge components below minimum occ/wgh",PARAMETER_TYPE_BOOLEAN,true,NULL,"yes");	
	m_commandLineManager->defineParameter("-cov","covariance flooring ratio",PARAMETER_TYPE_FLOAT,true,NULL,"0.05");	
	m_commandLineManager->defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-vrb","verbose",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");	
	
	// parse command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
		
	// get the parameters
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *m_strFileAccList = m_commandLineManager->getParameterValue("-acc");	
	bool m_bMixtureDoubling = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-dbl"));
	int m_iMixtureIncrement = 0;
	if (m_commandLineManager->isParameterSet("-inc") == true) {	
		m_iMixtureIncrement = atoi(m_commandLineManager->getParameterValue("-inc"));
	}
	const char *m_strSplittingCriterion = m_commandLineManager->getParameterValue("-crt");
	int m_iSplittingCriterion = GMMEditor::getSplittingCriterion(m_strSplittingCriterion);
	float m_fMinimumComponentOccupation = atof(m_commandLineManager->getParameterValue("-occ"));
	float m_fMinimumComponentWeight = atof(m_commandLineManager->getParameterValue("-wgh"));
	float m_fEpsilon = atof(m_commandLineManager->getParameterValue("-eps"));
	bool m_bMixtureMerging = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-mrg"));
	float m_fCovarianceFlooringRatio = atof(m_commandLineManager->getStrParameterValue("-cov"));
	const char *m_strFileModelsOutput = m_commandLineManager->getParameterValue("-out");
	//bool m_bVerbose = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-vrb"));
	
	// load the phonetic symbol set
	PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
	m_phoneSet->load();
	
	// load the HMMs
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);	
	m_hmmManager->load(m_strFileModels);
	
	// load the accumulators
	AccMetadata metadata;
	MAccumulatorPhysical mAccumulatorPhysical;
	Accumulator::loadAccumulatorList(m_strFileAccList,mAccumulatorPhysical,metadata);
	
	// create the GMM editor
	GMMEditor *m_gmmEditor = new GMMEditor(m_hmmManager,&mAccumulatorPhysical);
	m_gmmEditor->initialize();	
	
	// mixture merging ?
	if (m_bMixtureMerging == true) {		
		m_gmmEditor->mixtureMerge(m_fMinimumComponentOccupation,m_fMinimumComponentWeight,m_fCovarianceFlooringRatio);
	}
	
	// mixture doubling
	if (m_bMixtureDoubling == true) {
		m_gmmEditor->mixtureDouble(m_fMinimumComponentOccupation,m_fEpsilon);
	}
	// mixture increment
	else {
		m_gmmEditor->mixtureIncrement(m_iMixtureIncrement,m_iSplittingCriterion,false,
			m_fMinimumComponentOccupation,m_fEpsilon);
	}
		
	// store the models after the Gaussian merging process is done	
	m_hmmManager->store(m_strFileModelsOutput);
		
	// clean-up
	//delete m_hmmManager;
	delete m_phoneSet;
	delete m_commandLineManager;	
	
	return -1;	
}
