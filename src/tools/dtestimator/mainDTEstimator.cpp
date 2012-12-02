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
#include "HMMManager.h"
#include "LexiconManager.h"
#include "DTEstimator.h"
#include "MLEstimator.h"
#include "MLFFile.h"
#include "PhoneSet.h"

using namespace Bavieca;
 
// main for the tool "DTEstimator"
int main(int argc, char *argv[]) {
	
	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("dtestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-accNum","input accumulator filelist (numerator)",
		PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-accDen","input accumulator filelist (denominator)",
		PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-cov","covariance flooring ratio",PARAMETER_TYPE_FLOAT,true,NULL,"0.05");
	m_commandLineManager->defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-E","learning rate constant",PARAMETER_TYPE_FILE,true,NULL,"2.0");
	m_commandLineManager->defineParameter("-I","I-smoothing",PARAMETER_TYPE_STRING,true,"none|prev","none");
	m_commandLineManager->defineParameter("-tau","I-smoothing constant",PARAMETER_TYPE_FLOAT,true,NULL,"100.0");
	
	// parse the command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get the parameters
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *m_strFileAccListNum = m_commandLineManager->getParameterValue("-accNum");
	const char *m_strFileAccListDen = m_commandLineManager->getParameterValue("-accDen");
	float m_fCovarianceFlooringRatio = atof(m_commandLineManager->getParameterValue("-cov"));
	const char *m_strFileModelsOutput = m_commandLineManager->getParameterValue("-out");
	float m_fE = atof(m_commandLineManager->getParameterValue("-E"));
	const char *m_strISmoothing = m_commandLineManager->getParameterValue("-I");
	float m_fTau = -1.0;
	if (m_commandLineManager->isParameterSet("-tau")) {
		m_fTau = atof(m_commandLineManager->getParameterValue("-tau"));
	}
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
   	
	// load the HMMs
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);	
	m_hmmManager->load(m_strFileModels);

	// initialize the HMMs
	m_hmmManager->initializeEstimation(ACCUMULATOR_TYPE_PHYSICAL,UCHAR_MAX,UCHAR_MAX);
	
	// create the discriminative training estimator
	DTEstimator *m_dtEstimator = new DTEstimator(m_hmmManager);
	
	// estimate parameters
	m_dtEstimator->estimateParameters(m_strFileAccListNum,m_strFileAccListDen,m_fE,m_strISmoothing,m_fTau,true);
	
	// floor covariances
	m_dtEstimator->floorCovariances(m_fCovarianceFlooringRatio);
	
	// create the output HMMs
	m_hmmManager->store(m_strFileModelsOutput);
	
	// clean-up 
	delete m_dtEstimator;
	delete m_hmmManager;
	delete m_phoneSet;
	delete m_commandLineManager;

	return 0;
}
