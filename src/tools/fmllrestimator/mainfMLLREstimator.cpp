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
#include "FMLLREstimator.h"
#include "PhoneSet.h"
#include "TimeUtils.h"
#include "Transform.h"

using namespace Bavieca;

// main for the tool "mapEstimator"
int main(int argc, char *argv[]) {
	
	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("fmllrestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","bootstrap acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-bat","batch file containing pairs (featureFile alignmentFile)",
		PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");
	m_commandLineManager->defineParameter("-tra","feature transform",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-bst","whether to assign all occupation to best scoring Gaussian component",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");
	
	// parse the parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// load the parameter values
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *m_strFileBatch = m_commandLineManager->getParameterValue("-bat");
	const char *m_strFormat = m_commandLineManager->getParameterValue("-for");
	const char *m_strFileTransform = m_commandLineManager->getParameterValue("-tra");
	bool m_bBestComponentOnly = m_commandLineManager->getBoolParameterValue("-bst");
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeDecoding();	
	
	// get starting time
   double dBegin = TimeUtils::getTimeMilliseconds();	
		
	// create the fMLLR estimator
	FMLLREstimator *m_fMLLREstimator = new FMLLREstimator(m_phoneSet,m_hmmManager,20,m_bBestComponentOnly);
	m_fMLLREstimator->initializeEstimation();
	double dLikelihood = 0.0;
	m_fMLLREstimator->feedAdaptationData(m_strFileBatch,m_strFormat,&dLikelihood,true);	
	
	double dEndFeed = TimeUtils::getTimeMilliseconds();
	
	Transform *transform = m_fMLLREstimator->estimateTransform(NULL);
	assert(transform);
	
	// store the transform
	transform->store(m_strFileTransform);
	delete transform;
	
	m_fMLLREstimator->uninitializeEstimation();
	
	double dEndComputation = TimeUtils::getTimeMilliseconds();
	
   // show the processing time
   printf("collecting data:      %6.2f seconds\n",(dEndFeed-dBegin)/1000.0);
   printf("computing transforms: %6.2f seconds\n",(dEndComputation-dEndFeed)/1000.0);
	
	// clean-up
	delete m_phoneSet;
	delete m_hmmManager;
	delete m_fMLLREstimator;
	delete m_commandLineManager;
	
	return 0;
}
