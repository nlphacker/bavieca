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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "CommandLineManager.h"
#include "HMMManager.h"
#include "HLDAEstimator.h"
#include "PhoneSet.h"
#include "RegressionTree.h"
#include "Transform.h"

using namespace std;

#include <string>
#include <map>

using namespace Bavieca;

// main for the tool "hmmx"
int main(int argc, char *argv[])
{	
	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("hmmx",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-tra","model transform",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-rgt","regression-tree",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-in","inout acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
	
	// (2) process command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// load parameters
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileTransform = m_commandLineManager->getParameterValue("-tra");
	const char *m_strFileRegressionTree = m_commandLineManager->getParameterValue("-rgt");
	const char *m_strFileModelsInput = m_commandLineManager->getParameterValue("-in");
	const char *m_strFileModelsOutput = m_commandLineManager->getParameterValue("-out");
	
	// load the phone set
	PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
	m_phoneSet->load();
	
	HMMManager *m_hmmManager = NULL;
		
	// regression tree based transform?	
	if (m_strFileRegressionTree != NULL) {
	
		// load the acoustic models
		m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
		m_hmmManager->load(m_strFileModelsInput);
		m_hmmManager->initializeDecoding();	
			
		// load the regression tree
		RegressionTree *m_regressionTree = new RegressionTree(m_hmmManager);
		m_regressionTree->load(m_strFileRegressionTree);
		
		// load the transforms
		m_regressionTree->loadTransforms(m_strFileTransform);
		
		// apply the transforms
		m_regressionTree->applyTransforms();
		
		delete m_regressionTree;
	}
	// regular model transform
	else {
	
		// load the input acoustic models
		m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);
		m_hmmManager->load(m_strFileModelsInput);
		
		// load the transformation
		Transform *m_transform = new Transform();
		m_transform->load(m_strFileTransform);
		
		// apply the actual transformation
		HLDAEstimator::applyTransform(m_transform->getTransform(),m_hmmManager);
	
		delete m_transform;	
	}
	
	// store the output acoustic models to disk
	m_hmmManager->store(m_strFileModelsOutput);
	
	// clean-up	
	delete m_hmmManager;
	delete m_phoneSet;
	delete m_commandLineManager;	
	
	return 0;
}



