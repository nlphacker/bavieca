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
#include "RegressionTree.h"
#include "TimeUtils.h"

using namespace Bavieca;

// main for the tool "regtree"
int main(int argc, char *argv[]) {
	
	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("regtree",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-met","clustering method",PARAMETER_TYPE_STRING,true,"kMeans|EM","EM");
	m_commandLineManager->defineParameter("-rgc","number of regression classes (base-classes)",PARAMETER_TYPE_INTEGER,true,"[1|1000]","50");
	m_commandLineManager->defineParameter("-gau","minimum number of Gaussian components per base-class",PARAMETER_TYPE_INTEGER,true,NULL,"50");
	m_commandLineManager->defineParameter("-out","file to store the regression tree",PARAMETER_TYPE_FILE,false);
	
	// (2) parse command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
	// clustering method
	unsigned char m_iClusteringMethod;
	if (strcmp(m_commandLineManager->getParameterValue("-met"),"kMeans") == 0) {
		 m_iClusteringMethod = CLUSTERING_METHOD_KMEANS;
	} else {
		 m_iClusteringMethod = CLUSTERING_METHOD_EM;
	}	
	int m_iRegressionClasses = atoi(m_commandLineManager->getParameterValue("-rgc"));
	int m_iMinimumGaussianComponentsBaseClass = atoi(m_commandLineManager->getParameterValue("-gau"));
	const char *m_strFileRegressionTree = m_commandLineManager->getParameterValue("-out");

   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeDecoding();
	
	// build the regression tree
	RegressionTree *m_regressionTree = new RegressionTree(m_hmmManager);
	m_regressionTree->build(m_iRegressionClasses,m_iClusteringMethod,
		m_iMinimumGaussianComponentsBaseClass);
		
	// store the regression tree to disk
	m_regressionTree->store(m_strFileRegressionTree);
	
	// clean-up
	delete m_regressionTree;
	delete m_hmmManager;
	delete m_phoneSet;
	delete m_commandLineManager;
		
	return 0;
}
