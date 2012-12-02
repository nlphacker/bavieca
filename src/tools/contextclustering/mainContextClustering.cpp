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
#include "ContextModeling.h"
#include "Global.h"
#include "PhoneSet.h"
#include "PhoneticRulesManager.h"
#include "TimeUtils.h"

using namespace std;

#include <string>
#include <map>

using namespace Bavieca;

// main for the "contextclustering" tool
int main(int argc, char *argv[]) {

	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("contextclustering",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-rul","phonetic rules",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-ww","within-word context modeling order",PARAMETER_TYPE_STRING,true,NULL,HMM_CONTEXT_MODELING_TRIPHONES_STR);
	m_commandLineManager->defineParameter("-cw","cross-word context modeling order",PARAMETER_TYPE_STRING,true,NULL,HMM_CONTEXT_MODELING_TRIPHONES_STR);
	m_commandLineManager->defineParameter("-met","clustering method",PARAMETER_TYPE_STRING,true,"local|global",CLUSTERING_METHOD_LOCAL_STR);
	m_commandLineManager->defineParameter("-mrg","whether to perform bottom up merging",PARAMETER_TYPE_BOOLEAN,true,NULL,"yes");
	m_commandLineManager->defineParameter("-acc","logical accumulator list",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-occ","minimum cluster occupation",PARAMETER_TYPE_FLOAT,true,NULL,"200");
	m_commandLineManager->defineParameter("-gan","minimum likelihood gain",PARAMETER_TYPE_FLOAT,true,NULL,"2000");
	m_commandLineManager->defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
	
	// parse the parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}

	// get the parameter values	
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");	
	const char *m_strFilePhoneticRules = m_commandLineManager->getParameterValue("-rul");
	unsigned char m_iContextModelingOrderWW = Accumulator::getContextModelingOrder(m_commandLineManager->getParameterValue("-ww"));
	if (m_iContextModelingOrderWW == UCHAR_MAX) {
		return -1;
	}
	unsigned char m_iContextModelingOrderCW = Accumulator::getContextModelingOrder(m_commandLineManager->getParameterValue("-cw"));
	if (m_iContextModelingOrderCW == UCHAR_MAX) {
		return -1;
	}
	bool m_bGlobalClustering = false;
	if (strcmp(m_commandLineManager->getParameterValue("-met"),"global") == 0) {
		m_bGlobalClustering = true;
	}	
	bool m_bBottomUpMerging = true;
	if (strcmp(m_commandLineManager->getParameterValue("-mrg"),"no") == 0) {
		m_bBottomUpMerging = false;
	}	
	const char *m_strFileAccList = m_commandLineManager->getParameterValue("-acc");
	float m_fMinimumClusterOccupation = atof(m_commandLineManager->getParameterValue("-occ"));
	float m_fMinimumLikelihoodGain = atof(m_commandLineManager->getParameterValue("-gan"));
	const char *m_strFileModelsOutput = m_commandLineManager->getParameterValue("-out");
	
	// load the phonetic symbol set
	PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
	m_phoneSet->load();
	
	// load the HMMs
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_ESTIMATION);	
	m_hmmManager->load(m_strFileModels);
	
	double dTimeStart = TimeUtils::getTimeMilliseconds();
	
	// load the accumulators
	AccMetadata metadata;
	MAccumulatorLogical mAccumulatorLogical;
	Accumulator::loadAccumulatorList(m_strFileAccList,mAccumulatorLogical,metadata);

	double dTimeEndLoading = TimeUtils::getTimeMilliseconds();
	
	// do the clustering
	ContextModeling *m_contextModeling = new ContextModeling(m_hmmManager->getFeatureDimensionality(),
		COVARIANCE_MODELLING_TYPE_DIAGONAL,
		m_phoneSet,m_iContextModelingOrderWW,m_iContextModelingOrderCW,m_strFilePhoneticRules,m_bGlobalClustering,
		m_bBottomUpMerging,m_fMinimumClusterOccupation,m_fMinimumLikelihoodGain);

	if (m_contextModeling->clusterContextDependentUnits(mAccumulatorLogical,m_hmmManager) == false) {
		BVC_ERROR << "problem found during context dependent clustering";
	}
	
	double dTimeEndClustering = TimeUtils::getTimeMilliseconds();
	
	printf("loading time:    %.2f seconds\n",(dTimeEndLoading-dTimeStart)/1000.0);
	printf("clustering time: %.2f seconds\n",(dTimeEndClustering-dTimeEndLoading)/1000.0);
	
	// create the output HMMs
	m_hmmManager->store(m_strFileModelsOutput);
	
	// clean-up
	Accumulator::destroy(mAccumulatorLogical);
	delete m_hmmManager;
	delete m_contextModeling;
	delete m_phoneSet;
	delete m_commandLineManager;
	
	return 0;
}
