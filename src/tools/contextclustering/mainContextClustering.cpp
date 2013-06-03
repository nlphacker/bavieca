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

#include <stdexcept>

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

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("contextclustering",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-rul","phonetic rules",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-ww","within-word context modeling order",PARAMETER_TYPE_STRING,true,NULL,HMM_CONTEXT_MODELING_TRIPHONES_STR);
		commandLineManager.defineParameter("-cw","cross-word context modeling order",PARAMETER_TYPE_STRING,true,NULL,HMM_CONTEXT_MODELING_TRIPHONES_STR);
		commandLineManager.defineParameter("-met","clustering method",PARAMETER_TYPE_STRING,true,"local|global",CLUSTERING_METHOD_LOCAL_STR);
		commandLineManager.defineParameter("-mrg","whether to perform bottom up merging",PARAMETER_TYPE_BOOLEAN,true,NULL,"yes");
		commandLineManager.defineParameter("-acc","logical accumulator list",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-occ","minimum cluster occupation",PARAMETER_TYPE_FLOAT,true,NULL,"200");
		commandLineManager.defineParameter("-gan","minimum likelihood gain",PARAMETER_TYPE_FLOAT,true,NULL,"2000");
		commandLineManager.defineParameter("-out","output acoustic models",PARAMETER_TYPE_FILE,false);
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
	
		// get the parameter values	
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");	
		const char *strFilePhoneticRules = commandLineManager.getParameterValue("-rul");
		unsigned char iContextModelingOrderWW = Accumulator::getContextModelingOrder(commandLineManager.getParameterValue("-ww"));
		if (iContextModelingOrderWW == UCHAR_MAX) {
			BVC_ERROR << "wrong within-word context";
		}
		unsigned char iContextModelingOrderCW = Accumulator::getContextModelingOrder(commandLineManager.getParameterValue("-cw"));
		if (iContextModelingOrderCW == UCHAR_MAX) {
			BVC_ERROR << "wrong cross-word context";
		}
		bool bGlobalClustering = false;
		if (strcmp(commandLineManager.getParameterValue("-met"),"global") == 0) {
			bGlobalClustering = true;
		}	
		bool bBottomUpMerging = true;
		if (strcmp(commandLineManager.getParameterValue("-mrg"),"no") == 0) {
			bBottomUpMerging = false;
		}	
		const char *strFileAccList = commandLineManager.getParameterValue("-acc");
		float fMinimumClusterOccupation = atof(commandLineManager.getParameterValue("-occ"));
		float fMinimumLikelihoodGain = atof(commandLineManager.getParameterValue("-gan"));
		const char *strFileModelsOutput = commandLineManager.getParameterValue("-out");
		
		// load the phonetic symbol set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the HMMs
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_ESTIMATION);	
		hmmManager.load(strFileModels);
		
		double dTimeStart = TimeUtils::getTimeMilliseconds();
		
		// load the accumulators
		AccMetadata metadata;
		MAccumulatorLogical mAccumulatorLogical;
		Accumulator::loadAccumulatorList(strFileAccList,mAccumulatorLogical,metadata);
	
		double dTimeEndLoading = TimeUtils::getTimeMilliseconds();
		
		//printf("ww: %d cw: %d\n",iContextModelingOrderWW,iContextModelingOrderCW);
		
		// do the clustering
		ContextModeling contextModeling(hmmManager.getFeatureDim(),
			COVARIANCE_MODELLING_TYPE_DIAGONAL,
			&phoneSet,iContextModelingOrderWW,iContextModelingOrderCW,strFilePhoneticRules,bGlobalClustering,
			bBottomUpMerging,fMinimumClusterOccupation,fMinimumLikelihoodGain);
	
		if (contextModeling.clusterContextDependentUnits(mAccumulatorLogical,&hmmManager) == false) {
			BVC_ERROR << "problem found during context dependent clustering";
		}
		
		double dTimeEndClustering = TimeUtils::getTimeMilliseconds();
		
		printf("loading time:    %.2f seconds\n",(dTimeEndLoading-dTimeStart)/1000.0);
		printf("clustering time: %.2f seconds\n",(dTimeEndClustering-dTimeEndLoading)/1000.0);
		
		// create the output HMMs
		hmmManager.store(strFileModelsOutput);
		
		// clean-up
		Accumulator::destroy(mAccumulatorLogical);
	
	} catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;
}
