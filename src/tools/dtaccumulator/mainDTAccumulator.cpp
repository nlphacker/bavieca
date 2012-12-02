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
#include "DTAccumulator.h"
#include "MLFFile.h"
#include "PhoneSet.h"

using namespace Bavieca;
 
// main for the tool "dtaccumulator"
int main(int argc, char *argv[]) {
	
	// define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("mlaccumulator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","file containing the phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","input acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lex","pronunciation lexicon",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-fea","feature folder",PARAMETER_TYPE_FOLDER,false);
	m_commandLineManager->defineParameter("-cfg","feature configuration file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mlf","master label file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lat","directory containing hmm-marked lattices",PARAMETER_TYPE_FOLDER,false);
	m_commandLineManager->defineParameter("-ams","acoustic model scale factor",PARAMETER_TYPE_FLOAT,false);
	m_commandLineManager->defineParameter("-lms","language model scale factor",PARAMETER_TYPE_FLOAT,true,NULL,"1.0");	
	m_commandLineManager->defineParameter("-opt","file containing optional symbols that can be inserted at word boundaries",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-pro","allow multiple pronunciations",PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
	m_commandLineManager->defineParameter("-fwd","forward pruning",PARAMETER_TYPE_FLOAT,true,"[-100.0|-10.0]","-20");
	m_commandLineManager->defineParameter("-bwd","backward pruning",PARAMETER_TYPE_FLOAT,true,"[100.0|10000.0]","800");
	m_commandLineManager->defineParameter("-tre","maximum trellis size (MB)",PARAMETER_TYPE_INTEGER,true,NULL,"500");
	m_commandLineManager->defineParameter("-dAccNum","file to dump accumulators (numerator)",
		PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-dAccDen","file to dump accumulators (denominator)",
		PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-obj","objective function",PARAMETER_TYPE_STRING,true,"MMI|bMMI","MMI");	
	m_commandLineManager->defineParameter("-bst","boosting factor for bMMI",PARAMETER_TYPE_FLOAT,true,NULL,"0.5");
	m_commandLineManager->defineParameter("-can","statistics cancelation",PARAMETER_TYPE_BOOLEAN,true,NULL,"yes");	
	
	// parse the command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// load the parameters
	const char *strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *strFileModels = m_commandLineManager->getParameterValue("-mod");
	const char *strFileLexicon = m_commandLineManager->getParameterValue("-lex");
	const char *strFolderFeatures = m_commandLineManager->getParameterValue("-fea"); 
	const char *strFileFeatureConfiguration = m_commandLineManager->getParameterValue("-cfg"); 
	const char *strFileMLF = m_commandLineManager->getParameterValue("-mlf"); 
	const char *strFolderLattices = m_commandLineManager->getParameterValue("-lat");
	float fScaleAM = atof(m_commandLineManager->getParameterValue("-ams"));
	float fScaleLM = atof(m_commandLineManager->getParameterValue("-lms"));	 
	const char *strFileOptionalSymbols = m_commandLineManager->getParameterValue("-opt");
	bool bMultiplePronunciations = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-pro"));
	float fForwardPruningBeam = atof(m_commandLineManager->getParameterValue("-fwd")); 
	float fBackwardPruningBeam = atof(m_commandLineManager->getParameterValue("-bwd")); 
	int iTrellisMaxSize = atoi(m_commandLineManager->getParameterValue("-tre")); 
	bool bTrellisCache = true;
	int iTrellisCacheMaxSize = atoi(m_commandLineManager->getParameterValue("-tre"));	
	const char *strFileAccumulatorsNum = m_commandLineManager->getParameterValue("-dAccNum");
	const char *strFileAccumulatorsDen = m_commandLineManager->getParameterValue("-dAccDen");
	const char *strObjectiveFunction = m_commandLineManager->getParameterValue("-obj");
	float fBoostingFactor = atof(m_commandLineManager->getParameterValue("-bst"));
	bool bStatisticsCancelation = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-can"));
		
	// create the accumulator object
	DTAccumulator dtAccumulator(strFilePhoneSet,strFileFeatureConfiguration,
		strFolderFeatures,strFileModels,strFileOptionalSymbols,bMultiplePronunciations,strFileLexicon,strFileMLF,
		strFolderLattices,fScaleAM,fScaleLM,strFileAccumulatorsNum,strFileAccumulatorsDen,strObjectiveFunction,
		fBoostingFactor,bStatisticsCancelation,fForwardPruningBeam,fBackwardPruningBeam,iTrellisMaxSize,bTrellisCache,
		iTrellisCacheMaxSize);
			
	dtAccumulator.initialize();
	dtAccumulator.accumulate();
		
	// clean-up	
	delete m_commandLineManager;
	
	return 0;
}

