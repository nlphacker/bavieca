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
#include <iostream>
#include <cstdlib>

#include "Viterbi.h"
#include "AlignmentFile.h"
#include "AudioFile.h"
#include "BatchFile.h"
#include "BestPath.h"
#include "CommandLineManager.h"
#include "ConfigurationDynamicDecoder.h"
#include "ConfigurationFeatures.h"
#include "DynamicNetworkX.h"
#include "DynamicDecoderX.h"
#include "NetworkBuilderX.h"
#include "FeatureExtractor.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "FillerManager.h"
#include "HMMManager.h"
#include "LexiconManager.h"
#include "LexUnitsFile.h"
#include "LMARPA.h"
#include "LMFSM.h"
#include "LMManager.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the tool "lmfsm"
int main(int argc, char *argv[]) {

	try {

		// define command line parameters
		CommandLineManager commandLineManager("lmfsm",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);	
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-lm","input language model",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-fsm","finite state machine",PARAMETER_TYPE_FILE,false);
		
		// process command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get command line parameters
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileLexicon = commandLineManager.getParameterValue("-lex");
		const char *strFileLM = commandLineManager.getParameterValue("-lm");
		const char *strFileFSM = commandLineManager.getParameterValue("-fsm");
	
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
	
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet); 
		lexiconManager.load();
		
		// load the language model in ARPA format
		LMARPA lmARPA(&lexiconManager,strFileLM);
		lmARPA.load();
		
		// build the FSM from the lm in ARPA format
		LMFSM lmFSM(&lexiconManager);
		lmFSM.build(&lmARPA);
		//cout << "likelihood: " << lmFSM.computeLikelihood("THE MAGNETS STICK TO THE WIRE") << endl;
		//cout << "likelihood: " << lmFSM.computeLikelihood("THE INVESTMENT IS GOOD BUT NOT GREAT YET") << endl;
		lmFSM.store(strFileFSM);
	
		// load the FSM from the lm in ARPA format
		LMFSM lmFSM2(&lexiconManager);
		lmFSM2.load(strFileFSM);
	} 
	catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;
}
