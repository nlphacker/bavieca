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
#include "FillerManager.h"
#include "HMMManager.h"
#include "LMManager.h"
#include "TimeUtils.h"
#include "WFSABuilder.h"

using namespace Bavieca;

// main for the tool "wfsbuilder"
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("wfsabuilder",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lm","language model",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-ngram","n-gram",PARAMETER_TYPE_STRING,false,"zerogram|unigram|bigram|trigram");
		commandLineManager.defineParameter("-scl","language model scaling factor",PARAMETER_TYPE_FLOAT,false);
		commandLineManager.defineParameter("-ip","insertion penalty (standard lexical units)",PARAMETER_TYPE_FLOAT,false);
		commandLineManager.defineParameter("-ips","insertion penalty (silence and filler lexical units)",PARAMETER_TYPE_FLOAT,false);	
		commandLineManager.defineParameter("-ipf","filler specific insertion penalties",PARAMETER_TYPE_FILE,true);	
		commandLineManager.defineParameter("-srg","semiring used to do weight pushing",PARAMETER_TYPE_STRING,true,"none|tropical|log","log");
		commandLineManager.defineParameter("-net","decoding network to build",PARAMETER_TYPE_FILE,false);
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the parameters
		const char *strFilePhoneticSet = commandLineManager.getStrParameterValue("-pho");
		const char *strFileLexicon = commandLineManager.getStrParameterValue("-lex");
		const char *strFileLanguageModel = commandLineManager.getStrParameterValue("-lm");
		const char *strLanguageModelFormat = "ARPA";
		const char *strLanguageModelType = "ngram";
		const char *strLanguageModelNGram = commandLineManager.getStrParameterValue("-ngram");
		const char *strFileModels = commandLineManager.getStrParameterValue("-mod");
		float fLMScalingFactor = commandLineManager.getFloatParameterValue("-scl");
		float fInsertionPenaltyStandard = commandLineManager.getFloatParameterValue("-ip");
		float fInsertionPenaltyFiller = commandLineManager.getFloatParameterValue("-ips");
		const char *strFileInsertionPenaltyFiller = NULL;
		if (commandLineManager.isParameterSet("-ipf")) {
			strFileInsertionPenaltyFiller = commandLineManager.getStrParameterValue("-ipf");
		}
		const char *strFileDecodingNetwork = commandLineManager.getStrParameterValue("-net");
	
		// parameter injection
		
		// wsj
		/*string strPhoneticSetFile = "/home/hmmdecoder/src/config/wsj/phoneset.txt";
		string strLexiconFile = "/home/hmmdecoder/src/config/WFST/wsj-5k.lex";
		string strLanguageModelFile = "/home/hmmdecoder/src/config/WFST/wsj-5k-cnp.arpa";
		//string strLanguageModelFile = "/home/hmmdecoder/src/config/wsj/wsj-5k-cnp-bigram.arpa";	
		//string strFileModels = "/home/hmmtrainer/src/data/models80h/july26/gi_500_2_100_sil/models28.bin";
		//string strFileModels = "/home/hmmtrainer/src/data/models3|30min/(6)/models/models05.bin";
		//string strFileModels = "/home/hmmtrainer/src/data/models3|30min/(6)/models/models16.bin";
		string strFileModels = "/home/hmmdecoder/src/config/WFST/wsj/models25.bin";
		string strFileList = "/home/hmmdecoder/src/config/wsj/list.txt";*/
		
		/*string strPhoneticSetFile = "/home/speech/wsj/scripts/test/config/phoneset.txt";
		string strLexiconFile = "/home/speech/wsj/scripts/test/lm/wsj-5k.lex";
		string strLanguageModelFile = "/home/speech/wsj/scripts/test/lm/wsj-5k-cnp.arpa";
		string strFileModels = "/home/speech/wsj/models/may18th/1400_400_0.05_paramMyFixedHamming/gi/models25.bin";
		string strFileNetwork = "/home/speech/wsj/scripts/test/wfsa/unigram.bin";
		//float fLMScalingFactor = 0.0;
		float fLMScalingFactor = 25.0;
		//float fInsertionPenalty = 0.0;
		float fInsertionPenalty = -18.0;	*/
		
		// ScienceTutor
		/*string strPhoneticSetFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/phoneset.txt";
		//string strLexiconFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/lexicon.lex";
		string strLexiconFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/lexicon2Filler.lex";
		//string strLexiconFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/lexicon1Filler.lex";
		//string strLexiconFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/lexiconNoFiller.lex";
		string strLanguageModelFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/lm.arpa";
		//string strLexiconFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/MS/lexicon.lex";
		//string strLanguageModelFile = "/home/hmmdecoder/src/config/WFST/ScienceTutor/MS/lm.arpa";
		//string strFileModels = "/home/hmmdecoder/src/config/WFST/ScienceTutor/models05.bin";
		string strFileModels = "/home/hmmdecoder/src/config/WFST/ScienceTutor/models20.bin";
		float fLMScalingFactor = 25.0;
		float fInsertionPenalty = -20.0;*/
		
		// kids
		/*string strPhoneticSetFile = "/home/hmmdecoder/src/config/WFST/kids/phoneset.txt";
		//string strLexiconFile = "/home/hmmdecoder/src/config/WFST/kids/story.lex";
		string strLexiconFile = "/home/hmmdecoder/src/config/WFST/kids/storyNoFiller.lex";
		string strLanguageModelFile = "/home/hmmdecoder/src/config/WFST/kids/racer-lm.arpa.gen";
		//string strFileModels = "/home/hmmdecoder/src/config/WFST/kids/models05.bin";
		string strFileModels = "/home/hmmdecoder/src/config/WFST/kids/models25.bin";
		float fLMScalingFactor = 40.0;
		float fInsertionPenalty = 0.0;*/
		
		/*string strLanguageModelFileFormat = "text";
		string strLanguageModelFormat = "CMU";
		string strLanguageModelType = "ngram";
		string strLanguageModelNGram = "unigram";*/
		//string strLanguageModelNGram = "bigram";
		//string strLanguageModelNGram = "trigram";
	
		// load the phone set
		PhoneSet phoneSet(strFilePhoneticSet);
		phoneSet.load();
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		hmmManager.initializeDecoding();	
		
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet);
		lexiconManager.load();
		// set default insertion penalty to each lexical unit in the lexicon
		lexiconManager.attachLexUnitPenalties(fInsertionPenaltyStandard,fInsertionPenaltyFiller);
		// set specific insertion penalties if available
		if (strFileInsertionPenaltyFiller != NULL) {
			FillerManager fillerManager(strFileInsertionPenaltyFiller);	
			fillerManager.load();
			fillerManager.attachInsertionPenaltyFillers(&lexiconManager);
		}
		lexiconManager.print();  
		
		// load the language model
		LMManager lmManager(&lexiconManager,strFileLanguageModel,strLanguageModelFormat,
									strLanguageModelType,strLanguageModelNGram); 
		lmManager.load();
		lmManager.print();
	
		// build the decoding network
		int iNGram = LMManager::getNGram(strLanguageModelNGram);
		WFSABuilder wfsaBuilder(&phoneSet,&hmmManager,&lexiconManager,&lmManager,iNGram,fLMScalingFactor);
		
		double dBegin = TimeUtils::getTimeMilliseconds();
		
		WFSAcceptor *wfsAcceptor = wfsaBuilder.build();
		if (!wfsAcceptor) {
			BVC_ERROR << "unable to create the WFSA";
		}
		wfsAcceptor->print();
		
		// store the acceptor to disk
		wfsAcceptor->store(&lexiconManager,strFileDecodingNetwork);
		
		double dEnd = TimeUtils::getTimeMilliseconds();
		double dMillisecondsInterval = dEnd - dBegin;
		printf("building time: %.2f\n",dMillisecondsInterval/1000.0);
	
		delete wfsAcceptor;
	
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;
}

