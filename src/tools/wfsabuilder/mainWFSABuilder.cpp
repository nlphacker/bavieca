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

	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("wfsabuilder",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lm","language model",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-ngram","n-gram",PARAMETER_TYPE_STRING,false,"zerogram|unigram|bigram|trigram");
	m_commandLineManager->defineParameter("-scl","language model scaling factor",PARAMETER_TYPE_FLOAT,false);
	m_commandLineManager->defineParameter("-ip","insertion penalty (standard lexical units)",PARAMETER_TYPE_FLOAT,false);
	m_commandLineManager->defineParameter("-ips","insertion penalty (silence and filler lexical units)",PARAMETER_TYPE_FLOAT,false);	
	m_commandLineManager->defineParameter("-ipf","filler specific insertion penalties",PARAMETER_TYPE_FILE,true);	
	m_commandLineManager->defineParameter("-srg","semiring used to do weight pushing",PARAMETER_TYPE_STRING,true,"none|tropical|log","log");
	m_commandLineManager->defineParameter("-net","decoding network to build",PARAMETER_TYPE_FILE,false);
	
	// parse the parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get the parameters
	const char *m_strFilePhoneticSet = m_commandLineManager->getStrParameterValue("-pho");
	const char *m_strFileLexicon = m_commandLineManager->getStrParameterValue("-lex");
	const char *m_strFileLanguageModel = m_commandLineManager->getStrParameterValue("-lm");
	const char *m_strLanguageModelFormat = "ARPA";
	const char *m_strLanguageModelType = "ngram";
	const char *m_strLanguageModelNGram = m_commandLineManager->getStrParameterValue("-ngram");
	const char *m_strFileModels = m_commandLineManager->getStrParameterValue("-mod");
	float m_fLMScalingFactor = m_commandLineManager->getFloatParameterValue("-scl");
	float m_fInsertionPenaltyStandard = m_commandLineManager->getFloatParameterValue("-ip");
	float m_fInsertionPenaltyFiller = m_commandLineManager->getFloatParameterValue("-ips");
	const char *m_strFileInsertionPenaltyFiller = NULL;
	if (m_commandLineManager->isParameterSet("-ipf")) {
		m_strFileInsertionPenaltyFiller = m_commandLineManager->getStrParameterValue("-ipf");
	}
	const char *m_strFileDecodingNetwork = m_commandLineManager->getStrParameterValue("-net");

   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneticSet);
   m_phoneSet->load();
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeDecoding();	
   
   // load the lexicon
   LexiconManager *m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet);
   m_lexiconManager->load();
   // set default insertion penalty to each lexical unit in the lexicon
   m_lexiconManager->attachLexUnitPenalties(m_fInsertionPenaltyStandard,m_fInsertionPenaltyFiller);
   // set specific insertion penalties if available
   if (m_strFileInsertionPenaltyFiller != NULL) {
		FillerManager m_fillerManager(m_strFileInsertionPenaltyFiller);	
		m_fillerManager.load();
		m_fillerManager.attachInsertionPenaltyFillers(m_lexiconManager);
   }
   m_lexiconManager->print();  
   
   // load the language model
	LMManager *m_lmManager = new LMManager(m_lexiconManager,
											m_strFileLanguageModel,
											m_strLanguageModelFormat,
											m_strLanguageModelType,
											m_strLanguageModelNGram); 
	m_lmManager->load();
	m_lmManager->print();

	// build the decoding network
	int m_iNGram = LMManager::getNGram(m_strLanguageModelNGram);
	WFSABuilder *m_wfsaBuilder = new WFSABuilder(m_phoneSet,m_hmmManager,m_lexiconManager,m_lmManager,m_iNGram,m_fLMScalingFactor);
	
	double dBegin = TimeUtils::getTimeMilliseconds();
	
	WFSAcceptor *m_wfsAcceptor = m_wfsaBuilder->build();
	if (m_wfsAcceptor == NULL) {
		return -1;
	}
	m_wfsAcceptor->print();
	
	// store the acceptor to disk
	m_wfsAcceptor->store(m_lexiconManager,m_strFileDecodingNetwork);
	
   double dEnd = TimeUtils::getTimeMilliseconds();
   double dMillisecondsInterval = dEnd - dBegin;
	printf("building time: %.2f\n",dMillisecondsInterval/1000.0);

	delete m_wfsaBuilder;
	delete m_wfsAcceptor;
	delete m_commandLineManager;
	delete m_phoneSet;
	delete m_hmmManager;
	delete m_lexiconManager;
	delete m_lmManager;
	
	return 0;
}

