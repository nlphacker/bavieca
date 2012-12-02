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


#include "BatchFile.h"
#include "BestPath.h"
#include "CommandLineManager.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "FillerManager.h"
#include "Global.h"
#include "HMMManager.h"
#include "HypothesisLattice.h"
#include "LMManager.h"
#include "Mappings.h"
#include "TimeUtils.h"
#include "TrnFile.h"
#include "Viterbi.h"
#include "ViterbiX.h"

using namespace Bavieca;

// main for the tool "latticeeditor"
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager *m_commandLineManager = new CommandLineManager("latticeeditor",
			SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		m_commandLineManager->defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
		m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,true);
		m_commandLineManager->defineParameter("-lm","language model",PARAMETER_TYPE_FILE,true);
		m_commandLineManager->defineParameter("-ngram","language model",PARAMETER_TYPE_STRING,
			true,"zerogram|unigram|bigram|trigram");
		m_commandLineManager->defineParameter("-bat","lattices to process",PARAMETER_TYPE_FILE,false);
		m_commandLineManager->defineParameter("-act","action to perform",PARAMETER_TYPE_FOLDER,
			false,"wer|pp|align|compact|rescore|lm|addpath");
		//m_commandLineManager->defineParameter("-for","output format",PARAMETER_TYPE_STRING,true,"binary|text","binary");
		m_commandLineManager->defineParameter("-trn","transcription file",PARAMETER_TYPE_FILE,true);
		m_commandLineManager->defineParameter("-hyp","hypotheses file",PARAMETER_TYPE_FILE,true);
		m_commandLineManager->defineParameter("-hypf","hypotheses format",PARAMETER_TYPE_FILE,true,"trn|ctm","trn");
		m_commandLineManager->defineParameter("-ip","insertion penalty",PARAMETER_TYPE_FLOAT,true);
		m_commandLineManager->defineParameter("-ipf","file containing insertion penalties for filler symbols",PARAMETER_TYPE_FILE,true);
		m_commandLineManager->defineParameter("-ams","acoustic model scale factor",PARAMETER_TYPE_FLOAT,true);
		m_commandLineManager->defineParameter("-lms","language model scale factor",PARAMETER_TYPE_FLOAT,true);
		m_commandLineManager->defineParameter("-res","rescoring method",PARAMETER_TYPE_STRING,true,"likelihood|pp","likelihood");
		m_commandLineManager->defineParameter("-conf","confidence annotation method",PARAMETER_TYPE_STRING,true,"posteriors|accumulated|maximum","maximum");
		m_commandLineManager->defineParameter("-map","file containing word mappings for WER computation",PARAMETER_TYPE_FILE,true);	
		m_commandLineManager->defineParameter("-vrb","verbose output",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");			
	
		// parse the parameters
		if (m_commandLineManager->parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the parameter values
		const char *m_strAction = m_commandLineManager->getParameterValue("-act");
		const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
		const char *m_strFileLexicon = m_commandLineManager->getParameterValue("-lex");
		const char *m_strFileLattice = NULL;
		const char *m_strFileBatch = NULL;
		bool m_bVerbose = m_commandLineManager->getBoolParameterValue("-vrb");
		
		// load the phone set
		PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
		m_phoneSet->load();
		
		// load the lexicon
		LexiconManager *m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
		m_lexiconManager->load();
		if (m_commandLineManager->isParameterSet("-ip")) {
		
			// global insertion penalty
			float fInsertionPenalty = atof(m_commandLineManager->getParameterValue("-ip")); 
			m_lexiconManager->attachLexUnitPenalties(fInsertionPenalty,fInsertionPenalty);
			
			// insertion penalty for fillers
			if (m_commandLineManager->isParameterSet("-ipf")) {
			
				const char *strFileIPFiller = m_commandLineManager->getParameterValue("-ipf");
				FillerManager m_fillerManager(strFileIPFiller);	
				m_fillerManager.load();
				m_fillerManager.attachInsertionPenaltyFillers(m_lexiconManager);
			}
		}
		
		if (m_commandLineManager->isParameterSet("-lat")) {
			m_strFileLattice = m_commandLineManager->getParameterValue("-lat");
		} else {
			assert(m_commandLineManager->isParameterSet("-bat") == true);
			m_strFileBatch = m_commandLineManager->getParameterValue("-bat");
		}
		
		LMManager *m_lmManager = NULL;
		if (m_commandLineManager->isParameterSet("-lm")) {
		
			const char *strFileLanguageModel = m_commandLineManager->getParameterValue("-lm");
			const char *strLanguageModelFormat = "ARPA";
			const char *strLanguageModelType = "ngram";
			const char *strLanguageModelNGram = m_commandLineManager->getParameterValue("-ngram");
			
			// load the language model
			m_lmManager = new LMManager(m_lexiconManager,
													strFileLanguageModel,
													strLanguageModelFormat,
													strLanguageModelType,
													strLanguageModelNGram); 
			m_lmManager->load();
			
			// build the lm-graph
			m_lmManager->buildLMGraph();
		}
		
		// lattice WER computation
		if (strcmp(m_strAction,"wer") == 0) {
			
			const char *strFileTrn = m_commandLineManager->getParameterValue("-trn");
			const char *strFileMappings = m_commandLineManager->getParameterValue("-map");
			
			// load the mappings if any
			Mappings *mappings = NULL;
			if (strFileMappings) {
				mappings = new Mappings(strFileMappings);
				if (mappings->load() == false) {
					BVC_ERROR << "unable to load the mappings";
				}
			}
			
			// batch mode
			if (m_strFileBatch != NULL) {
			
				// load the batch file
				BatchFile *m_batchFile = new BatchFile(m_strFileBatch,"lattice|utteranceID");
				m_batchFile->load();
				
				// load the transcription file
				TrnFile *trnFile = new TrnFile(strFileTrn);
				trnFile->load();
				
				LatticeWER latticeWERAll;
				LatticeDepth latticeDepthAll;
				HypothesisLattice::reset(&latticeWERAll);
				HypothesisLattice::reset(&latticeDepthAll);
				
				for(unsigned int i=0 ; i<m_batchFile->size() ; ++i) {
					
					const char *strFileLatticeInput = m_batchFile->getField(i,"lattice");
					const char *strUtteranceId = m_batchFile->getField(i,"utteranceID");
					
					printf("id: %s\n",strUtteranceId);
					
					// load the lattice
					HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
					hypothesisLattice->load(strFileLatticeInput);	
					hypothesisLattice->check();
					
					hypothesisLattice->forwardEdgeMerge();
					hypothesisLattice->backwardEdgeMerge();
						
					// get the transcription
					const char *strTranscription = trnFile->getTranscription(strUtteranceId);
					if (strTranscription == NULL) {
						BVC_ERROR << "no transcription for utterance: \"" << strUtteranceId << "\" was found";
					}
					
					// extract lexical units from the transcription
					VLexUnit vLexUnitsTranscription;
					bool bAllKnown;
					m_lexiconManager->getLexUnits(strTranscription,vLexUnitsTranscription,bAllKnown);
					
					// compute the lattice WER
					LatticeWER *latticeWER = hypothesisLattice->computeWER(vLexUnitsTranscription,mappings);
					if (latticeWER == NULL) {
						BVC_ERROR << "unable to compute the WER for the utterance \"" << strUtteranceId << "\"";
					}
					HypothesisLattice::add(&latticeWERAll,latticeWER);
					
					// compute the lattice depth
					LatticeDepth *latticeDepth = hypothesisLattice->computeDepth();
					HypothesisLattice::add(&latticeDepthAll,latticeDepth);
						
					delete latticeWER;
					delete latticeDepth;
					delete hypothesisLattice;
				}
				HypothesisLattice::print(&latticeWERAll);
				HypothesisLattice::print(&latticeDepthAll);
				delete trnFile;	
				delete m_batchFile;	
			}
		}
		// lattice alignment (and HMM-state marking)
		else if (strcmp(m_strAction,"align") == 0) {
		
			const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");
			float m_fBeamWidth = 2000;
			
			// load the acoustic models
			HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
			m_hmmManager->load(m_strFileModels);
			m_hmmManager->initializeDecoding();
			
			// create the aligner object
			Viterbi *m_viterbi = new Viterbi(m_phoneSet,m_hmmManager,m_lexiconManager,m_fBeamWidth);
			
			// load the batch file
			BatchFile *m_batchFile = new BatchFile(m_strFileBatch,"latticeIn|features|latticeOut");
			m_batchFile->load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {	
				
				const char *strFileLatticeInput = m_batchFile->getField(i,"latticeIn");
				const char *strFileFeatures = m_batchFile->getField(i,"features");
				const char *strFileLatticeOutput = m_batchFile->getField(i,"latticeOut");
				
				// load the lattice
				HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
				hypothesisLattice->load(strFileLatticeInput);
				
				hypothesisLattice->check();
				
				// load the features
				FeatureFile *featureFile = new FeatureFile(strFileFeatures,MODE_READ);
				featureFile->load();
				int iFeatures = -1;
				float *fFeatures = featureFile->getFeatureVectors(&iFeatures);
				if (iFeatures != hypothesisLattice->getFrames()) {
					BVC_ERROR << "features and lattice do not match";
				}
				
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
					
				// mark the lattice with acoustic scores, phone-boundaries and HMM-states
				hypothesisLattice->hmmMarking(m_hmmManager);
				
				// mark the lattice with phone-alignments
				if (m_viterbi->align(fFeatures,iFeatures,hypothesisLattice) == false) {
					BVC_ERROR << "unable to generate phone-level alignments for the lattice";
				}
				
				// write the lattice to disk
				hypothesisLattice->store(strFileLatticeOutput);
				
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice->getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice->getFrames());
					
				delete [] fFeatures;
				delete featureFile;
				delete hypothesisLattice;
				
				//break;
			}
			delete m_batchFile;	
			
			delete m_viterbi;
			delete m_hmmManager;	
		}
		// attach language model log-likelihoods and insertion penalties
		else if (strcmp(m_strAction,"lm") == 0) {
		
			// load the batch file
			BatchFile *m_batchFile = new BatchFile(m_strFileBatch,"latticeIn|latticeOut");
			m_batchFile->load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {	
				
				const char *strFileLatticeInput = m_batchFile->getField(i,"latticeIn");
				const char *strFileLatticeOutput = m_batchFile->getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
				hypothesisLattice->load(strFileLatticeInput);
				
				// attach language model scores to the edges in the lattice
				hypothesisLattice->attachLMProbabilities(m_lmManager);
				
				// attach insertion penalties
				if (m_commandLineManager->isParameterSet("-ip")) {
					hypothesisLattice->attachInsertionPenalty(m_lexiconManager);
				}
				
				// write the lattice to disk
				hypothesisLattice->store(strFileLatticeOutput);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice->getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice->getFrames());
				
				delete hypothesisLattice;
			}
		
			delete m_batchFile;
			delete m_lmManager;
		}
		// add a path to the lattice in case it is not already there
		else if (strcmp(m_strAction,"addpath") == 0) {
			
			// load the batch file
			BatchFile *m_batchFile = new BatchFile(m_strFileBatch,"latticeIn|alignment|latticeOut");
			m_batchFile->load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {	
				
				const char *strFileLatticeInput = m_batchFile->getField(i,"latticeIn");
				const char *strFileAlignment = m_batchFile->getField(i,"alignment");
				const char *strFileLatticeOutput = m_batchFile->getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
				hypothesisLattice->load(strFileLatticeInput);
				
				// load the alignment
				Alignment *alignment = Alignment::load(strFileAlignment,m_lexiconManager);
				assert(alignment);
				
				// get the word sequence from the alignment
				VLexUnit vLexUnits;
				VWordAlignment *vWordAlignment = alignment->getWordAlignment();
				for(VWordAlignment::iterator it = vWordAlignment->begin() ; it != vWordAlignment->end() ; ++it) {
					LexUnit *lexUnit = m_lexiconManager->getLexUnitPron((*it)->iLexUnitPron);
					if (m_lexiconManager->isStandard(lexUnit)) {
						vLexUnits.push_back(m_lexiconManager->getLexUnitPron((*it)->iLexUnitPron));
					}
				}
				
				// get the lattice WER taking the path to be added as the reference
				LatticeWER *latticeWER = hypothesisLattice->computeWER(vLexUnits,NULL);
				if (latticeWER == NULL) {
					BVC_ERROR << "unable to compute the WER";
				}
				// if the path is not in the lattice then we need to insert it
				if (latticeWER->iErrors != 0) {
					hypothesisLattice->addPath(alignment,true);
				}
				delete latticeWER;	
				
				// write the lattice to disk
				hypothesisLattice->store(strFileLatticeOutput);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice->getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice->getFrames());
				
				delete alignment;
				delete hypothesisLattice;
			}
		
			delete m_batchFile;
			delete m_lmManager;
		}
		// confidence/posterior probabilities
		else if (strcmp(m_strAction,"pp") == 0) {
		
			assert(m_commandLineManager->isParameterSet("-ams"));
			assert(m_commandLineManager->isParameterSet("-lms"));
		
			// get scale factors
			float fScaleAM = atof(m_commandLineManager->getParameterValue("-ams"));
			float fScaleLM = atof(m_commandLineManager->getParameterValue("-lms"));
			
			// load the batch file
			BatchFile *m_batchFile = new BatchFile(m_strFileBatch,"latticeIn|latticeOut");
			m_batchFile->load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {	
				
				const char *strFileLatticeInput = m_batchFile->getField(i,"latticeIn");
				const char *strFileLatticeOutput = m_batchFile->getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
				hypothesisLattice->load(strFileLatticeInput);
				
				// attach insertion penalties
				if (m_commandLineManager->isParameterSet("-ip")) {
					hypothesisLattice->attachInsertionPenalty(m_lexiconManager);
				}
				
				// compute forward-backward scores
				hypothesisLattice->computeForwardBackwardScores(fScaleAM,fScaleLM);
				hypothesisLattice->computePosteriorProbabilities();
				
				// compute confidence estimates?
				if (m_commandLineManager->isParameterSet("-conf")) {
					hypothesisLattice->computeConfidenceScore(CONFIDENCE_MEASURE_MAXIMUM);	
				}
					
				// write the lattice to disk
				hypothesisLattice->store(strFileLatticeOutput);
				
				// text format
				char strFileLatticeOutputTxt[1024];
				sprintf(strFileLatticeOutputTxt,"%s.txt",strFileLatticeOutput);
				hypothesisLattice->store(strFileLatticeOutputTxt,FILE_FORMAT_TEXT);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice->getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice->getFrames());
				
				delete hypothesisLattice;
			}
		
			delete m_batchFile;
			delete m_lmManager;
		}
		// compacting
		else if (strcmp(m_strAction,"compact") == 0) {
		
			// load the batch file
			BatchFile *m_batchFile = new BatchFile(m_strFileBatch,"latticeIn|latticeOut");
			m_batchFile->load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {	
				
				const char *strFileLatticeInput = m_batchFile->getField(i,"latticeIn");
				const char *strFileLatticeOutput = m_batchFile->getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
				hypothesisLattice->load(strFileLatticeInput);
				
				// forward/backward compacting
				hypothesisLattice->forwardEdgeMerge();
				hypothesisLattice->backwardEdgeMerge();	
				
				// write the lattice to disk
				hypothesisLattice->store(strFileLatticeOutput);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice->getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice->getFrames());
				
				delete hypothesisLattice;
			}
		
			delete m_batchFile;
			delete m_lmManager;
		}
		// rescoring
		else if (strcmp(m_strAction,"rescore") == 0) {
		
			// get the rescoring method
			const char *strRescoringMethod = m_commandLineManager->getParameterValue("-res");
			
			float fScaleAM = 0.0;
			float fScaleLM = 0.0;
			if (strcmp(strRescoringMethod,RESCORING_METHOD_LIKELIHOOD) == 0) {
				
				assert(m_commandLineManager->isParameterSet("-ams"));
				assert(m_commandLineManager->isParameterSet("-lms"));	
			
				// get scale factors
				fScaleAM = atof(m_commandLineManager->getParameterValue("-ams"));
				fScaleLM = atof(m_commandLineManager->getParameterValue("-lms"));
			}
				
			// get the hypothesis format
			const char *strFileHypFormat = m_commandLineManager->getParameterValue("-hypf");
			const char *strFileHypothesis = NULL;
			bool bTrn = true;
			const char *strBatchType = "lattice|utteranceId";
			FileOutput *fileHyp = NULL;
			if ((!strFileHypFormat) || strcmp(strFileHypFormat,"trn") == 0) {	
				strFileHypothesis = m_commandLineManager->getParameterValue("-hyp");
				fileHyp = new FileOutput(strFileHypothesis,false);
				fileHyp->open();	
			} else {
				assert(strcmp(strFileHypFormat,"ctm") == 0);
				strBatchType = "lattice|utteranceId|hypothesis";
				bTrn = false;
			}
				
			// load the batch file
			BatchFile *m_batchFile = new BatchFile(m_strFileBatch,strBatchType);
			m_batchFile->load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {	
				
				const char *strFileLatticeInput = m_batchFile->getField(i,"lattice");
				const char *strUtteranceId = m_batchFile->getField(i,"utteranceId");
				if (bTrn == false) {
					strFileHypothesis = m_batchFile->getField(i,"hypothesis");
				}
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice *hypothesisLattice = new HypothesisLattice(m_phoneSet,m_lexiconManager,m_bVerbose);
				hypothesisLattice->load(strFileLatticeInput);
				
				// attach insertion penalties
				if (m_commandLineManager->isParameterSet("-ip")) {
					hypothesisLattice->attachInsertionPenalty(m_lexiconManager);
				}	
				
				// likelihood based rescoring: set scaling factors
				if (strcmp(strRescoringMethod,RESCORING_METHOD_LIKELIHOOD) == 0) {	
					hypothesisLattice->setScalingFactors(fScaleAM,fScaleLM);
				} 
				
				// lattice rescoring
				BestPath *bestPath = hypothesisLattice->rescore(strRescoringMethod);
				if (bestPath != NULL) {
					if (bTrn) {
						bestPath->write(fileHyp->getStream(),strUtteranceId);
					} else {	
						FileOutput fileHyp(strFileHypothesis,false);
						fileHyp.open();	
						bestPath->write(fileHyp.getStream(),strUtteranceId,strUtteranceId,0.0,false,true,true);
						fileHyp.close();
					}
					delete bestPath;
				}
				
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice->getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice->getFrames());
				
				delete hypothesisLattice;
			}
			if (bTrn) {
				fileHyp->close();
				delete fileHyp;
			}	
		
			delete m_batchFile;
		}
		// unsupported action
		else {	
			printf("action not supported\n");
			return -1;
		}
		
		// clean-up
		delete m_lexiconManager;
		delete m_phoneSet;
		delete m_commandLineManager;
	
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}

	return 0;
}
