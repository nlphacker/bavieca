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
		CommandLineManager commandLineManager("latticeeditor",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-lm","language model",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-ngram","language model",PARAMETER_TYPE_STRING,
			true,"zerogram|unigram|bigram|trigram");
		commandLineManager.defineParameter("-bat","lattices to process",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-act","action to perform",PARAMETER_TYPE_FOLDER,
			false,"wer|pp|align|compact|rescore|lm|addpath");
		//commandLineManager.defineParameter("-for","output format",PARAMETER_TYPE_STRING,true,"binary|text","binary");
		commandLineManager.defineParameter("-trn","transcription file",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-hyp","hypotheses file",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-hypf","hypotheses format",PARAMETER_TYPE_FILE,true,"trn|ctm","trn");
		commandLineManager.defineParameter("-ip","insertion penalty",PARAMETER_TYPE_FLOAT,true);
		commandLineManager.defineParameter("-ipf","file containing insertion penalties for filler symbols",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-ams","acoustic model scale factor",PARAMETER_TYPE_FLOAT,true);
		commandLineManager.defineParameter("-lms","language model scale factor",PARAMETER_TYPE_FLOAT,true);
		commandLineManager.defineParameter("-res","rescoring method",PARAMETER_TYPE_STRING,true,"likelihood|pp","likelihood");
		commandLineManager.defineParameter("-conf","confidence annotation method",PARAMETER_TYPE_STRING,true,"posteriors|accumulated|maximum","maximum");
		commandLineManager.defineParameter("-map","file containing word mappings for WER computation",PARAMETER_TYPE_FILE,true);	
		commandLineManager.defineParameter("-vrb","verbose output",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");			
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the parameter values
		const char *strAction = commandLineManager.getParameterValue("-act");
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileLexicon = commandLineManager.getParameterValue("-lex");
		const char *strFileLattice = NULL;
		const char *strFileBatch = NULL;
		bool bVerbose = commandLineManager.getBoolParameterValue("-vrb");
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet); 
		lexiconManager.load();
		if (commandLineManager.isParameterSet("-ip")) {
		
			// global insertion penalty
			float fInsertionPenalty = atof(commandLineManager.getParameterValue("-ip")); 
			lexiconManager.attachLexUnitPenalties(fInsertionPenalty,fInsertionPenalty);
			
			// insertion penalty for fillers
			if (commandLineManager.isParameterSet("-ipf")) {
			
				const char *strFileIPFiller = commandLineManager.getParameterValue("-ipf");
				FillerManager fillerManager(strFileIPFiller);	
				fillerManager.load();
				fillerManager.attachInsertionPenaltyFillers(&lexiconManager);
			}
		}
		
		if (commandLineManager.isParameterSet("-lat")) {
			strFileLattice = commandLineManager.getParameterValue("-lat");
		} else {
			assert(commandLineManager.isParameterSet("-bat"));
			strFileBatch = commandLineManager.getParameterValue("-bat");
		}
		
		LMManager *lmManager = NULL;
		if (commandLineManager.isParameterSet("-lm")) {
		
			const char *strFileLanguageModel = commandLineManager.getParameterValue("-lm");
			const char *strLanguageModelFormat = "ARPA";
			const char *strLanguageModelType = "ngram";
			const char *strLanguageModelNGram = commandLineManager.getParameterValue("-ngram");
			
			// load the language model
			lmManager = new LMManager(&lexiconManager,strFileLanguageModel,strLanguageModelFormat,
													strLanguageModelType,strLanguageModelNGram); 
			lmManager->load();
			lmManager->buildLMGraph();
		}
		
		// lattice WER computation
		if (strcmp(strAction,"wer") == 0) {
			
			const char *strFileTrn = commandLineManager.getParameterValue("-trn");
			const char *strFileMappings = commandLineManager.getParameterValue("-map");
			
			// load the mappings if any
			Mappings *mappings = NULL;
			if (strFileMappings) {
				mappings = new Mappings(strFileMappings);
				mappings->load();
			}
			
			// batch mode
			if (strFileBatch) {
			
				// load the batch file
				BatchFile batchFile(strFileBatch,"lattice|utteranceID");
				batchFile.load();
				
				// load the transcription file
				TrnFile trnFile(strFileTrn);
				trnFile.load();
				
				LatticeWER latticeWERAll;
				LatticeDepth latticeDepthAll;
				HypothesisLattice::reset(&latticeWERAll);
				HypothesisLattice::reset(&latticeDepthAll);
				
				for(unsigned int i=0 ; i<batchFile.size() ; ++i) {
					
					const char *strFileLatticeInput = batchFile.getField(i,"lattice");
					const char *strUtteranceId = batchFile.getField(i,"utteranceID");
					
					printf("id: %s\n",strUtteranceId);
					
					// load the lattice
					HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
					hypothesisLattice.load(strFileLatticeInput);	
					hypothesisLattice.check();
					
					hypothesisLattice.forwardEdgeMerge();
					/*char strFile[2000];
					sprintf(strFile,"%s_out.txt",strFileLatticeInput);
					hypothesisLattice.storeTextFormat(strFile);*/
					hypothesisLattice.backwardEdgeMerge();
						
					// get the transcription
					const char *strTranscription = trnFile.getTranscription(strUtteranceId);
					if (strTranscription == NULL) {
						BVC_ERROR << "no transcription for utterance: \"" << strUtteranceId << "\" was found";
					}
					
					// extract lexical units from the transcription
					VLexUnit vLexUnitsTranscription;
					bool bAllKnown;
					lexiconManager.getLexUnits(strTranscription,vLexUnitsTranscription,bAllKnown);
					//lexiconManager->print(vLexUnitsTranscription);
					
		
					
					// compute the lattice WER
					LatticeWER *latticeWER = hypothesisLattice.computeWER(vLexUnitsTranscription,mappings);
					if (latticeWER == NULL) {
						BVC_ERROR << "unable to compute the WER for the utterance \"" << strUtteranceId << "\"";
					}
					HypothesisLattice::add(&latticeWERAll,latticeWER);
					
					// compute the lattice depth
					LatticeDepth *latticeDepth = hypothesisLattice.computeDepth();
					HypothesisLattice::add(&latticeDepthAll,latticeDepth);
						
					delete latticeWER;
					delete latticeDepth;
				}
				HypothesisLattice::print(&latticeWERAll);
				HypothesisLattice::print(&latticeDepthAll);
			}
			if (mappings) {
				delete mappings;
			}
		}
		// lattice alignment (and HMM-state marking)
		else if (strcmp(strAction,"align") == 0) {
		
			const char *strFileModels = commandLineManager.getParameterValue("-mod");
			float fBeamWidth = 2000;
			
			// load the acoustic models
			HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
			hmmManager.load(strFileModels);
			hmmManager.initializeDecoding();
			
			// create the aligner object
			Viterbi viterbi(&phoneSet,&hmmManager,&lexiconManager,fBeamWidth);
			
			// load the batch file
			BatchFile batchFile(strFileBatch,"latticeIn|features|latticeOut");
			batchFile.load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {	
				
				const char *strFileLatticeInput = batchFile.getField(i,"latticeIn");
				const char *strFileFeatures = batchFile.getField(i,"features");
				const char *strFileLatticeOutput = batchFile.getField(i,"latticeOut");
				
				printf("lattice #: %d (%s)\n",i,strFileLatticeInput);
				
				// load the lattice
				HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
				hypothesisLattice.load(strFileLatticeInput);	
				hypothesisLattice.check();
				
				// load the features
				FeatureFile featureFile(strFileFeatures,MODE_READ);
				featureFile.load();
				int iFeatures = -1;
				float *fFeatures = featureFile.getFeatureVectors(&iFeatures);
				if (iFeatures != hypothesisLattice.getFrames()) {
					BVC_ERROR << "features and lattice do not match";
				}
				
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
					
				// mark the lattice with acoustic scores, phone-boundaries and HMM-states
				hypothesisLattice.hmmMarking(&hmmManager);
				
				// mark the lattice with phone-alignments
				if (viterbi.align(fFeatures,iFeatures,&hypothesisLattice) == false) {
					BVC_ERROR << "unable to generate phone-level alignments for the lattice";
				}
				
				// write the lattice to disk
				hypothesisLattice.store(strFileLatticeOutput);
				
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice.getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice.getFrames());
					
				delete [] fFeatures;
			}
		}
		// attach language model log-likelihoods and insertion penalties
		else if (strcmp(strAction,"lm") == 0) {
		
			// load the batch file
			BatchFile batchFile(strFileBatch,"latticeIn|latticeOut");
			batchFile.load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {	
				
				const char *strFileLatticeInput = batchFile.getField(i,"latticeIn");
				const char *strFileLatticeOutput = batchFile.getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
				hypothesisLattice.load(strFileLatticeInput);
				
				// attach language model scores to the edges in the lattice
				hypothesisLattice.attachLMProbabilities(lmManager);
				
				// attach insertion penalties
				if (commandLineManager.isParameterSet("-ip")) {
					hypothesisLattice.attachInsertionPenalty(&lexiconManager);
				}
				
				// write the lattice to disk
				hypothesisLattice.store(strFileLatticeOutput);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice.getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice.getFrames());
			}
		}
		// add a path to the lattice in case it is not already there
		else if (strcmp(strAction,"addpath") == 0) {
			
			// load the batch file
			BatchFile batchFile(strFileBatch,"latticeIn|alignment|latticeOut");
			batchFile.load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {	
				
				const char *strFileLatticeInput = batchFile.getField(i,"latticeIn");
				const char *strFileAlignment = batchFile.getField(i,"alignment");
				const char *strFileLatticeOutput = batchFile.getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
				hypothesisLattice.load(strFileLatticeInput);
				
				// load the alignment
				Alignment *alignment = Alignment::load(strFileAlignment,&lexiconManager);
				assert(alignment);
				
				// get the word sequence from the alignment
				VLexUnit vLexUnits;
				VWordAlignment *vWordAlignment = alignment->getWordAlignment();
				for(VWordAlignment::iterator it = vWordAlignment->begin() ; it != vWordAlignment->end() ; ++it) {
					LexUnit *lexUnit = lexiconManager.getLexUnitPron((*it)->iLexUnitPron);
					if (lexiconManager.isStandard(lexUnit)) {
						vLexUnits.push_back(lexiconManager.getLexUnitPron((*it)->iLexUnitPron));
					}
				}
				
				// get the lattice WER taking the path to be added as the reference
				LatticeWER *latticeWER = hypothesisLattice.computeWER(vLexUnits,NULL);
				if (latticeWER == NULL) {
					BVC_ERROR << "unable to compute the WER";
				}
				// if the path is not in the lattice then we need to insert it
				if (latticeWER->iErrors != 0) {
					hypothesisLattice.addPath(alignment,true);
				}
				delete latticeWER;	
				
				// write the lattice to disk
				hypothesisLattice.store(strFileLatticeOutput);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice.getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice.getFrames());
				
				delete alignment;
			}
		}
		// confidence/posterior probabilities
		else if (strcmp(strAction,"pp") == 0) {
		
			assert(commandLineManager.isParameterSet("-ams"));
			assert(commandLineManager.isParameterSet("-lms"));
		
			// get scale factors
			float fScaleAM = atof(commandLineManager.getParameterValue("-ams"));
			float fScaleLM = atof(commandLineManager.getParameterValue("-lms"));
			
			// load the batch file
			BatchFile batchFile(strFileBatch,"latticeIn|latticeOut");
			batchFile.load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {	
				
				const char *strFileLatticeInput = batchFile.getField(i,"latticeIn");
				const char *strFileLatticeOutput = batchFile.getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
				hypothesisLattice.load(strFileLatticeInput);
				
				// attach insertion penalties
				if (commandLineManager.isParameterSet("-ip")) {
					hypothesisLattice.attachInsertionPenalty(&lexiconManager);
				}
				
				// compute forward-backward scores
				hypothesisLattice.computeForwardBackwardScores(fScaleAM,fScaleLM);
				hypothesisLattice.computePosteriorProbabilities();
				
				// compute confidence estimates?
				if (commandLineManager.isParameterSet("-conf")) {
					hypothesisLattice.computeConfidenceScore(CONFIDENCE_MEASURE_MAXIMUM);	
				}
					
				// write the lattice to disk
				hypothesisLattice.store(strFileLatticeOutput);
				
				// text format
				ostringstream strFileLatticeOutputTxt;
				strFileLatticeOutputTxt << strFileLatticeOutput << ".txt";
				hypothesisLattice.store(strFileLatticeOutputTxt.str().c_str(),FILE_FORMAT_TEXT);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice.getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice.getFrames());
			}
		}
		// compacting
		else if (strcmp(strAction,"compact") == 0) {
		
			// load the batch file
			BatchFile batchFile(strFileBatch,"latticeIn|latticeOut");
			batchFile.load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {	
				
				const char *strFileLatticeInput = batchFile.getField(i,"latticeIn");
				const char *strFileLatticeOutput = batchFile.getField(i,"latticeOut");
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
				hypothesisLattice.load(strFileLatticeInput);
				
				// forward/backward compacting
				hypothesisLattice.forwardEdgeMerge();
				hypothesisLattice.backwardEdgeMerge();	
				
				// write the lattice to disk
				hypothesisLattice.store(strFileLatticeOutput);
		
				double dTimeEnd = TimeUtils::getTimeMilliseconds();
				double dTime = (dTimeEnd-dTimeBegin)/1000.0;
				double dRTF = dTime/(hypothesisLattice.getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice.getFrames());
			}
		}
		// rescoring
		else if (strcmp(strAction,"rescore") == 0) {
		
			// get the rescoring method
			const char *strRescoringMethod = commandLineManager.getParameterValue("-res");
			
			float fScaleAM = 0.0;
			float fScaleLM = 0.0;
			if (strcmp(strRescoringMethod,RESCORING_METHOD_LIKELIHOOD) == 0) {
				
				assert(commandLineManager.isParameterSet("-ams"));
				assert(commandLineManager.isParameterSet("-lms"));	
			
				// get scale factors
				fScaleAM = atof(commandLineManager.getParameterValue("-ams"));
				fScaleLM = atof(commandLineManager.getParameterValue("-lms"));
			}
				
			// get the hypothesis format
			const char *strFileHypFormat = commandLineManager.getParameterValue("-hypf");
			const char *strFileHypothesis = NULL;
			bool bTrn = true;
			const char *strBatchType = "lattice|utteranceId";
			FileOutput *fileHyp = NULL;
			if ((!strFileHypFormat) || strcmp(strFileHypFormat,"trn") == 0) {	
				strFileHypothesis = commandLineManager.getParameterValue("-hyp");
				fileHyp = new FileOutput(strFileHypothesis,false);
				fileHyp->open();	
			} else {
				assert(strcmp(strFileHypFormat,"ctm") == 0);
				strBatchType = "lattice|utteranceId|hypothesis";
				bTrn = false;
			}
				
			// load the batch file
			BatchFile batchFile(strFileBatch,strBatchType);
			batchFile.load();
			
			// process the batch file	
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {	
				
				const char *strFileLatticeInput = batchFile.getField(i,"lattice");
				const char *strUtteranceId = batchFile.getField(i,"utteranceId");
				if (bTrn == false) {
					strFileHypothesis = batchFile.getField(i,"hypothesis");
				}
		
				double dTimeBegin = TimeUtils::getTimeMilliseconds();
				
				// load the lattice
				HypothesisLattice hypothesisLattice(&phoneSet,&lexiconManager,bVerbose);
				hypothesisLattice.load(strFileLatticeInput);
				
				// attach insertion penalties
				if (commandLineManager.isParameterSet("-ip")) {
					hypothesisLattice.attachInsertionPenalty(&lexiconManager);
				}	
				
				// likelihood based rescoring: set scaling factors
				if (strcmp(strRescoringMethod,RESCORING_METHOD_LIKELIHOOD) == 0) {	
					hypothesisLattice.setScalingFactors(fScaleAM,fScaleLM);
				} 
				
				// lattice rescoring
				BestPath *bestPath = hypothesisLattice.rescore(strRescoringMethod);
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
				double dRTF = dTime/(hypothesisLattice.getFrames()/100.0);
				
				printf("Lattice processing time: %.4fs (RTF= %5.4f) frames: %d\n",dTime,dRTF,hypothesisLattice.getFrames());
			}
			if (bTrn) {
				fileHyp->close();
				delete fileHyp;
			}	
			if (lmManager) {
				delete lmManager;
			}
		}
		// unsupported action
		else {	
			printf("action not supported\n");
			return -1;
		}
	
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}

	return 0;
}
