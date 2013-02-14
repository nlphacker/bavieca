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


#include "Global.h"

#include "AlignmentFile.h"
#include "BatchFile.h"
#include "FeatureFile.h"
#include "FeatureExtractor.h"
#include "FileUtils.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "HMMManager.h"
#include "HMMInitializer.h"
#include "MLEstimator.h"
#include "ReferenceText.h"
#include "CommandLineManager.h"
#include "AudioFile.h"
#include "LexUnitsFile.h"
#include "TimeUtils.h"
#include "ViterbiX.h"

#include <iostream>
#include <cstdlib>

using namespace std;

#include <string>
#include <map>
#include <vector>

using namespace Bavieca;

// main for the aligner tool: "aligner"
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("aligner",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-fea","feature vectors",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-txt","text to align to the audio",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");
		commandLineManager.defineParameter("-out","output alignment file",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-fof","folder containing the features to align (mlf mode)",
			PARAMETER_TYPE_FOLDER,true);	
		commandLineManager.defineParameter("-mlf","master label file containing data to align",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-dir","output directory to store the alignments (-mlf active)",
			PARAMETER_TYPE_FOLDER,true);	
		commandLineManager.defineParameter("-bat","batch file containing entries (featuresFile txtFile alignmentFile)",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-opt","file containing optional symbols that can be inserted at word-boundaries",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-pro","whether to allow multiple pronunciations",
			PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
		commandLineManager.defineParameter("-bea","beam width used for pruning",PARAMETER_TYPE_FLOAT,true,NULL,"1000.0");
		commandLineManager.defineParameter("-hlt","whether to halt the batch processing if an error is found",
			PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");
		
		// parse the command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileLexicon = commandLineManager.getParameterValue("-lex");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");	
		const char *strFileAlignment = commandLineManager.getParameterValue("-out");	
		const char *strAlignmentFormat = commandLineManager.getParameterValue("-for");
		
		// allow multiple pronunciations?
		bool bMultiplePronunciations = CommandLineManager::str2bool(commandLineManager.getParameterValue("-pro"));
		
		// insert optional symbols?
		const char *strFileOptionalSymbols = commandLineManager.getParameterValue("-opt");
	
		// get the beam width
		float fBeamWidth = atof(commandLineManager.getParameterValue("-bea"));
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet); 
		lexiconManager.load();
		
		// load the optional symbols
		VLexUnit vLexUnitOptional;
		if (strFileOptionalSymbols) {
			// get the optional symbols
			LexUnitsFile lexUnitsFile(&lexiconManager,strFileOptionalSymbols);
			lexUnitsFile.load();
			lexUnitsFile.getLexUnits(vLexUnitOptional);
		} else {
			vLexUnitOptional.push_back(lexiconManager.getLexUnitSilence());
		}
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		hmmManager.initializeDecoding();
		
		// create the aligner
		double dUtteranceLikelihood;
		int iErrorCode;
		fBeamWidth = FLT_MAX;
		ViterbiX viterbiX(&phoneSet,&lexiconManager,&hmmManager,fBeamWidth,1000,true,1000);
		
		// case 1: feature vectors vs text (features are assumed to be already normalized)
		if (commandLineManager.isParameterSet("-fea") == true) {
		
			string strFileFeatures = commandLineManager.getParameterValue("-fea");
			string strFileText = commandLineManager.getParameterValue("-txt"); 
		
			// load the features
			FeatureFile featureFile(strFileFeatures.c_str(),MODE_READ);
			featureFile.load();
			int iFeatureVectors = 0;
			float *fFeatureVectors = (float*)featureFile.getFeatureVectors(&iFeatureVectors);	
			
			// load the lexical units
			LexUnitsFile lexUnitsFile(&lexiconManager,strFileText.c_str());
			lexUnitsFile.load();
			VLexUnit *vLexUnitText = lexUnitsFile.getLexUnits();
		
			Alignment *alignment = viterbiX.processUtterance(*vLexUnitText,bMultiplePronunciations,
				vLexUnitOptional,fFeatureVectors,iFeatureVectors,&dUtteranceLikelihood,iErrorCode);	
			if (alignment == NULL) {
				BVC_ERROR << "unable to perform the alignment";
			}
			
			// write the alignment to the output file
			// binary format
			if (strcmp(strAlignmentFormat,"binary") == 0) {
				alignment->store(strFileAlignment);
			} 
			// text format
			else {
				VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(&lexiconManager);
				AlignmentFile alignmentFile(&phoneSet,&lexiconManager);
				alignmentFile.store(*vPhoneAlignment,strFileAlignment);
				AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
			}	
			
			// clean-up
			delete alignment;
			delete [] fFeatureVectors;
		} 
		// case 2: batch mode
		else if (commandLineManager.isParameterSet("-bat") == true) {
		
			// load the batch file
			string strFileBatch = commandLineManager.getParameterValue("-bat");
			BatchFile batchFile(strFileBatch.c_str(),"features|transcription|alignment");
			batchFile.load();		
			
			// check that there is at least one entry in the batch file
			if (batchFile.size() <= 0) {
				BVC_ERROR << "empty batch file";
			}		
			
			// determine whether to stop if an error is found
			bool bStopIfError = false;
			if (strcmp(commandLineManager.getParameterValue("-stp"),"yes") == 0) {
				bStopIfError = true;
			}
			
			// iterate through the batch file elements
			for(unsigned int i = 0 ; i < batchFile.size() ; ++i) {
			
				const char *strFileFeatures = batchFile.getField(i,"features");
				const char *strFileText = batchFile.getField(i,"transcription");
				const char *strFileAlignment = batchFile.getField(i,"alignment");
			
				// load the features
				FeatureFile featureFile(strFileFeatures,MODE_READ);
				try {
					featureFile.load();
				} catch (ExceptionBase &e) {
					if (bStopIfError) {
						BVC_ERROR << "unable to load the feature file " << strFileFeatures;
					} else {
						BVC_WARNING << "unable to load the feature file " << strFileFeatures;
						continue;
					}	
				}			
				int iFeatureVectors = 0;
				float *fFeatureVectors = featureFile.getFeatureVectors(&iFeatureVectors);	
						
				// load the lexical units
				LexUnitsFile lexUnitsFile(&lexiconManager,strFileText);
				try {
					lexUnitsFile.load();
				} catch (ExceptionBase &e) {
					if (bStopIfError) {
						BVC_ERROR << "unable to load the text file " << strFileText;
					} else {
						BVC_WARNING << "unable to load the text file " << strFileText;
						continue;
					}	
				}			
				VLexUnit *vLexUnitText = lexUnitsFile.getLexUnits();
				
				Alignment *alignment = viterbiX.processUtterance(*vLexUnitText,bMultiplePronunciations,
					vLexUnitOptional,fFeatureVectors,iFeatureVectors,&dUtteranceLikelihood,iErrorCode);				
				if (!alignment) {
					if (bStopIfError) {
						BVC_ERROR << "aligning features file: " << strFileFeatures << " to the corresponding transcription";
					} else {
						BVC_WARNING << "aligning features file: " << strFileFeatures << " to the corresponding transcription";
						delete [] fFeatureVectors;
						continue;
					}
				}
						
				// write the alignment into the output file
				try {
					// binary format
					if (strcmp(strAlignmentFormat,"binary") == 0) {
						alignment->store(strFileAlignment);
					}
					// text format
					else {
						VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(&lexiconManager);
						AlignmentFile alignmentFile(&phoneSet,&lexiconManager);
						alignmentFile.store(*vPhoneAlignment,strFileAlignment);
						AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);	
					}	
				} catch (ExceptionBase &e) {
					if (bStopIfError) {
						BVC_ERROR << "unable to create the alignment file: " << strFileAlignment;
					} else {
						BVC_WARNING << "unable to create the alignment file: " << strFileAlignment;
						delete alignment;
						delete [] fFeatureVectors;
						continue;
					}	
				}
				// clean-up
				delete alignment;	
				delete [] fFeatureVectors;	
			}
		}
		// case 3: master label file
		else {
			
			assert(commandLineManager.isParameterSet("-mlf") == true);
		
			const char *strFileMLF = commandLineManager.getParameterValue("-mlf");
		
			// get the output directory (base directory)
			const char *strDirectoryOutput = NULL;
			if (commandLineManager.isParameterSet("-dir") == true) {
				strDirectoryOutput = commandLineManager.getParameterValue("-dir");
			}
			
			// load the master label file
			MLFFile mlfFile(&lexiconManager,strFileMLF,MODE_READ);
			mlfFile.load();
			
			VMLFUtterance *vMLFUtterance = mlfFile.getUtterances();
			
			// determine whether to stop if an error is found
			bool bStopIfError = false;
			if (strcmp(commandLineManager.getParameterValue("-stp"),"yes") == 0) {
				bStopIfError = true;
			}
			
			// get the folder containing the features
			const char *strFolderFeatures = commandLineManager.getParameterValue("-fof");
			char strFileFeatures[1024];
			char strFileAlignment[1024];
			
			// iterate through all the utterances in the Master Label File
			for(VMLFUtterance::iterator it = vMLFUtterance->begin() ; it != vMLFUtterance->end() ; ++it) {	
			
				// create the feature file
				sprintf(strFileFeatures,"%s%c%s",strFolderFeatures,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
				
				// build the alignment filename
				sprintf(strFileAlignment,"%s%c%s",strDirectoryOutput,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
				FileUtils::replaceExtension(strFileAlignment,strFileAlignment,"ali");
				// get the folder name and create it if necessary
				char strPathAlignment[1024];
				FileUtils::getFolder(strPathAlignment,strFileAlignment);
				if (FileUtils::createPath(strPathAlignment) != RETURN_CODE_SUCCESS) {
					if (bStopIfError) {
						BVC_ERROR << "unable to create the path: " << strPathAlignment;
					} else {
						BVC_WARNING << "unable to create the path: " << strPathAlignment;
						continue;
					}
				}
				
				// create the output directory if it does not exist
			
				// load the features
				FeatureFile featureFile(strFileFeatures,MODE_READ);
				try {
					featureFile.load();
				} catch (ExceptionBase &e) {
					if (bStopIfError) {
						BVC_ERROR << "unable to load the features file: " << strFileFeatures; 
					} else {
						BVC_WARNING << "unable to load the features file: " << strFileFeatures; 
						continue;
					}	
				}
				int iFeatureVectors = 0;
				float *fFeatureVectors = (float*)featureFile.getFeatureVectors(&iFeatureVectors);	
				
				Alignment *alignment = viterbiX.processUtterance((*it)->vLexUnit,bMultiplePronunciations,
					vLexUnitOptional,fFeatureVectors,iFeatureVectors,&dUtteranceLikelihood,iErrorCode);				
				if (!alignment) {
					if (bStopIfError) {
						BVC_ERROR << "aligning features file: " << strFileFeatures 
							<< " to the corresponding transcription in the MLF";
					} else {
						BVC_WARNING << "aligning features file: " << strFileFeatures 
							<< " to the corresponding transcription in the MLF";
						delete [] fFeatureVectors;
						continue;
					}
				}
						
				// write the alignment into the output file
				try {
					// binary format
					if (strcmp(strAlignmentFormat,"binary") == 0) {
						alignment->store(strFileAlignment);
					}
					// text format
					else {
						VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(&lexiconManager);
						AlignmentFile alignmentFile(&phoneSet,&lexiconManager);
						alignmentFile.store(*vPhoneAlignment,strFileAlignment);
						AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);	
					}
				} catch (ExceptionBase &e) {
					if (bStopIfError) {
						BVC_ERROR << "unable to create the alignment file: " << strFileAlignment;
					} else {
						BVC_WARNING << "unable to create the alignment file: " << strFileAlignment;
						delete alignment;
						delete [] fFeatureVectors;
						continue;
					}	
				}
				// clean-up
				delete alignment;	
				delete [] fFeatureVectors;	
			}
		}
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;
}

