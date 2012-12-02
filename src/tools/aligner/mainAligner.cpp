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

	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("aligner",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-lex","pronunciation dictionary (lexicon)",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-fea","feature vectors",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-txt","text to align to the audio",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-for","alignment file format",PARAMETER_TYPE_STRING,false,"binary|text");
 	m_commandLineManager->defineParameter("-out","output alignment file",PARAMETER_TYPE_FILE,true);
 	m_commandLineManager->defineParameter("-fof","folder containing the features to align (mlf mode)",
 		PARAMETER_TYPE_FOLDER,true);	
 	m_commandLineManager->defineParameter("-mlf","master label file containing data to align",PARAMETER_TYPE_FILE,true);
 	m_commandLineManager->defineParameter("-dir","output directory to store the alignments (-mlf active)",
 		PARAMETER_TYPE_FOLDER,true);	
	m_commandLineManager->defineParameter("-bat","batch file containing entries (featuresFile txtFile alignmentFile)",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-opt","file containing optional symbols that can be inserted at word-boundaries",PARAMETER_TYPE_FILE,true);
	m_commandLineManager->defineParameter("-pro","whether to allow multiple pronunciations",
		PARAMETER_TYPE_BOOLEAN,true,NULL,"no");
	m_commandLineManager->defineParameter("-bea","beam width used for pruning",PARAMETER_TYPE_FLOAT,true,NULL,"1000.0");
	m_commandLineManager->defineParameter("-hlt","whether to halt the batch processing if an error is found",
		PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");
	
	// parse the command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileLexicon = m_commandLineManager->getParameterValue("-lex");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");	
	const char *m_strFileAlignment = m_commandLineManager->getParameterValue("-out");	
	const char *m_strAlignmentFormat = m_commandLineManager->getParameterValue("-for");
	
	// allow multiple pronunciations?
	bool m_bMultiplePronunciations = CommandLineManager::str2bool(m_commandLineManager->getParameterValue("-pro"));
	
	// insert optional symbols?
	const char *m_strFileOptionalSymbols = m_commandLineManager->getParameterValue("-opt");

	// get the beam width
	float m_fBeamWidth = atof(m_commandLineManager->getParameterValue("-bea"));
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
   
   // load the lexicon
   LexiconManager *m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
   m_lexiconManager->load();
   
   // load the optional symbols
	VLexUnit m_vLexUnitOptional;
	if (m_strFileOptionalSymbols) {
		// get the optional symbols
		LexUnitsFile lexUnitsFile(m_lexiconManager,m_strFileOptionalSymbols);
		lexUnitsFile.load();
		lexUnitsFile.getLexUnits(m_vLexUnitOptional);
	} else {
		m_vLexUnitOptional.push_back(m_lexiconManager->getLexUnitSilence());
	}
	
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeDecoding();
	
	// create the aligner
	double dUtteranceLikelihood;
	int iErrorCode;
	m_fBeamWidth = FLT_MAX;
	ViterbiX *m_viterbiX = new ViterbiX(m_phoneSet,m_lexiconManager,m_hmmManager,m_fBeamWidth,1000,true,1000);
	
   // case 1: feature vectors vs text (features are assumed to be already normalized)
	if (m_commandLineManager->isParameterSet("-fea") == true) {
	
		string m_strFileFeatures = m_commandLineManager->getParameterValue("-fea");
		string m_strFileText = m_commandLineManager->getParameterValue("-txt"); 
	
		// load the features
		FeatureFile *featureFile = new FeatureFile(m_strFileFeatures.c_str(),MODE_READ);
		featureFile->load();
		int iFeatureVectors = 0;
		float *fFeatureVectors = (float*)featureFile->getFeatureVectors(&iFeatureVectors);	
		delete featureFile;
		
		// load the lexical units
		LexUnitsFile *lexUnitsFile = new LexUnitsFile(m_lexiconManager,m_strFileText.c_str());
		lexUnitsFile->load();
		VLexUnit *vLexUnitText = lexUnitsFile->getLexUnits();
	
		Alignment *alignment = m_viterbiX->processUtterance(*vLexUnitText,m_bMultiplePronunciations,
			m_vLexUnitOptional,fFeatureVectors,iFeatureVectors,&dUtteranceLikelihood,iErrorCode);	
		if (alignment == NULL) {
			BVC_ERROR << "unable to perform the alignment";
		}
		
		// write the alignment to the output file
		// binary format
		if (strcmp(m_strAlignmentFormat,"binary") == 0) {
			alignment->store(m_strFileAlignment);
		} 
		// text format
		else {
			VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(m_lexiconManager);
			AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,m_lexiconManager);
			alignmentFile->store(*vPhoneAlignment,m_strFileAlignment);
			delete alignmentFile;
			AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
		}	
		
		// clean-up
		delete alignment;
		delete lexUnitsFile;
   	delete [] fFeatureVectors;
	} 
	// case 2: batch mode
	else if (m_commandLineManager->isParameterSet("-bat") == true) {
	
		// load the batch file
		string m_strFileBatch = m_commandLineManager->getParameterValue("-bat");
		BatchFile batchFile(m_strFileBatch.c_str(),"features|transcription|alignment");
		batchFile.load();		
		
		// check that there is at least one entry in the batch file
		if (batchFile.size() <= 0) {
			printf("Error: no entries found in the batch file\n");
			return -1;		
		}		
		
		// determine whether to stop if an error is found
		bool m_bStopIfError = false;
		if (strcmp(m_commandLineManager->getParameterValue("-stp"),"yes") == 0) {
			m_bStopIfError = true;
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
				if (m_bStopIfError) {
					return -1;
				} else {
					continue;
				}	
			}			
			int iFeatureVectors = 0;
			float *fFeatureVectors = featureFile.getFeatureVectors(&iFeatureVectors);	
					
			// load the lexical units
			LexUnitsFile lexUnitsFile(m_lexiconManager,strFileText);
			try {
				lexUnitsFile.load();
			} catch (ExceptionBase &e) {
				if (m_bStopIfError) {
					return -1;
				} else {
					continue;
				}	
			}			
			VLexUnit *vLexUnitText = lexUnitsFile.getLexUnits();
			
			Alignment *alignment = m_viterbiX->processUtterance(*vLexUnitText,m_bMultiplePronunciations,
				m_vLexUnitOptional,fFeatureVectors,iFeatureVectors,&dUtteranceLikelihood,iErrorCode);				
			if (!alignment) {
				printf("Error aligning features file: %s to the corresponding transcription\n",strFileFeatures);
				if (m_bStopIfError == true) {
					return -1;
				} else {
					delete [] fFeatureVectors;
					continue;
				}
			}
					
			// write the alignment into the output file
			try {
				// binary format
				if (strcmp(m_strAlignmentFormat,"binary") == 0) {
					alignment->store(strFileAlignment);
				}
				// text format
				else {
					VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(m_lexiconManager);
					AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,m_lexiconManager);
					alignmentFile->store(*vPhoneAlignment,strFileAlignment);
					delete alignmentFile;
					AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);	
				}	
			} catch (ExceptionBase &e) {
				if (m_bStopIfError == true) {
					return -1;
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
		
		assert(m_commandLineManager->isParameterSet("-mlf") == true);
	
		const char *m_strFileMLF = m_commandLineManager->getParameterValue("-mlf");
	
		// get the output directory (base directory)
		const char *m_strDirectoryOutput = NULL;
		if (m_commandLineManager->isParameterSet("-dir") == true) {
			m_strDirectoryOutput = m_commandLineManager->getParameterValue("-dir");
		}
		
		// load the master label file
		MLFFile *m_mlfFile = new MLFFile(m_lexiconManager,m_strFileMLF,MODE_READ);
		m_mlfFile->load();
		
		VMLFUtterance *m_vMLFUtterance = m_mlfFile->getUtterances();
		
		// determine whether to stop if an error is found
		bool m_bStopIfError = false;
		if (strcmp(m_commandLineManager->getParameterValue("-stp"),"yes") == 0) {
			m_bStopIfError = true;
		}
		
		// get the folder containing the features
		const char *strFolderFeatures = m_commandLineManager->getParameterValue("-fof");
		char strFileFeatures[1024];
		char strFileAlignment[1024];
		
		// iterate through all the utterances in the Master Label File
		for(VMLFUtterance::iterator it = m_vMLFUtterance->begin() ; it != m_vMLFUtterance->end() ; ++it) {	
		
			// create the feature file
			sprintf(strFileFeatures,"%s%c%s",strFolderFeatures,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
			
			// build the alignment filename
			sprintf(strFileAlignment,"%s%c%s",m_strDirectoryOutput,PATH_SEPARATOR,(*it)->strFilePattern.c_str());
			FileUtils::replaceExtension(strFileAlignment,strFileAlignment,"ali");
			// get the folder name and create it if necessary
			char strPathAlignment[1024];
			FileUtils::getFolder(strPathAlignment,strFileAlignment);
			if (FileUtils::createPath(strPathAlignment) != RETURN_CODE_SUCCESS) {
				printf("Error creating the path: %s\n",strPathAlignment);
				if (m_bStopIfError) {
					return -1;
				} else {
					continue;
				}
			}
			
			// create the output directory if it does not exist
		
			// load the features
			FeatureFile *featureFile = new FeatureFile(strFileFeatures,MODE_READ);
			try {
				featureFile->load();
			} catch (ExceptionBase &e) {
				if (m_bStopIfError) {
					return -1;
				} else {
					continue;
				}	
			}
			int iFeatureVectors = 0;
			float *fFeatureVectors = (float*)featureFile->getFeatureVectors(&iFeatureVectors);	
			delete featureFile;
			
			Alignment *alignment = m_viterbiX->processUtterance((*it)->vLexUnit,m_bMultiplePronunciations,
				m_vLexUnitOptional,fFeatureVectors,iFeatureVectors,&dUtteranceLikelihood,iErrorCode);				
			if (!alignment) {
				printf("Error aligning features file: %s to the corresponding transcription in the MLF\n",strFileFeatures);
				if (m_bStopIfError == true) {
					return -1;
				} else {
					delete [] fFeatureVectors;
					continue;
				}
			}
					
			// write the alignment into the output file
			try {
				// binary format
				if (strcmp(m_strAlignmentFormat,"binary") == 0) {
					alignment->store(strFileAlignment);
				}
				// text format
				else {
					VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(m_lexiconManager);
					AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,m_lexiconManager);
					alignmentFile->store(*vPhoneAlignment,strFileAlignment);
					delete alignmentFile;
					AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);	
				}	
			} catch (ExceptionBase &e) {
				if (m_bStopIfError == true) {
					return -1;
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

		// clean-up
		delete m_mlfFile;
	}

	// clean-up	
	delete m_viterbiX;
	delete m_hmmManager;
	delete m_lexiconManager;
	delete m_phoneSet;
	delete m_commandLineManager;
	
	return 0;
}

