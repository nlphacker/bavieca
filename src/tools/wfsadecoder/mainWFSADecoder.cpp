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


#include "AudioFile.h"
#include "BatchFile.h"
#include "BestPath.h"
#include "ConfigurationFeatures.h"
#include "ConfigurationWFSADecoder.h"
#include "CommandLineManager.h"
#include "FeatureExtractor.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "LogMessage.h"
#include "SADModule.h"
#include "Transform.h"
#include "TimeUtils.h"
#include "Viterbi.h"
#include "WFSADecoder.h"

using namespace Bavieca;

// main for the tool "wfsdecoder"
int main(int argc, char *argv[]) {

	try {

		// define command line parameters
		CommandLineManager commandLineManager("wfsadecoder",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-cfg","configuration",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-bat","batch file containing entries [rawFile/featureFile utteranceId]",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineIncompatibility("-list|-audio",true);
		commandLineManager.defineParameter("-hyp","hypothesis",PARAMETER_TYPE_FILE,false);
		
		// parse command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get the command-line parameters
		const char *strFileConfiguration = commandLineManager.getParameterValue("-cfg");
		
		// load the configuration file
		ConfigurationWFSADecoder *configuration = new ConfigurationWFSADecoder(strFileConfiguration);
		configuration->load();
		
		// get the configuration parameters
		
		// decoding network
		const char *strFileNetwork = configuration->getStrParameterValue("decodingNetwork.file");
		
		// input
		const char *strInputType = configuration->getStrParameterValue("input.type");	
	
		// phone-set
		const char *strFilePhoneticSymbolSet = configuration->getStrParameterValue("phoneticSymbolSet.file");	
		
		// feature extraction
		const char *strFileConfigurationFeatures = 
			configuration->getStrParameterValue("feature.configurationFile");
		const char *strCepstralNormalizationMode =
			configuration->getStrParameterValue("feature.cepstralNormalization.mode");
		const char *strCepstralNormalizationMethod =
			configuration->getStrParameterValue("feature.cepstralNormalization.method");
		int iCepstralNormalizationBufferSize = 
			configuration->getIntParameterValue("feature.cepstralNormalization.bufferSize");
		float fWarpFactor = 1.0;
		if (configuration->isParameterSet("feature.warpFactor")) {
			fWarpFactor = atof(configuration->getParameterValue("feature.warpFactor"));
		}
		const char *strFileFeatureTransform = NULL;
		if (configuration->isParameterSet("feature.transformFile")) {
			strFileFeatureTransform = configuration->getParameterValue("feature.transformFile");
		}
		
		// acoustic models
		const char *strFileModels = configuration->getStrParameterValue("acousticModels.file");
	
		// lexicon	
		const char *strFileLexicon = configuration->getStrParameterValue("lexicon.file");
		
		// pruning
		int iMaxActiveStates = configuration->getIntParameterValue("pruning.maxActiveStates");
		float fLikelihoodBeam = configuration->getFloatParameterValue("pruning.likelihoodBeam");	
		
		// output lattice?
		bool bLatticeGeneration = configuration->isParameterSet("output.lattice.folder");	
		const char *strFolderLattices = NULL;
		int iMaxWordSequencesState = -1;
		if (bLatticeGeneration) {
			strFolderLattices = configuration->getStrParameterValue("output.lattice.folder"); 
			iMaxWordSequencesState = 
				configuration->getIntParameterValue("output.lattice.maxWordSequencesState");
		}
		
		// output features?
		bool bOutputFeatures = configuration->isParameterSet("output.features.folder");
		const char *strFolderFeatures = NULL;
		if (bOutputFeatures) {
			strFolderFeatures = configuration->getStrParameterValue("output.features.folder"); ;
		}
		
		// output alignment?
		bool bOutputAlignment = configuration->isParameterSet("output.alignment.folder");
		const char *strFolderAlignments = NULL;
		if (bOutputAlignment) {
			strFolderAlignments = configuration->getStrParameterValue("output.alignment.folder"); ;
		}
	
		// output audio?
		bool bOutputAudio = configuration->isParameterSet("output.audio.folder");
		const char *strFolderAudio = NULL;
		if (bOutputAudio) {
			strFolderAudio = configuration->getStrParameterValue("output.audio.folder");
		}
		
		// (3) initialize the system
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneticSymbolSet);
		phoneSet.load();
		
		// load the feature transforms
		VTransform vTransformFeatures;
		if (strFileFeatureTransform) {
			BatchFile batchFile(strFileFeatureTransform,"transform");
			batchFile.load();
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {
				Transform *transform = new Transform();
				transform->load(batchFile.getField(i,0u));
				vTransformFeatures.push_back(transform);
			}	
		}   
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		hmmManager.initializeDecoding();	
		
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet);
		lexiconManager.load();
		lexiconManager.print(false);
		
		// create the aligner object (if needed)
		Viterbi *viterbi = NULL;
		if (bOutputAlignment == true) {
			viterbi = new Viterbi(&phoneSet,&hmmManager,&lexiconManager);
		}
		
		// load the feature configuration
		ConfigurationFeatures configurationFeatures(strFileConfigurationFeatures);
		configurationFeatures.load();
		
		// get the feature normalization mode and method
		int iCepstralNormalizationMode = FeatureExtractor::getNormalizationMode(strCepstralNormalizationMode);
		int iCepstralNormalizationMethod = FeatureExtractor::getNormalizationMethod(strCepstralNormalizationMethod);
		
		// create the feature extractor
		FeatureExtractor featureExtractor(&configurationFeatures,fWarpFactor,
			iCepstralNormalizationBufferSize,iCepstralNormalizationMode,iCepstralNormalizationMethod);
		featureExtractor.initialize();
	
		// load the decoding network
		WFSAcceptor *wfsAcceptor = WFSAcceptor::load(&lexiconManager,strFileNetwork);
		assert(wfsAcceptor);
	
		// create the decoding object
		WFSADecoder wfsaDecoder(&phoneSet,&hmmManager,&lexiconManager,
			wfsAcceptor,iMaxActiveStates,fLikelihoodBeam,bLatticeGeneration,iMaxWordSequencesState);
		wfsaDecoder.initialize();
		
		string strFileHypothesis = commandLineManager.getParameterValue("-hyp");
		
		double dLikelihoodTotal = 0.0;
		unsigned int iFeatureVectorsTotal = 0;
		unsigned int iUtterances = 0;
		double dTimeBegin = TimeUtils::getTimeMilliseconds();
		
		// BATCH MODE (batch file containing pairs [ featureFile  utteranceID ])
		if (commandLineManager.isParameterSet("-bat")) {
		
			// best single path options
			bool bBestSinglePathOutputFillers = false;
			bool bBestSinglePathOutputConfidenceValues = false;
			bool bBestSinglePathOutputSentenceDelimiters = false;
		
			string strFileList = commandLineManager.getParameterValue("-bat");
			VFeaturesUtterance vFeaturesUtterance;
			BatchFile *batchFile = NULL;
			
			// audio input
			if (strcmp(strInputType,"audio") == 0) {
			
				// load the batch file
				batchFile = new BatchFile(strFileList.c_str(),"audio|id");
				batchFile->load();
				
				// (4) extract session-data
				VUtteranceData vUtteranceData;
				for(unsigned int iUtterance = 0 ; iUtterance < batchFile->size() ; ++iUtterance) {
				
					// load the raw audio
					int iSamples = -1;
					short int *sSamples = AudioFile::load(batchFile->getField(iUtterance,"audio"),&iSamples);
				
					UtteranceData utteranceData;
					utteranceData.samples.sSamples = sSamples;
					utteranceData.samples.iSamples = iSamples;
					utteranceData.features.fFeatures = NULL;
					utteranceData.features.iFeatures = -1;
					vUtteranceData.push_back(utteranceData);
				}	
				
				// extract features
				featureExtractor.extractFeaturesSession(vUtteranceData,true);
				
				// apply feture transforms
				int iDimFea = featureExtractor.getFeatureDimensionality();
				for(VTransform::iterator it = vTransformFeatures.begin() ; it != vTransformFeatures.end() ; ++it) {
					printf("%d -> %d\n",iDimFea,(*it)->getRows());
					for(VUtteranceData::iterator jt = vUtteranceData.begin() ; jt != vUtteranceData.end() ; ++jt) {
						float *fFeaturesX = new float[jt->features.iFeatures*(*it)->getRows()];
						for(int i=0 ; i < jt->features.iFeatures ; ++i) {
							(*it)->apply(jt->features.fFeatures+(i*iDimFea),fFeaturesX+(i*(*it)->getRows()));
						}
						delete [] jt->features.fFeatures;
						jt->features.fFeatures = fFeaturesX;
					}
					iDimFea = (*it)->getRows();
				}
				
				FeaturesUtterance featuresUtterance;
				for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it) {
					featuresUtterance.fFeatures = it->features.fFeatures;
					featuresUtterance.iFeatures = it->features.iFeatures;
					delete [] it->samples.sSamples;
					vFeaturesUtterance.push_back(featuresUtterance);	
				}
			}
			// feature input
			else {
			
				// load the batch file
				batchFile = new BatchFile(strFileList.c_str(),"features|id");
				batchFile->load();
			
				iUtterances = batchFile->size();
				for(unsigned int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) {
				
					string strUtteranceId = batchFile->getField(iUtterance,"id");
					string strFile = batchFile->getField(iUtterance,"features");
				
					FeatureFile featureFile(strFile.c_str(),MODE_READ);
					featureFile.load();
					FeaturesUtterance featuresUtterance;
					featuresUtterance.fFeatures = featureFile.getFeatureVectors(&featuresUtterance.iFeatures);
					vFeaturesUtterance.push_back(featuresUtterance);	
				}
			}
			
			// actual decoding
			FileOutput fileHypothesis(strFileHypothesis.c_str(),false);
			fileHypothesis.open();
			iUtterances = batchFile->size();
			for(unsigned int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) {
			
				string strUtteranceId = batchFile->getField(iUtterance,"id");
				printf("Processing utterance: %s\n",strUtteranceId.c_str());
				
				int iFeatures = vFeaturesUtterance[iUtterance].iFeatures;
				float *fFeatures = vFeaturesUtterance[iUtterance].fFeatures;;
				wfsaDecoder.viterbi(fFeatures,iFeatures);
				iFeatureVectorsTotal += iFeatures;
				
				// hypothesis lattice
				if (bLatticeGeneration) {
					HypothesisLattice *hypothesisLattice = wfsaDecoder.getHypothesisLattice();
					if (hypothesisLattice) {
						ostringstream ossText,ossBin;
						ossText << strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".txt";
						hypothesisLattice->store(ossText.str().c_str(),FILE_FORMAT_TEXT);
						ossBin << strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".bin";
						hypothesisLattice->store(ossBin.str().c_str(),FILE_FORMAT_BINARY);
						delete hypothesisLattice;	
					}	
				}
				
				// best path
				BestPath *bestPathUtterance = wfsaDecoder.getBestPath();		
				if (bestPathUtterance) {
					dLikelihoodTotal += bestPathUtterance->getPathScore();
					bestPathUtterance->print();
					bestPathUtterance->write(fileHypothesis.getStream(),strUtteranceId.c_str(),
						bBestSinglePathOutputSentenceDelimiters,bBestSinglePathOutputFillers,bBestSinglePathOutputConfidenceValues);
					// state-alignment
					if (bOutputAlignment) {
						VPhoneAlignment *vPhoneAlignment = viterbi->align(fFeatures,iFeatures,bestPathUtterance);
						if (vPhoneAlignment) {
							AlignmentFile alignmentFile(&phoneSet,&lexiconManager);
							char strFileAlignment[1024+1];
							sprintf(strFileAlignment,"%s%c%s.ali",strFolderAlignments,PATH_SEPARATOR,strUtteranceId.c_str());
							alignmentFile.store(*vPhoneAlignment,strFileAlignment);
							AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
						}	
					}
					delete bestPathUtterance;
				}
				
				// clean-up
				if (fFeatures != NULL) {
					delete [] fFeatures;	
				}
			}	
			fileHypothesis.close();
			delete batchFile;
		} 
		
		double dTimeEnd = TimeUtils::getTimeMilliseconds();
		double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
		double dRTF = dTimeSeconds/(((float)iFeatureVectorsTotal)/100.0); 
		
		printf("- summary ------------------------------------\n");
		printf("# utterances: %d speech time: %.2f seconds\n",iUtterances,((float)iFeatureVectorsTotal)/100.0);
		printf("decoding time: %.2f seconds (RTF: %5.2f)\n",dTimeSeconds,dRTF);
		printf("likelihood: %.4f (per frame: %8.4f)\n",dLikelihoodTotal,dLikelihoodTotal/((float)iFeatureVectorsTotal));
		printf("----------------------------------------------\n");	
	
		// clean up
		if (viterbi) {
			delete viterbi;
		}
		delete wfsAcceptor;
		
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	

	return 0;
}

