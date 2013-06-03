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
#include "LMManager.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

using namespace std;

#include <string>

using namespace Bavieca;

// main for the tool "dynamicdecoder"
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("dynamicdecoder",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);	
		commandLineManager.defineParameter("-cfg","configuration file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-hyp","hypothesis file",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-bat","batch file with entries [rawFile/featureFile utteranceId]",
			PARAMETER_TYPE_FILE,false);
			
		// (2) process command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get command line parameters
		const char *strFileConfiguration = commandLineManager.getParameterValue("-cfg");
		const char *strFileControl = commandLineManager.getParameterValue("-bat");
		const char *strFileHypothesis = commandLineManager.getParameterValue("-hyp");
		
		// load the configuration file
		ConfigurationDynamicDecoder configuration(strFileConfiguration);
		configuration.load();
		
		// configuration.ConfigurationFile::print();
		
		// get parameters from the configuration file
	
		// feature extraction
		const char *strFileConfigurationFeatures = 
			configuration.getStrParameterValue("feature.configurationFile");
		const char *strCepstralNormalizationMode =
			configuration.getStrParameterValue("feature.cepstralNormalization.mode");
		const char *strCepstralNormalizationMethod =
			configuration.getStrParameterValue("feature.cepstralNormalization.method");
		int iCepstralNormalizationBufferSize = 
			configuration.getIntParameterValue("feature.cepstralNormalization.bufferSize");
		float fWarpFactor = 1.0;
		if (configuration.isParameterSet("feature.warpFactor")) {
			fWarpFactor = atof(configuration.getParameterValue("feature.warpFactor"));
		}
		const char *strFileFeatureTransform = NULL;
		if (configuration.isParameterSet("feature.transformFile")) {
			strFileFeatureTransform = configuration.getParameterValue("feature.transformFile");
		}
			
		// phone set
		const char *strFilePhoneSet = 
			configuration.getStrParameterValue("phoneticSymbolSet.file");
		
		// acoustic models
		const char *strFileAcousticModels = 
			configuration.getStrParameterValue("acousticModels.file"); 
		
		// language model
		const char *strLanguageModelFile = 
			configuration.getStrParameterValue("languageModel.file"); 
		const char *strLanguageModelFormat = 
			configuration.getStrParameterValue("languageModel.format"); 
		const char *strLanguageModelType = 
			configuration.getStrParameterValue("languageModel.type"); 
		float fLanguageModelScalingFactor = 
			configuration.getFloatParameterValue("languageModel.scalingFactor"); 
		//bool bLanguageCrossUtterance = 
		//	configuration.getBoolParameterValue("languageModel.crossUtterance");
			
		// lexicon
		const char *strFileLexicon = 
			configuration.getStrParameterValue("lexicon.file"); 
			
		// insertion penalty
		float fInsertionPenaltyStandard = 
			configuration.getFloatParameterValue("insertionPenalty.standard"); 
		float fInsertionPenaltyFiller = 
			configuration.getFloatParameterValue("insertionPenalty.filler"); 
		const char *strFileInsertionPenaltyFiller = 
			configuration.getStrParameterValue("insertionPenalty.filler.file"); 
		
		// pruning parameters
		int iMaxActiveArcs = configuration.getIntParameterValue("pruning.maxActiveArcs");
		int iMaxActiveArcsWE = configuration.getIntParameterValue("pruning.maxActiveArcsWE");
		int iMaxActiveTokensArc = configuration.getIntParameterValue("pruning.maxActiveTokensArc");	
		float fBeamWidthArcs = configuration.getFloatParameterValue("pruning.likelihoodBeam");
		float fBeamWidthArcsWE = configuration.getFloatParameterValue("pruning.likelihoodBeamWE");
		float fBeamWidthTokensArc = configuration.getFloatParameterValue("pruning.likelihoodBeamTokensArc");
		
		// output lattice?
		bool bLatticeGeneration = configuration.isParameterSet("output.lattice.folder");	
		const char *strFolderLattices = NULL;
		int iMaxWordSequencesState = -1;
		if (bLatticeGeneration) {
			strFolderLattices = configuration.getStrParameterValue("output.lattice.folder"); 
			iMaxWordSequencesState = 
				configuration.getIntParameterValue("output.lattice.maxWordSequencesState");
		}
		
		// output features?
		bool bOutputFeatures = configuration.isParameterSet("output.features.folder");
		const char *strFolderFeatures = NULL;
		if (bOutputFeatures) {
			strFolderFeatures = configuration.getStrParameterValue("output.features.folder"); ;
		}
		
		// output alignment?
		bool bOutputAlignment = configuration.isParameterSet("output.alignment.folder");
		const char *strFolderAlignments = NULL;
		if (bOutputAlignment) {
			strFolderAlignments = configuration.getStrParameterValue("output.alignment.folder"); ;
		}
	
		// output audio?
		bool bOutputAudio = configuration.isParameterSet("output.audio.folder");
		const char *strFolderAudio = NULL;
		if (bOutputAudio) {
			strFolderAudio = configuration.getStrParameterValue("output.audio.folder");
		}
	
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
	
		// load the lexicon
		LexiconManager lexiconManager(strFileLexicon,&phoneSet); 
		lexiconManager.load();
		lexiconManager.attachLexUnitPenalties(fInsertionPenaltyStandard,fInsertionPenaltyFiller);
		
		FillerManager fillerManager(strFileInsertionPenaltyFiller);	
		fillerManager.load();
		fillerManager.attachInsertionPenaltyFillers(&lexiconManager);
		
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
			
		// load the HMMs used for the estimation
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileAcousticModels);
		hmmManager.initializeDecoding();
		
		// create the aligner object?
		Viterbi *viterbi = NULL;
		if (bOutputAlignment) {
			viterbi = new Viterbi(&phoneSet,&hmmManager,&lexiconManager);	
		}
		
		// load the language model
		LMManager lmManager(&lexiconManager,
									strLanguageModelFile,
									strLanguageModelFormat,
									strLanguageModelType); 
		lmManager.load();
		
		NetworkBuilderX networkBuilder(&phoneSet,&hmmManager,&lexiconManager);
		
		// build the decoding network
		DynamicNetworkX *network = networkBuilder.build();
		if (!network) {
			BVC_ERROR << "unable to build the decoding network";
		}
	
		DynamicDecoderX decoder(&phoneSet,&hmmManager,&lexiconManager,
				&lmManager,fLanguageModelScalingFactor,network,iMaxActiveArcs,
				iMaxActiveArcsWE,iMaxActiveTokensArc,fBeamWidthArcs,fBeamWidthArcsWE,fBeamWidthTokensArc,
				bLatticeGeneration,iMaxWordSequencesState);
	
		// initialize the decoder
		decoder.initialize();
		
		double dTimeBegin = TimeUtils::getTimeMilliseconds();	
		
		double dLikelihoodTotal = 0.0;
		int iFeatureVectorsTotal = 0;
		
		// load the batch file
		BatchFile batchFile(strFileControl,"audio|id");
		batchFile.load();
		
		// (4) extract session-data
		VUtteranceData vUtteranceData;
		for(unsigned int iUtterance = 0 ; iUtterance < batchFile.size() ; ++iUtterance) {
		
			// load the raw audio
			int iSamples = -1;
			short int *sSamples = AudioFile::load(batchFile.getField(iUtterance,"audio"),&iSamples);
		
			UtteranceData utteranceData;
			utteranceData.samples.sSamples = sSamples;
			utteranceData.samples.iSamples = iSamples;
			utteranceData.mFeatures = NULL;
			vUtteranceData.push_back(utteranceData);
		}	
		
		// extract features
		featureExtractor.extractFeaturesSession(vUtteranceData,true);
		
		// apply feature transforms
		for(VTransform::iterator it = vTransformFeatures.begin() ; it != vTransformFeatures.end() ; ++it) {
			BVC_VERB << "feature transform: (input dim: " << featureExtractor.getFeatureDim() 
				<< ") -> (output dim: " << (*it)->getRows() << ")";
			for(VUtteranceData::iterator jt = vUtteranceData.begin() ; jt != vUtteranceData.end() ; ++jt) {
				Matrix<float> *mFeaturesX = (*it)->apply(*jt->mFeatures);
				delete jt->mFeatures;
				jt->mFeatures = mFeaturesX;
			}
		}
		
		FileOutput fileHypothesis(strFileHypothesis,false);
		fileHypothesis.open();
		int iUtterance = 0;
		for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it, ++iUtterance) {
		
			const char *strUtteranceId = batchFile.getField(iUtterance,"id");
		
			cout << "processing utterance: " << strUtteranceId << endl;
			
			iFeatureVectorsTotal += it->mFeatures->getRows();
			Matrix<float> *mFeatures = it->mFeatures;
			
			if (hmmManager.getFeatureDim() != mFeatures->getCols()) {
				BVC_ERROR << "inconsistent feature dimensionality, HMMs: " << hmmManager.getFeatureDim() 
					<< ", features: " << mFeatures->getCols();
			}	
			
			decoder.beginUtterance();
			decoder.process(*mFeatures);
			
			// best path
			BestPath *bestPath = decoder.getBestPath();
			if (bestPath) {	
				// append the best path to a file (trn format)
				bestPath->write(fileHypothesis.getStream(),strUtteranceId);	
				bestPath->print(true);
				dLikelihoodTotal += bestPath->getPathScore();
			} else {
				cout << "no best path!!\n";
			}
			
			// hypothesis lattice
			if (bLatticeGeneration) {
				HypothesisLattice *hypothesisLattice = decoder.getHypothesisLattice();
				if (hypothesisLattice) {
					ostringstream ossText,ossBin;
					ossText << strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".txt";
					hypothesisLattice->store(ossText.str().c_str(),FILE_FORMAT_TEXT);
					ossBin << strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".bin";
					hypothesisLattice->store(ossBin.str().c_str(),FILE_FORMAT_BINARY);
					delete hypothesisLattice;
				} else {
					cout << "no hypothesis lattice!!\n";
				}
			}
			
			decoder.endUtterance();	
			
			// output features?
			if (bOutputFeatures) {
				ostringstream ossFileFeatures;
				ossFileFeatures << strFolderFeatures << PATH_SEPARATOR << strUtteranceId << ".fea"; 
				FeatureFile featureFile(ossFileFeatures.str().c_str(),MODE_WRITE);
				featureFile.store(*mFeatures);
			}
			
			// output alignment?
			if (bOutputAlignment && bestPath) {
			
				// create the state-level alignment and dump it to disk
				VPhoneAlignment *vPhoneAlignment = viterbi->align(*mFeatures,bestPath);
				if (vPhoneAlignment) {
					AlignmentFile alignmentFile(&phoneSet,&lexiconManager);
					ostringstream ossFileAlignment;
					ossFileAlignment << strFolderAlignments << PATH_SEPARATOR << strUtteranceId << ".ali";
					alignmentFile.store(*vPhoneAlignment,ossFileAlignment.str().c_str());
					AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
				} else {
					BVC_WARNING << "unable to perform the best-path alignment";
				}			
			}
			
			// clean-up
			if (bestPath) {
				delete bestPath;
			}	
			delete [] it->samples.sSamples;
			delete it->mFeatures;
		}
		fileHypothesis.close();
		
		double dTimeEnd = TimeUtils::getTimeMilliseconds();
		double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
		double dRTF = dTimeSeconds/(((float)iFeatureVectorsTotal)/100.0);
		
		BVC_INFORMATION << "- summary ------------------------------------";
		BVC_INFORMATION << "# utterances: " << iUtterance << " speech time: " << FLT(8,2) << 
			((float)iFeatureVectorsTotal)/100.0 << " seconds";
		BVC_INFORMATION << "decoding time: " << FLT(8,2) << dTimeSeconds << " seconds (RTF: " << 
			FLT(5,2) << dRTF << ")";
		BVC_INFORMATION << "likelihood: " << FLT(12,4) << dLikelihoodTotal << " (per frame: " << 
			FLT(8,4) << dLikelihoodTotal/((float)iFeatureVectorsTotal) << ")";
		BVC_INFORMATION << "----------------------------------------------";
		
		// uninitialize the decoder
		decoder.uninitialize();
		
		delete network;
		if (bOutputAlignment) {
			delete viterbi;
		}
	} 
	catch (std::runtime_error &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	
	return 0;
}
