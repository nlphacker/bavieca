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

void getHMMStateDecodingComposite(PhoneSet *m_phoneSet, HMMManager *m_hmmManager, VLexUnit &vLexUnitText, VHMMStateDecoding &vHMMStateDecodingComposite, LexUnit *lexUnitLeft, LexUnit *lexUnitRight);

// main for the tool "dynamicdecoder"
int main(int argc, char *argv[]) {

	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("dynamicdecoder",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);	
	m_commandLineManager->defineParameter("-cfg","configuration file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-hyp","hypothesis file",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-bat","batch file with entries [rawFile/featureFile utteranceId]",
		PARAMETER_TYPE_FILE,false);
	
	// (2) process command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get command line parameters
	const char *m_strFileConfiguration = m_commandLineManager->getParameterValue("-config");
	const char *m_strFileControl = m_commandLineManager->getParameterValue("-list");
	const char *m_strFileHypothesis = m_commandLineManager->getParameterValue("-hyp");

	PhoneSet *m_phoneSet = NULL;
	HMMManager *m_hmmManager = NULL;
	LexiconManager *m_lexiconManager = NULL;
	LMManager *m_lmManager = NULL;
	
	// load the configuration file
	ConfigurationDynamicDecoder *m_configuration = new ConfigurationDynamicDecoder(m_strFileConfiguration);
	m_configuration->load();
	
	// get parameters from the configuration file

	// feature extraction
	const char *m_strFileConfigurationFeatures = 
		m_configuration->getStrParameterValue("feature.configurationFile");
	const char *m_strCepstralNormalizationMode =
		m_configuration->getStrParameterValue("feature.cepstralNormalization.mode");
	const char *m_strCepstralNormalizationMethod =
		m_configuration->getStrParameterValue("feature.cepstralNormalization.method");
	int m_iCepstralNormalizationBufferSize = 
		m_configuration->getIntParameterValue("feature.cepstralNormalization.bufferSize");
	float m_fWarpFactor = 1.0;
	if (m_configuration->isParameterSet("feature.warpFactor")) {
		m_fWarpFactor = atof(m_configuration->getParameterValue("feature.warpFactor"));
	}
	const char *m_strFileFeatureTransform = NULL;
	if (m_configuration->isParameterSet("feature.transformFile")) {
		m_strFileFeatureTransform = m_configuration->getParameterValue("feature.transformFile");
	}
		
	// phone set
	const char *m_strFilePhoneSet = 
		m_configuration->getStrParameterValue("phoneticSymbolSet.file");
	
	// acoustic models
	const char *m_strFileAcousticModels = 
		m_configuration->getStrParameterValue("acousticModels.file"); 
	
	// language model
	const char *m_strLanguageModelFile = 
		m_configuration->getStrParameterValue("languageModel.file"); 
	const char *m_strLanguageModelFormat = 
		m_configuration->getStrParameterValue("languageModel.format"); 
	const char *m_strLanguageModelType = 
		m_configuration->getStrParameterValue("languageModel.type"); 
	float m_fLanguageModelScalingFactor = 
		m_configuration->getFloatParameterValue("languageModel.scalingFactor"); 
	const char *m_strLanguageModelNGram = 
		m_configuration->getStrParameterValue("languageModel.ngram"); 
	int m_iNGram = LMManager::getNGram(m_strLanguageModelNGram);
	//bool m_bLanguageCrossUtterance = 
	//	m_configuration->getBoolParameterValue("languageModel.crossUtterance");
		 
	// lexicon
	const char *m_strFileLexicon = 
		m_configuration->getStrParameterValue("lexicon.file"); 
		
	// insertion penalty
	float m_fInsertionPenaltyStandard = 
		m_configuration->getFloatParameterValue("insertionPenalty.standard"); 
	float m_fInsertionPenaltyFiller = 
		m_configuration->getFloatParameterValue("insertionPenalty.filler"); 
	const char *m_strFileInsertionPenaltyFiller = 
		m_configuration->getStrParameterValue("insertionPenalty.filler.file"); 
	
	// pruning parameters
	int m_iMaxActiveArcs = m_configuration->getIntParameterValue("pruning.maxActiveArcs");
	int m_iMaxActiveArcsWE = m_configuration->getIntParameterValue("pruning.maxActiveArcsWE");
	int m_iMaxActiveTokensArc = m_configuration->getIntParameterValue("pruning.maxActiveTokensArc");	
	float m_fBeamWidthArcs = m_configuration->getFloatParameterValue("pruning.likelihoodBeam");
	float m_fBeamWidthArcsWE = m_configuration->getFloatParameterValue("pruning.likelihoodBeamWE");
	float m_fBeamWidthTokensArc = m_configuration->getFloatParameterValue("pruning.likelihoodBeamTokensArc");
	
	// output lattice?
	bool m_bLatticeGeneration = m_configuration->isParameterSet("output.lattice.folder");	
	const char *m_strFolderLattices = NULL;
	int m_iMaxWordSequencesState = -1;
	if (m_bLatticeGeneration) {
		m_strFolderLattices = m_configuration->getStrParameterValue("output.lattice.folder"); 
		m_iMaxWordSequencesState = 
			m_configuration->getIntParameterValue("output.lattice.maxWordSequencesState");
	}
	
	// output features?
	bool m_bOutputFeatures = m_configuration->isParameterSet("output.features.folder");
	const char *m_strFolderFeatures = NULL;
	if (m_bOutputFeatures) {
		m_strFolderFeatures = m_configuration->getStrParameterValue("output.features.folder"); ;
	}
	
	// output alignment?
	bool m_bOutputAlignment = m_configuration->isParameterSet("output.alignment.folder");
	const char *m_strFolderAlignments = NULL;
	if (m_bOutputAlignment) {
		m_strFolderAlignments = m_configuration->getStrParameterValue("output.alignment.folder"); ;
	}

	// output audio?
	bool m_bOutputAudio = m_configuration->isParameterSet("output.audio.folder");
	const char *m_strFolderAudio = NULL;
	if (m_bOutputAudio) {
		m_strFolderAudio = m_configuration->getStrParameterValue("output.audio.folder");
	}

   // load the phone set
   m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();

   // load the lexicon
   m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
   m_lexiconManager->load();
   m_lexiconManager->attachLexUnitPenalties(m_fInsertionPenaltyStandard,m_fInsertionPenaltyFiller);
	
	FillerManager m_fillerManager(m_strFileInsertionPenaltyFiller);	
	m_fillerManager.load();
	m_fillerManager.attachInsertionPenaltyFillers(m_lexiconManager);
	
   // load the feature configuration
   ConfigurationFeatures *configurationFeatures = new ConfigurationFeatures(m_strFileConfigurationFeatures);
   configurationFeatures->load();
   
   // get the feature normalization mode and method
   int m_iCepstralNormalizationMode = FeatureExtractor::getNormalizationMode(m_strCepstralNormalizationMode);
   int m_iCepstralNormalizationMethod = FeatureExtractor::getNormalizationMethod(m_strCepstralNormalizationMethod);
   
   // create the feature extractor
   FeatureExtractor *m_featureExtractor = new FeatureExtractor(configurationFeatures,m_fWarpFactor,
   	m_iCepstralNormalizationBufferSize,m_iCepstralNormalizationMode,m_iCepstralNormalizationMethod);
   m_featureExtractor->initialize();
   
   delete configurationFeatures;	
   
   // load the feature transforms
   VTransform m_vTransformFeatures;
   if (m_strFileFeatureTransform) {
		BatchFile *m_batchFile = new BatchFile(m_strFileFeatureTransform,"transform");
		m_batchFile->load();
		for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {
			Transform *transform = new Transform();
			transform->load(m_batchFile->getField(i,0u));
			m_vTransformFeatures.push_back(transform);
		}	
		delete m_batchFile;	
   }
		
	// load the HMMs used for the estimation
	m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileAcousticModels);
	m_hmmManager->initializeDecoding();
	
	// create the aligner object?
	Viterbi *m_viterbi = NULL;
	if (m_bOutputAlignment) {
		m_viterbi = new Viterbi(m_phoneSet,m_hmmManager,m_lexiconManager);	
	}
	
   // load the language model
	m_lmManager = new LMManager(m_lexiconManager,
										m_strLanguageModelFile,
										m_strLanguageModelFormat,
										m_strLanguageModelType,
										m_strLanguageModelNGram); 
	m_lmManager->load();
	m_lmManager->buildLMGraph();
	
	NetworkBuilderX *m_networkBuilder = new NetworkBuilderX(m_phoneSet,m_hmmManager,m_lexiconManager);
	
	// build the decoding network
	DynamicNetworkX *m_network = m_networkBuilder->build();
	if (m_network == NULL) {
		BVC_ERROR << "unable to build the network";
	}

	DynamicDecoderX *m_decoder = new DynamicDecoderX(m_phoneSet,m_hmmManager,m_lexiconManager, 
			m_lmManager,m_fLanguageModelScalingFactor,m_iNGram,m_network,m_iMaxActiveArcs,
			m_iMaxActiveArcsWE,m_iMaxActiveTokensArc,m_fBeamWidthArcs,m_fBeamWidthArcsWE,m_fBeamWidthTokensArc,
			m_bLatticeGeneration,m_iMaxWordSequencesState);

	// initialize the decoder
	m_decoder->initialize();
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();	
	
	double dLikelihoodTotal = 0.0;
	int iFeatureVectorsTotal = 0;
	int iUtterances = 0;
	
	// load the batch file
	BatchFile batchFile(m_strFileControl,"audio|id");
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
		utteranceData.features.fFeatures = NULL;
		utteranceData.features.iFeatures = -1;
		vUtteranceData.push_back(utteranceData);
	}	
	
	// extract features
	m_featureExtractor->extractFeaturesSession(vUtteranceData,true);
	
	// apply feture transforms
	int iDimFea = m_featureExtractor->getFeatureDimensionality();
	for(VTransform::iterator it = m_vTransformFeatures.begin() ; it != m_vTransformFeatures.end() ; ++it) {
		printf("%d -> %d\n",iDimFea,(*it)->getRows());
		for(VUtteranceData::iterator jt = vUtteranceData.begin() ; jt != vUtteranceData.end() ; ++jt) {
			float *fFeaturesX = new float[jt->features.iFeatures*(*it)->getRows()];
			for(int i=0 ; i < jt->features.iFeatures ; ++i) {
				assert(0);
				//(*it)->apply(jt->features.fFeatures+(i*iDimFea),fFeaturesX+(i*(*it)->getRows()));
				// this needs to be fixed
				assert(0);
			}
			delete [] jt->features.fFeatures;
			jt->features.fFeatures = fFeaturesX;
		}
		iDimFea = (*it)->getRows();
	}	
	
	FileOutput fileHypothesis(m_strFileHypothesis,false);
	fileHypothesis.open();
	int iUtterance = 0;
	for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it, ++iUtterance) {
	
		const char *strUtteranceId = batchFile.getField(iUtterance,"id");
	
		cout << "processing utterance: " << strUtteranceId << endl;
		
		iFeatureVectorsTotal += it->features.iFeatures;
		float *fFeatureVectors = it->features.fFeatures;
		int iFeatureVectors = it->features.iFeatures;
		
		m_decoder->beginUtterance();
		m_decoder->process(fFeatureVectors,iFeatureVectors);
		
		// best path
		BestPath *bestPath = m_decoder->getBestPath();
		if (bestPath != NULL) {	
			// append the best path to a file (trn format)
			bestPath->write(fileHypothesis.getStream(),strUtteranceId);	
			bestPath->print(true);
			dLikelihoodTotal += bestPath->getPathScore();
		} else {
			cout << "no best path!!\n";
		}
		
		// hypothesis lattice
		if (m_bLatticeGeneration) {
			HypothesisLattice *hypothesisLattice = m_decoder->getHypothesisLattice();
			if (hypothesisLattice) {
				ostringstream ossText,ossBin;
				ossText << m_strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".txt";
				hypothesisLattice->store(ossText.str().c_str(),FILE_FORMAT_TEXT);
				ossBin << m_strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".bin";
				hypothesisLattice->store(ossBin.str().c_str(),FILE_FORMAT_BINARY);
				delete hypothesisLattice;
			} else {
				cout << "no hypothesis lattice!!\n";
			}
		}
		
		m_decoder->endUtterance();	
		
		// output features?
		if (m_bOutputFeatures) {
			ostringstream ossFileFeatures;
			ossFileFeatures << m_strFolderFeatures << PATH_SEPARATOR << strUtteranceId << ".fea"; 
			FeatureFile featureFile(ossFileFeatures.str().c_str(),MODE_WRITE);
			featureFile.store(fFeatureVectors,iFeatureVectors);
		}
		
		// output alignment?
		if (m_bOutputAlignment && bestPath) {
		
			// create the state-level alignment and dump it to disk
			VPhoneAlignment *vPhoneAlignment = m_viterbi->align(fFeatureVectors,iFeatureVectors,bestPath);
			if (vPhoneAlignment) {
				AlignmentFile alignmentFile(m_phoneSet,m_lexiconManager);
				ostringstream ossFileAlignment;
				ossFileAlignment << m_strFolderAlignments << PATH_SEPARATOR << strUtteranceId << ".ali";
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
		delete [] it->features.fFeatures;
	}
	fileHypothesis.close();
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	double dRTF = dTimeSeconds/(((float)iFeatureVectorsTotal)/100.0);
	
	printf("- summary ------------------------------------\n");
	printf("# utterances: %d speech time: %.2f seconds\n",iUtterances,((float)iFeatureVectorsTotal)/100.0);
	printf("decoding time: %.2f seconds (RTF: %5.2f)\n",dTimeSeconds,dRTF);
	printf("likelihood: %.4f (per frame: %8.4f)\n",dLikelihoodTotal,dLikelihoodTotal/((float)iFeatureVectorsTotal));
	printf("----------------------------------------------\n");
	
	// uninitialize the decoder
	m_decoder->uninitialize();
	
	delete m_decoder;
	delete m_network;
	delete m_networkBuilder;
	delete m_lmManager;
	delete m_lexiconManager;
	delete m_featureExtractor;
	delete m_hmmManager;
	delete m_phoneSet;
	delete m_configuration;
	delete m_commandLineManager;
	if (m_bOutputAlignment) {
		delete m_viterbi;
	}
	
	return 0;
}
