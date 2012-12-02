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


#include "BaviecaAPI.h"

#include "ViterbiX.h"
#include "AlignmentFile.h"
#include "AudioFile.h"
#include "BatchFile.h"
#include "BestPath.h"
#include "ConfigurationBavieca.h"
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
#include "LogMessage.h"
#include "PhoneSet.h"
#include "SADModule.h"
#include "TextAligner.h"
#include "TextAlignment.h"

namespace Bavieca {

// constructor (receives default configuration parameters)
BaviecaAPI::BaviecaAPI(const char *strFileConfiguration) {

	assert(strFileConfiguration);
	m_strFileConfiguration = strFileConfiguration;
	
	m_configuration = NULL;
	m_featureExtractor = NULL;
	m_phoneSet = NULL;
	m_lexiconManager = NULL;
	m_sadModule = NULL;
	m_lmManager = NULL;
	m_viterbiX = NULL;	
	m_hmmManager = NULL;
	m_network = NULL;
	m_networkBuilder = NULL;
	m_dynamicDecoder = NULL;
	m_bInitialized = false;
}

// destructor
BaviecaAPI::~BaviecaAPI() {
}

// initialize API (overriding parameters as needed)
bool BaviecaAPI::initialize(unsigned char iFlags, ParamValueI *paramValue, int iParameters) {

	assert(m_bInitialized == false);

	try {
	
		m_iFlags = iFlags;
	
		// (1) load configuration parameters
		m_configuration = new ConfigurationBavieca(m_strFileConfiguration.c_str());
		m_configuration->load();
		
		// (2) override parameters if needed
		if (iParameters > 0) {	
			assert(paramValue);
			for(int i=0 ; i < iParameters ; ++i) {
				assert(paramValue[i].strParameter && paramValue[i].strValue);
				m_configuration->setParameterValue(paramValue[i].strParameter,paramValue[i].strValue);
			}
		}
		
		// (3) get configuration parameters
		
		// phone set
		const char *m_strFilePhoneSet = 
			m_configuration->getStrParameterValue("phoneticSymbolSet.file");
		
		// acoustic models
		const char *m_strFileAcousticModels = 
			m_configuration->getStrParameterValue("acousticModels.file"); 
		
		// features
		const char *m_strCepstralNormalizationMode =
			m_configuration->getStrParameterValue("feature.cepstralNormalization.mode");
		const char *m_strCepstralNormalizationMethod =
			m_configuration->getStrParameterValue("feature.cepstralNormalization.method");
		int m_iCepstralNormalizationBufferSize = 
			m_configuration->getIntParameterValue("feature.cepstralNormalization.bufferSize");
		float m_fWarpFactor = 1.0;
		if (m_configuration->isParameterSet("feature.warpFactor")) {
			m_fWarpFactor = (float)atof(m_configuration->getParameterValue("feature.warpFactor"));
		}
		
		// load the phone set
		m_phoneSet = new PhoneSet(m_strFilePhoneSet); 
		m_phoneSet->load();
	
		// load the HMMs used for the estimation
		m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
		m_hmmManager->load(m_strFileAcousticModels);
		m_hmmManager->initializeDecoding();
		
		// load the feature configuration
		const char *m_strFileConfigurationFeatures = m_configuration->getStrParameterValue("feature.configurationFile");
		ConfigurationFeatures *configurationFeatures = new ConfigurationFeatures(m_strFileConfigurationFeatures);
		configurationFeatures->load();	
		
		// get the feature normalization mode and method
		int m_iCepstralNormalizationMode = FeatureExtractor::getNormalizationMode(m_strCepstralNormalizationMode);
		int m_iCepstralNormalizationMethod = FeatureExtractor::getNormalizationMethod(m_strCepstralNormalizationMethod);
		
		// create the feature extractor
		m_featureExtractor = new FeatureExtractor(configurationFeatures,m_fWarpFactor,
			m_iCepstralNormalizationBufferSize,m_iCepstralNormalizationMode,m_iCepstralNormalizationMethod);
		m_featureExtractor->initialize();
		delete configurationFeatures;
		
		// speech activity detection
		if (m_iFlags & INIT_SAD) {
		
			// make sure required configuration parameters are defined
			if (!m_configuration->areSADParametersSet()) {
				BVC_ERROR << "wrong configuration parameters: not all the required parameters are set";
			}	
		
			// speech activity detection
			int m_iMaxGaussianComponentsSilence = m_configuration->getIntParameterValue("sad.maxGaussianSilence");
			int m_iMaxGaussianComponentsSpeech = m_configuration->getIntParameterValue("sad.maxGaussianSpeech");
			float m_fPenaltySilenceToSpeech = m_configuration->getFloatParameterValue("sad.speechPenalty");
			int m_iFramesPadding = m_configuration->getIntParameterValue("sad.speechPadding");
		
			m_sadModule = new SADModule(m_phoneSet,m_hmmManager,m_iMaxGaussianComponentsSilence,
				m_iMaxGaussianComponentsSpeech,m_fPenaltySilenceToSpeech,m_iFramesPadding);
			m_sadModule->initialize();
		}
		
		// aligner
		if (m_iFlags & INIT_ALIGNER) {
		
			// make sure required configuration parameters are defined
			if (!m_configuration->areAlignerParametersSet()) {
				BVC_ERROR << "wrong configuration parameters: not all the required parameters are set";
			}
			
			// lexicon
			const char *m_strFileLexicon = 
				m_configuration->getStrParameterValue("lexicon.file");	
		
			// load the lexicon
			m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
			m_lexiconManager->load();
			
			float fBeamWidth = 1000;
			m_viterbiX = new ViterbiX(m_phoneSet,m_lexiconManager,m_hmmManager,fBeamWidth,2000,true,2000);	
		}
		
		// decoder
		if (m_iFlags & INIT_DECODER) {
		
			// make sure required configuration parameters are defined
			if (!m_configuration->areDecodingParametersSet()) {
				BVC_ERROR << "wrong configuration parameters: not all the required parameters are set";
			}
			
			// lexicon
			if (!m_lexiconManager) {
				const char *m_strFileLexicon = 
					m_configuration->getStrParameterValue("lexicon.file");	
		
				// load the lexicon
				m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet); 
				m_lexiconManager->load();
			}
			
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
				//m_configuration->getBoolParameterValue("languageModel.crossUtterance");
					
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
			bool m_bLatticeGeneration = m_configuration->isParameterSet("output.lattice.maxWordSequencesState");	
			//m_bLatticeGeneration = false;
			int m_iMaxWordSequencesState = -1;
			if (m_bLatticeGeneration) {
				m_iMaxWordSequencesState = 
					m_configuration->getIntParameterValue("output.lattice.maxWordSequencesState");
			}
			
			// attach insertion penalties to lexical units
			m_lexiconManager->attachLexUnitPenalties(m_fInsertionPenaltyStandard,m_fInsertionPenaltyFiller);
			
			FillerManager m_fillerManager(m_strFileInsertionPenaltyFiller);	
			m_fillerManager.load();
			m_fillerManager.attachInsertionPenaltyFillers(m_lexiconManager);	
			
			// load the language model
			m_lmManager = new LMManager(m_lexiconManager,
												m_strLanguageModelFile,
												m_strLanguageModelFormat,
												m_strLanguageModelType,
												m_strLanguageModelNGram); 
			m_lmManager->load();
			m_lmManager->buildLMGraph();
			
			m_networkBuilder = new NetworkBuilderX(m_phoneSet,m_hmmManager,m_lexiconManager);
			
			// build the decoding network
			m_network = m_networkBuilder->build();
			if (m_network == NULL) {
				BVC_ERROR << "unable to build the network";
			}
		
			m_dynamicDecoder = new DynamicDecoderX(m_phoneSet,m_hmmManager,m_lexiconManager,
					m_lmManager,m_fLanguageModelScalingFactor,m_iNGram,m_network,m_iMaxActiveArcs,
					m_iMaxActiveArcsWE,m_iMaxActiveTokensArc,m_fBeamWidthArcs,m_fBeamWidthArcsWE,m_fBeamWidthTokensArc,
					m_bLatticeGeneration,m_iMaxWordSequencesState);
					
			m_dynamicDecoder->initialize();
		}
		
		// speaker adaptation
		if (m_iFlags & INIT_ADAPTATION) {
			
		}
			
	} catch (ExceptionBase) {	
		return false;
	}	

	m_bInitialized = true;
	return true;
}

// uninitialize the API	
void BaviecaAPI::uninitialize() {

	assert(m_bInitialized);
	
	delete m_configuration;
	delete m_phoneSet;
	delete m_lexiconManager;
	delete m_hmmManager;
	delete m_featureExtractor;
	
	if (m_sadModule) {
		delete m_sadModule;
	}
	if (m_viterbiX) {
		delete m_viterbiX;
	}
	if (m_dynamicDecoder) {
		m_dynamicDecoder->uninitialize();
		delete m_dynamicDecoder;
		delete m_network;
		delete m_networkBuilder;
	}
	if (m_lmManager) {
		delete m_lmManager;
	}
	
	m_configuration = NULL;
	m_featureExtractor = NULL;
	m_phoneSet = NULL;
	m_lexiconManager = NULL;
	m_sadModule = NULL;
	m_lmManager = NULL;
	m_viterbiX = NULL;	
	m_hmmManager = NULL;
	m_network = NULL;
	m_networkBuilder = NULL;
	m_dynamicDecoder = NULL;
	m_bInitialized = false;	
}

// extract features from the audio
float *BaviecaAPI::extractFeatures(short *sSamples, int iSamples, int *iFeatures) {

	assert(m_bInitialized);
	
	return m_featureExtractor->extractFeaturesStream(sSamples,iSamples,iFeatures);
}

// return feature dimensionality
int BaviecaAPI::getFeatureDim() {

	assert(m_bInitialized);
	assert(m_featureExtractor);	
	return m_featureExtractor->getFeatureDimensionality();
}

// free features extracted using extractFeatures(...)
void BaviecaAPI::free(float *fFeatures) {
	delete [] fFeatures;	
}

// start a SAD session
void BaviecaAPI::sadBeginSession() {
	
	assert(m_bInitialized);
	m_featureExtractor->resetStream();
	m_sadModule->beginSession();
}

// terminate a SAD session
void BaviecaAPI::sadEndSession() {

	assert(m_bInitialized);
	m_sadModule->endSession();
}

// proces the given features
void BaviecaAPI::sadFeed(float *fFeatures, int iFeatures) {

	assert(m_bInitialized);
	m_sadModule->processFeatures(fFeatures,iFeatures);
}

// recover speech segments by doing back-tracking on the grid
SpeechSegmentI *BaviecaAPI::sadRecoverSpeechSegments(int *iSegments) {

	//m_sadModule->printGrid();
	assert(m_bInitialized);
	SpeechSegmentI *speechSegmentI;
	VSpeechSegment vSpeechSegment;
	m_sadModule->recoverSpeechSegments(vSpeechSegment);
	if (vSpeechSegment.empty()) {
		*iSegments = 0;
		return NULL;
	}	
	speechSegmentI = new SpeechSegmentI[vSpeechSegment.size()];
	int i = 0;
	for(VSpeechSegment::iterator it = vSpeechSegment.begin() ; it != vSpeechSegment.end() ; ++it, ++i) {
		speechSegmentI[i].iFrameStart = (*it)->iFrameStart;
		speechSegmentI[i].iFrameEnd = (*it)->iFrameEnd;
		delete *it;
	}
	*iSegments = vSpeechSegment.size();
	return speechSegmentI;
}

// free speech segments returned by sarRecoverSpeechSegments(...)
void BaviecaAPI::free(SpeechSegmentI *speechSegments) {
	delete [] speechSegments;
}


// forced alignment between features and audio
WordAlignmentI *BaviecaAPI::align(float *fFeatures, int iFeatures, const char *strText, 
	bool bMultiplePronunciations, int *iWords) {

	assert(m_bInitialized);
	VLexUnit vLexUnit;
	bool bAllKnown;
	if (m_lexiconManager->getLexUnits(strText,vLexUnit,bAllKnown) == false) {
		return NULL;
	}
	if (bAllKnown == false) {
		return NULL;
	}	
	// create the aligner
	double dUtteranceLikelihood;
	int iErrorCode;
	
	VLexUnit vLexUnitOptional;
	vLexUnitOptional.push_back(m_lexiconManager->getLexUnitSilence());
	Alignment *alignment = m_viterbiX->processUtterance(vLexUnit,bMultiplePronunciations,
			vLexUnitOptional,fFeatures,iFeatures,&dUtteranceLikelihood,iErrorCode);
	if (alignment == NULL) {
		return NULL;
	}
	VPhoneAlignment *vPhoneAlignment = alignment->getPhoneAlignment(m_lexiconManager);
	if (vPhoneAlignment == NULL) {
		return NULL;
	}
	VLexUnitAlignment *vLexUnitAlignment = AlignmentFile::getVLexUnitAlignment(*vPhoneAlignment);
	assert(vLexUnitAlignment);
	assert(vLexUnitAlignment->size() >= vLexUnit.size());
	WordAlignmentI *wordAlignmentI = new WordAlignmentI[vLexUnitAlignment->size()];
	int iPhone = 0;
	int i=0;
	for(VLexUnitAlignment::iterator it = vLexUnitAlignment->begin() ; it != vLexUnitAlignment->end() ; ++it, ++i) {
		const char *strLexUnit = m_lexiconManager->getStrLexUnitPron((*it)->lexUnit->iLexUnitPron);
		wordAlignmentI[i].strWord = new char[strlen(strLexUnit)+1];
		strcpy(wordAlignmentI[i].strWord,strLexUnit);
		wordAlignmentI[i].iFrameStart = (*it)->iBegin;
		wordAlignmentI[i].iFrameEnd = (*it)->iEnd;
		// phone-alignment
		wordAlignmentI[i].iPhones = (*it)->lexUnit->vPhones.size();
		wordAlignmentI[i].phoneAlignment = new PhoneAlignmentI[wordAlignmentI[i].iPhones];
		for(int j=0 ; j<wordAlignmentI[i].iPhones; ++j) {
			assert(iPhone < (int)vPhoneAlignment->size());
			const char *strPhone = m_phoneSet->getStrPhone((*vPhoneAlignment)[iPhone]->iPhone);
			assert(strPhone);
			wordAlignmentI[i].phoneAlignment[j].strPhone = new char[strlen(strPhone)+1];
			strcpy(wordAlignmentI[i].phoneAlignment[j].strPhone,strPhone);
			wordAlignmentI[i].phoneAlignment[j].iFrameStart = (*vPhoneAlignment)[iPhone]->iStateBegin[0];
			wordAlignmentI[i].phoneAlignment[j].iFrameEnd = (*vPhoneAlignment)[iPhone]->iStateEnd[NUMBER_HMM_STATES-1];
			++iPhone;
		}
	}	
	
	*iWords = vLexUnitAlignment->size();
	AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
	AlignmentFile::destroyLexUnitAlignment(vLexUnitAlignment);
	delete alignment;
	
	return wordAlignmentI;
}

// free word alignments returned by align(...)
void BaviecaAPI::free(WordAlignmentI *wordAlignments, int iWords) {
	for(int i=0 ; i < iWords ; ++i) {
		delete [] wordAlignments[i].strWord;
		for(int j=0 ; j < wordAlignments[i].iPhones ; ++j) {
			delete [] wordAlignments[i].phoneAlignment[j].strPhone;
		}
		delete [] wordAlignments[i].phoneAlignment;
	}
	delete [] wordAlignments;
}

// DECODING ------------------------------------------------------------------------------------------

// signal beginning of utterance
void BaviecaAPI::decBeginUtterance() {

	assert(m_bInitialized);
	assert(m_iFlags & INIT_DECODER);
	m_dynamicDecoder->beginUtterance();
}

// process feature vectors from an utterance
void BaviecaAPI::decProcess(float *fFeatures, int iFeatures) {

	assert(m_bInitialized);
	assert(m_iFlags & INIT_DECODER);
	m_dynamicDecoder->process(fFeatures,iFeatures);	
}

// get decoding results
WordHypothesisI *BaviecaAPI::decGetHypothesis(int *iWords, const char *strFileHypothesisLattice) {

	assert(m_bInitialized);
	assert(m_iFlags & INIT_DECODER);
	
	// best path
	WordHypothesisI *wordHypothesisI = NULL;
	BestPath *bestPath = m_dynamicDecoder->getBestPath();
	if (bestPath != NULL) {
		VLexUnit vLexUnit;
		LBestPathElement *lBestPathElements = bestPath->getBestPathElements();
		// (1) count standard lexical units
		int iStandard = 0;
		for(LBestPathElement::iterator it = lBestPathElements->begin() ; it != lBestPathElements->end() ; ++it) {
			if (m_lexiconManager->isStandard((*it)->lexUnit)) {
				++iStandard;
			}
		}
		*iWords = iStandard;
		// (2) extract standard lexical units with alignment information
		wordHypothesisI = new WordHypothesisI[iStandard];
		int i=0;
		for(LBestPathElement::iterator it = lBestPathElements->begin() ; it != lBestPathElements->end() ; ++it) {
			if (m_lexiconManager->isStandard((*it)->lexUnit)) {
				const char *strLexUnit = m_lexiconManager->getStrLexUnitPron((*it)->lexUnit->iLexUnitPron);
				wordHypothesisI[i].strWord = new char[strlen(strLexUnit)+1];
				strcpy(wordHypothesisI[i].strWord,strLexUnit);
				wordHypothesisI[i].iFrameStart = (*it)->iFrameStart;
				wordHypothesisI[i].iFrameEnd = (*it)->iFrameEnd;	
				++i;
			}
		}
		delete bestPath;
	}
	
	// hypothesis lattice
	if (strFileHypothesisLattice) {
		assert(m_bLatticeGeneration);
		HypothesisLattice *hypothesisLattice = m_dynamicDecoder->getHypothesisLattice();
		if (hypothesisLattice != NULL) {
			//hypothesisLattice->store(strFileHypothesisLattice,FILE_FORMAT_TEXT);
			hypothesisLattice->store(strFileHypothesisLattice,FILE_FORMAT_BINARY);
			delete hypothesisLattice;
		}
	}

	return wordHypothesisI;
}

// signal end of utterance
void BaviecaAPI::decEndUtterance() {

	assert(m_bInitialized);
	assert(m_iFlags & INIT_DECODER);
	m_dynamicDecoder->endUtterance();
}

// free word hypotheses returned by getHypothesis(...)
void BaviecaAPI::free(WordHypothesisI *wordHypothesis, int iWords) {
	for(int i=0 ; i < iWords ; ++i) {
		delete [] wordHypothesis[i].strWord;
	}
	delete [] wordHypothesis;
}

// return a word-level assessment given a hypothesis and a reference text
TAElementI *BaviecaAPI::getAssessment(WordHypothesisI *wordHypothesis, int iWordsHyp, const char *strReference, int *iTAElements) {

	// get lexical units from reference (there might be out-of-vocabulary words <unk>)
	VLexUnit vLexUnitRef;
	bool bAllKnown;
	if (m_lexiconManager->getLexUnits(strReference,vLexUnitRef,bAllKnown) == false) {
		return NULL;
	}
	
	// get lexical units from the hypothesis
	VLexUnit vLexUnitHyp;
	for(int i=0 ; i < iWordsHyp ; ++i) {
		LexUnit *lexUnit = m_lexiconManager->getLexUnitPronunciation(wordHypothesis[i].strWord);
		assert(lexUnit);
		vLexUnitHyp.push_back(lexUnit);
	}	

	// perform the actual alignment
	TextAligner textAligner(m_lexiconManager);
	TextAlignment *textAlignment = textAligner.align(vLexUnitHyp,vLexUnitRef);
	VTAElement &vTAElement = textAlignment->getAlignment();
	TAElementI *taElementI = new TAElementI[vTAElement.size()]; 
	int i=0;
	for(VTAElement::iterator it = vTAElement.begin() ; it != vTAElement.end() ; ++it,++i) {
		if ((*it)->iAlignmentEvent == TEXT_ALIGNMENT_EVENT_CORRECT) {
			taElementI[i].iEvent = TAE_CORRECT;
		} else if ((*it)->iAlignmentEvent == TEXT_ALIGNMENT_EVENT_SUBSTITUTION) {		
			taElementI[i].iEvent = TAE_SUBSTITUTION;
		} else if ((*it)->iAlignmentEvent == TEXT_ALIGNMENT_EVENT_DELETION) {
			taElementI[i].iEvent = TAE_DELETION;
		} else if ((*it)->iAlignmentEvent == TEXT_ALIGNMENT_EVENT_INSERTION) {
			taElementI[i].iEvent = TAE_INSERTION;
		}
		taElementI[i].iIndexRef = (*it)->iIndexReference;
		if ((*it)->iIndexReference != -1) {
			const char *strLexUnit = m_lexiconManager->getStrLexUnit((*it)->iLexUnitReference);
			taElementI[i].strWordRef = new char[strlen(strLexUnit)+1];
			strcpy(taElementI[i].strWordRef,strLexUnit);
		} else {
			taElementI[i].strWordRef = NULL;
		}
		taElementI[i].iIndexHyp = (*it)->iIndexHypothesis;
		if ((*it)->iIndexHypothesis != -1) {
			const char *strLexUnit = m_lexiconManager->getStrLexUnit((*it)->iLexUnitHypothesis);
			taElementI[i].strWordHyp = new char[strlen(strLexUnit)+1];
			strcpy(taElementI[i].strWordHyp,strLexUnit);
		} else {
			taElementI[i].strWordHyp = NULL;
		}
	}	
	*iTAElements = vTAElement.size();
	delete textAlignment;
	
	return taElementI;
}

// free text alignment elements returned by getAssessment(...)
void BaviecaAPI::free(TAElementI *taElements, int iElements) {
	for(int i=0 ; i < iElements ; ++i) {
		if (taElements[i].strWordRef) {
			delete [] taElements[i].strWordRef;
		}
		if (taElements[i].strWordHyp) {
			delete [] taElements[i].strWordHyp;
		}
	}
	assert(taElements);
	delete [] taElements;
}

// feed data into speaker adaptation
void BaviecaAPI::mllrFeed(const char *strReference, float *fFeatures, int iFeatures) {


}

// adapt current acoustic models using fed adaptation data
void BaviecaAPI::mllrAdapt() {


}

};	// end-of-namespace









