/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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


#include "ConfigurationDecoder.h"

// static variables definition and initilization

// (1) private members

const char *ConfigurationDecoder::m_strParametersAll[] =	{
	"system.task",
	"system.mode", 
	"decodingNetwork.compact.forward",
	"decodingNetwork.compact.backward",
	"input.type",	
	"audio.format",
	"audio.samplingRate",
	"audio.sampleSize",
	"speechActivityDetection.autoEndPoint",
	"speechActivityDetection.endPointPadding",
	"speechActivityDetection.penalty.silenceToSpeech",
	"speechActivityDetection.penalty.speechToSilence",
	"speechActivityDetection.penalty.silence",
	"speechActivityDetection.minSpeechFrames",
	"speechActivityDetection.minSilenceFrames",
	"speechActivityDetection.speechMixtures",
	"feature.configurationFile",
	"feature.cepstralNormalization",
	"feature.cepstralNormalization.method",
	"feature.cepstralNormalization.bufferSize",
	"feature.transform.file",
	"phoneticSymbolSet.file",
	"acousticModels.file",
	"acousticModels.format",
	"acousticModels.contextDependent",
	"durationModeling",
	"durationModeling.scalingFactor",
	"lexicon.file",
	"lexicon.alternativePronunciations",
	"languageModel.file",
	"languageModel.format",
	"languageModel.type",
	"languageModel.lookAhead",
	"languageModel.lookAhead.order",
	"languageModel.lookAhead.maxCacheSize",
	"languageModel.scalingFactor",
	"languageModel.ngram",
	"languageModel.crossUtterance",
	"pruning.global.maxTokens",
	"pruning.global.likelihood",
	"pruning.wordEnd.maxTokens",
	"pruning.wordEnd.likelihood",
	"pruning.wordIdentity.maxTokens",
	"pruning.wordIdentity.likelihood",
	"pruning.fanOut.maxTokens",
	"pruning.fanOut.likelihood",
	"log.file",
	"log.verbosityLevel",
	"insertionPenalty.standard",
	"insertionPenalty.filler",
	"insertionPenalty.filler.file",
	"adaptation.vtln",													// whether online VTLN will be applied
	"adaptation.vtln.estimation",										// VTL estimation mode (online/offline)
	"adaptation.vtln.warpFactor",										// warp factor to be applied to all the utterances (batch/online)
	"adaptation.vtln.warpFactor.file",								// file containing warp factors for each utterance ID (batch)
	"adaptation.mllr",													// whether online MLLR will be applied
	"adaptation.mllr.adaptCovariance",								// whether to adapt the covariance in addition to the mean
	"adaptation.mllr.regressionTree.file",							// file containing the regression tree
	"adaptation.mllr.regressionTree.minOccupation",				// minimum occupation of a tree-node to compute a transform
	"output.bestSinglePath",
	"output.bestSinglePath.outputSentenceDelimiters",
	"output.bestSinglePath.outputFillers",
	"output.bestSinglePath.outputConfidenceValues",
	"output.nBest",
	"output.nBest.n",
	"output.nBest.alignments.unique",
	"output.nBest.discardFillers",
	"output.nBest.singlePronunciation",
	"output.nBest.output.format",
	"output.nBest.output.alignments",
	"output.nBest.output.likelihood",
	"output.nBest.output.folder",
	"output.graph",
	"output.graph.pruning.forwardBackward",
	"output.graph.output.folder",
	"output.audio",
	"output.audio.output.folder",
	"output.features",
	"output.features.output.folder",
	"output.alignment",
	"output.alignment.output.folder",
	"confidenceEstimation.bestPath",
	"confidenceEstimation.method",
	"confidenceEstimation.graphPosteriors.insertionPenalty.standard",
	"confidenceEstimation.graphPosteriors.insertionPenalty.filler",
	"confidenceEstimation.graphPosteriors.scalingFactor.languageModel",
	"confidenceEstimation.graphPosteriors.confidenceMeasure",	
	};

int ConfigurationDecoder::m_iParametersAll = 87;


const char *ConfigurationDecoder::m_strParametersRequired[] = { 
	"system.task",
	"system.mode", 
	"decodingNetwork.compact.forward",
	"decodingNetwork.compact.backward",
	"input.type",		
	"audio.format",
	"audio.samplingRate",
	"audio.sampleSize",
	"speechActivityDetection.autoEndPoint",
	"speechActivityDetection.endPointPadding",
	"speechActivityDetection.penalty.silenceToSpeech",
	"speechActivityDetection.penalty.speechToSilence",
	"speechActivityDetection.penalty.silence",
	"speechActivityDetection.minSpeechFrames",
	"speechActivityDetection.minSilenceFrames",
	"speechActivityDetection.speechMixtures",
	"feature.configurationFile",	
	"feature.cepstralNormalization",
	"feature.cepstralNormalization.method",
	"feature.cepstralNormalization.bufferSize",
	"phoneticSymbolSet.file",
	"acousticModels.file",
	"acousticModels.format",
	"acousticModels.contextDependent",
	"durationModeling",
	"durationModeling.scalingFactor",	
	"lexicon.file",
	"lexicon.alternativePronunciations",
	"languageModel.file",
	"languageModel.format",
	"languageModel.type",
	"languageModel.lookAhead",
	"languageModel.lookAhead.order",	
	"languageModel.lookAhead.maxCacheSize",	
	"languageModel.scalingFactor",
	"languageModel.ngram",
	"languageModel.crossUtterance",
	"pruning.global.maxTokens",
	"pruning.global.likelihood",
	"pruning.wordEnd.maxTokens",
	"pruning.wordEnd.likelihood",
	"pruning.wordIdentity.maxTokens",
	"pruning.wordIdentity.likelihood",
	"pruning.fanOut.maxTokens",
	"pruning.fanOut.likelihood",
	"log.file",
	"log.verbosityLevel",
	"insertionPenalty.standard",
	"insertionPenalty.filler",
	"adaptation.vtln",													// whether online VTLN will be applied
	"adaptation.vtln.estimation",										// VTL estimation mode (online/offline)
	"adaptation.vtln.warpFactor",										// warp factor to be applied to all the utterances (batch/online)
	"adaptation.vtln.warpFactor.file",								// file containing warp factors for each utterance ID (batch)
	"adaptation.mllr",													// whether online MLLR will be applied
	"adaptation.mllr.adaptCovariance",								// whether to adapt the covariance in addition to the mean
	"adaptation.mllr.regressionTree.file",							// file containing the regression tree
	"adaptation.mllr.regressionTree.minOccupation",				// minimum occupation of a tree-node to compute a transform
	"output.bestSinglePath",
	"output.bestSinglePath.outputSentenceDelimiters",
	"output.bestSinglePath.outputFillers",
	"output.bestSinglePath.outputConfidenceValues",
	"output.nBest",
	"output.nBest.n",
	"output.nBest.alignments.unique",
	"output.nBest.discardFillers",
	"output.nBest.singlePronunciation",
	"output.nBest.output.format",	
	"output.nBest.output.alignments",
	"output.nBest.output.likelihood",
	"output.nBest.output.folder",
	"output.graph",
	"output.graph.pruning.forwardBackward",	
	"output.graph.output.folder",
	"output.audio",
	"output.audio.output.folder",	
	"output.features",
	"output.features.output.folder",	
	"output.alignment",
	"output.alignment.output.folder",
	"confidenceEstimation.bestPath",
	"confidenceEstimation.method",
	"confidenceEstimation.graphPosteriors.insertionPenalty.standard",
	"confidenceEstimation.graphPosteriors.insertionPenalty.filler",
	"confidenceEstimation.graphPosteriors.scalingFactor.languageModel",
	"confidenceEstimation.graphPosteriors.confidenceMeasure",	
	};

int ConfigurationDecoder::m_iParametersRequired = 85;

// (3) read-only configuration parameters

// constructor
ConfigurationDecoder::ConfigurationDecoder(const char *strConfigurationFile) {

	m_strConfigurationFile = strConfigurationFile;
}

// destructor
ConfigurationDecoder::~ConfigurationDecoder() {
   
}

// set the default parameter values for the non-required parameters
void ConfigurationDecoder::setDefaultParameterValues() {

	return;
}

// check that the required parameters are defined in the configuration file


// load the configuration
bool ConfigurationDecoder::load() {

	setParametersAll(m_strParametersAll,m_iParametersAll);
	setParametersRequired(m_strParametersRequired,m_iParametersRequired);

	if (Configuration::load() == false) {
		return false;
	}
	
	if (checkParametersDefined() == false) {
		return false;
	}
   
   // get the parameters from the attribute-value pairs
      
	// task and mode
   getParameter("system.task",m_strSystemTask);
   getParameter("system.mode",m_strSystemMode);
   
   // decoding network
   getParameter("decodingNetwork.compact.forward",m_bDecodingNetworkCompactForward);
   getParameter("decodingNetwork.compact.backward",m_bDecodingNetworkCompactBackward);
   
   // input to the system
   getParameter("input.type",m_strInputType);
   
   // audio format
   getParameter("audio.format",m_strAudioFormat);
   getParameter("audio.samplingRate",m_iSamplingRate);
   getParameter("audio.sampleSize",m_iSampleSize);
   
   // speech activity detection
   getParameter("speechActivityDetection.autoEndPoint",m_bSpeechActivityDetectionAutoEndPoint);
   getParameter("speechActivityDetection.endPointPadding",m_iSpeechActivityDetectionEndPointPadding);
   getParameter("speechActivityDetection.penalty.silenceToSpeech",m_fSpeechActivityDetectionPenaltySilenceToSpeech);
   getParameter("speechActivityDetection.penalty.speechToSilence",m_fSpeechActivityDetectionPenaltySpeechToSilence);
   getParameter("speechActivityDetection.penalty.silence",m_fSpeechActivityDetectionPenaltySilence);
   getParameter("speechActivityDetection.minSpeechFrames",m_iSpeechActivityDetectionMinSpeechFrames);
   getParameter("speechActivityDetection.minSilenceFrames",m_iSpeechActivityDetectionMinSilenceFrames);
   getParameter("speechActivityDetection.speechMixtures",m_iSpeechActivityDetectionSpeechMixtures);
   
   // feature extraction
   getParameter("feature.configurationFile",m_strFileConfigurationFeatures);
   getParameter("feature.cepstralNormalization",m_bCepstralNormalization);
   getParameter("feature.cepstralNormalization.method",m_strCepstralNormalizationMethod);
   getParameter("feature.cepstralNormalization.bufferSize",m_iCepstralNormalizationBufferSize);
   if (getParameter("feature.transform.file",m_strFileFeatureTransform) == false) {
   	m_strFileFeatureTransform = "";
   }
   
   // phonetic symbol set
   getParameter("phoneticSymbolSet.file",m_strPhoneticSetFile);
	
	// acoustic models
	getParameter("acousticModels.file",m_strAcousticModelsFile);
	getParameter("acousticModels.format",m_strAcousticModelsFormat);
	getParameter("acousticModels.contextDependent",m_bAcousticModelsContextDependent);
	
	// transition probabilities
	getParameter("durationModeling",m_bDurationModeling);
	getParameter("durationModeling.scalingFactor",m_fDurationModelingScalingFactor);	
	
	// lexicon
	getParameter("lexicon.file",m_strLexiconFile);	
	getParameter("lexicon.alternativePronunciations",m_bLexiconAlternativePronunciations);	
		
	// language models
	getParameter("languageModel.file",m_strLanguageModelFile);
	getParameter("languageModel.format",m_strLanguageModelFormat);
	getParameter("languageModel.type",m_strLanguageModelType);
	getParameter("languageModel.lookAhead",m_bLanguageModelLookAhead);
	getParameter("languageModel.lookAhead.order",m_strLanguageModelLookAheadOrder);
	getParameter("languageModel.lookAhead.maxCacheSize",m_iLanguageModelLookAheadMaxCacheSize);
	getParameter("languageModel.scalingFactor",m_fLanguageModelScalingFactor);
	getParameter("languageModel.ngram",m_strLanguageModelNGram);
	getParameter("languageModel.crossUtterance",m_bLanguageModelCrossUtterance);
	
	// log
	getParameter("log.file",m_strLogFile);
	getParameter("log.verbosityLevel",m_iLogVerbosityLevel);	
	
	// beam pruning
	getParameter("pruning.global.maxTokens",m_iPruningGlobalMaxTokens);
	getParameter("pruning.global.likelihood",m_fPruningGlobalLikelihood);
	getParameter("pruning.wordEnd.maxTokens",m_iPruningWordEndMaxTokens);
	getParameter("pruning.wordEnd.likelihood",m_fPruningWordEndLikelihood);
	getParameter("pruning.wordIdentity.maxTokens",m_iPruningWordIdentityMaxTokens);
	getParameter("pruning.wordIdentity.likelihood",m_fPruningWordIdentityLikelihood);
	getParameter("pruning.fanOut.maxTokens",m_iPruningFanOutMaxTokens);
	getParameter("pruning.fanOut.likelihood",m_fPruningFanOutLikelihood);
	
	// insertion penalties
	getParameter("insertionPenalty.standard",m_fInsertionPenaltyStandard);
	getParameter("insertionPenalty.filler",m_fInsertionPenaltyFiller);
	getParameter("insertionPenalty.filler.file",m_strInsertionPenaltyFillerFile);	// optional
	
	// speaker adaptation (online)
	// (VTLN)
	getParameter("adaptation.vtln",m_bAdaptationVTLN);
	getParameter("adaptation.vtln.estimation",m_bAdaptationVTLNEstimation);
	getParameter("adaptation.vtln.warpFactor",m_fAdaptationVTLNWarpFactor);	
	getParameter("adaptation.vtln.warpFactor.file",m_strAdaptationVTLNWarpFactorFile);	
	// (MLLR)
	getParameter("adaptation.mllr",m_bAdaptationMLLR);
	getParameter("adaptation.mllr.adaptCovariance",m_bAdaptationMLLRAdaptCovariance);
	getParameter("adaptation.mllr.regressionTree.file",m_strAdaptationMLLRRegressionTreeFile);
	getParameter("adaptation.mllr.regressionTree.minOccupation",m_fAdaptationMLLRRegressionTreeMinOccupation);
	
	// decoder output
	// best single path
	getParameter("output.bestSinglePath",m_bOutputBestSinglePath);
	getParameter("output.bestSinglePath.outputSentenceDelimiters",m_bOutputBestSinglePathOutputSentenceDelimiters);
	getParameter("output.bestSinglePath.outputFillers",m_bOutputBestSinglePathOutputFillers);
	getParameter("output.bestSinglePath.outputConfidenceValues",m_bOutputBestSinglePathOutputConfidenceValues);
	// graph
	getParameter("output.graph",m_bOutputGraph);
	getParameter("output.graph.output.folder",m_strOutputGraphOutputFolder);
	// nbest
	getParameter("output.nBest",m_bOutputNBest);
	getParameter("output.nBest.n",m_iOutputNBestN);
	getParameter("output.nBest.alignments.unique",m_bOutputNBestAlignmentsUnique);
	getParameter("output.nBest.discardFillers",m_bOutputNBestDiscardFillers);
	getParameter("output.nBest.singlePronunciation",m_bOutputNBestSinglePronunciation);	
	getParameter("output.nBest.output.format",m_strOutputNBestOutputFormat);	
	getParameter("output.nBest.output.alignments",m_bOutputNBestOutputAlignments);
	getParameter("output.nBest.output.likelihood",m_bOutputNBestOutputLikelihood);	
	getParameter("output.nBest.output.folder",m_strOutputNBestOutputFolder);
	
	// audio
	getParameter("output.audio",m_bOutputAudio);
	getParameter("output.audio.output.folder",m_strOutputAudioOutputFolder);	
	// feature vectors
	getParameter("output.features",m_bOutputFeatureVectors);
	getParameter("output.features.output.folder",m_strOutputFeatureVectorsOutputFolder);
	// alignment
	getParameter("output.alignment",m_bOutputAlignment);
	getParameter("output.alignment.output.folder",m_strOutputAlignmentOutputFolder);
	
	// confidence estimation
	getParameter("confidenceEstimation.bestPath",m_bConfidenceEstimationBestPath);	
	getParameter("confidenceEstimation.method",m_strConfidenceEstimationMethod);	
	getParameter("confidenceEstimation.graphPosteriors.insertionPenalty.standard",m_fConfidenceEstimationGraphPosteriorsInsertionPenaltyStandard);
	getParameter("confidenceEstimation.graphPosteriors.insertionPenalty.filler",m_fConfidenceEstimationGraphPosteriorsInsertionPenaltyFiller);
	getParameter("confidenceEstimation.graphPosteriors.scalingFactor.languageModel",m_fConfidenceEstimationGraphPosteriorsScalingFactorLanguageModel);
	getParameter("confidenceEstimation.graphPosteriors.confidenceMeasure",m_strConfidenceEstimationGraphPosteriorsConfidenceMeasure);
	
	// check that the values of the read parameters are correct
	if (checkParameterValues() == true) {
		m_bConfigurationLoaded = true;
		return true;
	} else {
		m_bConfigurationLoaded = false;
		return false;
	}
}


// check that the values of the read parameters are correct
bool ConfigurationDecoder::checkParameterValues() {

	// imp: check that the maximum cache size is at least equal to the maximum number of tokens
	
	// batch file type
	if ((m_strInputType.compare(BATCH_FILE_TYPE_AUDIO) != 0) && (m_strInputType.compare(BATCH_FILE_TYPE_FEATURES) != 0) &&
	(m_strInputType.compare(BATCH_FILE_TYPE_AUDIO_SEGMENTS) != 0) &&
	(m_strInputType.compare(BATCH_FILE_TYPE_FEATURES_SEGMENTS) != 0)) {
		return false;
	}

	// feature extraction
	if (ConfigurationDecoder::m_iSampleSize != 16) {
		return false;
	}
	
	// language model parameters	
	
	if (m_strLanguageModelType.compare("ngram") != 0) {
		return false;
	}	
	
	if ((m_strLanguageModelNGram.compare("zerogram") != 0) && 
		(m_strLanguageModelNGram.compare("unigram") != 0) &&
		(m_strLanguageModelNGram.compare("bigram") != 0) &&
		(m_strLanguageModelNGram.compare("trigram") != 0)) {
		return false;
	}	
	
	// language model look-ahead
	if (m_strLanguageModelNGram.compare("zerogram") == 0) {
		// make sure that the language model look-ahead is disabled
		if (m_bLanguageModelLookAhead == true) {
			return false;
		}
	}
	
	// language model look-ahead order
	if ((m_strLanguageModelLookAheadOrder.compare("zerogram") != 0) && 
		(m_strLanguageModelLookAheadOrder.compare("unigram") != 0) &&
		(m_strLanguageModelLookAheadOrder.compare("bigram") != 0) &&
		(m_strLanguageModelLookAheadOrder.compare("trigram") != 0)) {
		return false;
	}	
	
	// n-best
	if ((m_strOutputNBestOutputFormat.compare(NBEST_OUTPUT_FORMAT_DEFAULT) != 0) && (m_strOutputNBestOutputFormat.compare(NBEST_OUTPUT_FORMAT_SRILM) != 0))	
	{
		return false;
	}
	
	// insertion penalty 
	/*if ((ConfigurationDecoder::m_fInsertionPenalty > 0) ||
		(ConfigurationDecoder::m_fInsertionPenaltySilence > 0) ||
		(ConfigurationDecoder::m_fInsertionPenaltyFiller > 0))	
	{
		return false;	
	}*/

	return true;
}

// set the value of a parameter (if the value is acceptable) overriding its current value
bool ConfigurationDecoder::setParameterValue(const char* strParameter, const char *strValue) {

	if ((strParameter == NULL) || (strValue == NULL)) {
		return false;
	}

	// task and mode
	if (strcmp(strParameter,"system.task") == 0) {
	
		if (setParameter(m_strSystemTask,strValue) == false) {
			return false;
		}
	} 
	else if (strcmp(strParameter,"system.mode") == 0) {
	
		if (setParameter(m_strSystemMode,strValue) == false) {
			return false;
		}	
	} 
   // decoding network	
	else if (strcmp(strParameter,"decodingNetwork.compact.forward") == 0) {
	
		if (setParameter(m_bDecodingNetworkCompactForward,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"decodingNetwork.compact.backward") == 0) {
	
		if (setParameter(m_bDecodingNetworkCompactBackward,strValue) == false) {
			return false;
		}	
	} 
	// audio format	
	else if (strcmp(strParameter,"audio.format") == 0) {
	
		if (setParameter(m_strAudioFormat,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"audio.samplingRate") == 0) {
	
		if (setParameter(m_iSamplingRate,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"audio.sampleSize") == 0) {
	
		if (setParameter(m_iSampleSize,strValue) == false) {
			return false;
		}	
	}
	// speech activity detection	 
	else if (strcmp(strParameter,"speechActivityDetection.autoEndPoint") == 0) {
	
		if (setParameter(m_bSpeechActivityDetectionAutoEndPoint,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.endPointPadding") == 0) {
	
		if (setParameter(m_iSpeechActivityDetectionEndPointPadding,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.penalty.silenceToSpeech") == 0) {
	
		if (setParameter(m_fSpeechActivityDetectionPenaltySilenceToSpeech,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.penalty.speechToSilence") == 0) {
	
		if (setParameter(m_fSpeechActivityDetectionPenaltySpeechToSilence,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.penalty.silence") == 0) {
	
		if (setParameter(m_fSpeechActivityDetectionPenaltySilence,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.minSpeechFrames") == 0) {
	
		if (setParameter(m_iSpeechActivityDetectionMinSpeechFrames,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.minSilenceFrames") == 0) {
	
		if (setParameter(m_iSpeechActivityDetectionMinSilenceFrames,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"speechActivityDetection.speechMixtures") == 0) {
	
		if (setParameter(m_iSpeechActivityDetectionSpeechMixtures,strValue) == false) {
			return false;
		}	
	} 
   // feature extraction
	else if (strcmp(strParameter,"feature.cepstralNormalization") == 0) {
	
		if (setParameter(m_bCepstralNormalization,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"feature.cepstralNormalization.method") == 0) {
	
		if (setParameter(m_strCepstralNormalizationMethod,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"feature.cepstralNormalization.bufferSize") == 0) {
	
		if (setParameter(m_iCepstralNormalizationBufferSize,strValue) == false) {
			return false;
		}	
	} 
	// phonetic symbol set
	else if (strcmp(strParameter,"phoneticSymbolSet.file") == 0) {
	
		if (setParameter(m_strPhoneticSetFile,strValue) == false) {
			return false;
		}	
	} 
	// acoustic models
	else if (strcmp(strParameter,"acousticModels.file") == 0) {
	
		if (setParameter(m_strAcousticModelsFile,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"acousticModels.format") == 0) {
	
		if (setParameter(m_strAcousticModelsFormat,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"acousticModels.contextDependent") == 0) {
	
		if (setParameter(m_bAcousticModelsContextDependent,strValue) == false) {
			return false;
		}	
	} 
	// transition probabilities
	else if (strcmp(strParameter,"durationModeling") == 0) {
	
		if (setParameter(m_bDurationModeling,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"durationModeling.scalingFactor") == 0) {
	
		if (setParameter(m_fDurationModelingScalingFactor,strValue) == false) {
			return false;
		}	
	} 
	// lexicon	
	else if (strcmp(strParameter,"lexicon.file") == 0) {
	
		if (setParameter(m_strLexiconFile,strValue) == false) {
			return false;
		}	
	} 
	// language models
	else if (strcmp(strParameter,"languageModel.file") == 0) {
	
		if (setParameter(m_strLanguageModelFile,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"languageModel.format") == 0) {
	
		if (setParameter(m_strLanguageModelFormat,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"languageModel.type") == 0) {
	
		if (setParameter(m_strLanguageModelType,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"languageModel.lookAhead") == 0) {
	
		if (setParameter(m_bLanguageModelLookAhead,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"languageModel.scalingFactor") == 0) {
	
		if (setParameter(m_fLanguageModelScalingFactor,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"languageModel.ngram") == 0) {
	
		if (setParameter(m_strLanguageModelNGram,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"languageModel.crossUtterance") == 0) {
	
		if (setParameter(m_bLanguageModelCrossUtterance,strValue) == false) {
			return false;
		}	
	} 
	// log
	else if (strcmp(strParameter,"log.file") == 0) {
	
		if (setParameter(m_strLogFile,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"log.verbosityLevel") == 0) {
	
		if (setParameter(m_iLogVerbosityLevel,strValue) == false) {
			return false;
		}	
	} 
	// beam pruning
	else if (strcmp(strParameter,"pruning.global.maxTokens") == 0) {
	
		if (setParameter(m_iPruningGlobalMaxTokens,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.global.likelihood") == 0) {
	
		if (setParameter(m_fPruningGlobalLikelihood,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.wordEnd.maxTokens") == 0) {
	
		if (setParameter(m_iPruningWordEndMaxTokens,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.wordEnd.likelihood") == 0) {
	
		if (setParameter(m_fPruningWordEndLikelihood,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.wordIdentity.maxTokens") == 0) {
	
		if (setParameter(m_iPruningWordIdentityMaxTokens,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.wordIdentity.likelihood") == 0) {
	
		if (setParameter(m_fPruningWordIdentityLikelihood,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.fanOut.maxTokens") == 0) {
	
		if (setParameter(m_iPruningFanOutMaxTokens,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"pruning.fanOut.likelihood") == 0) {
	
		if (setParameter(m_fPruningFanOutLikelihood,strValue) == false) {
			return false;
		}	
	} 
	// insertion penalties
	else if (strcmp(strParameter,"insertionPenalty.standard") == 0) {
	
		if (setParameter(m_fInsertionPenaltyStandard,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"insertionPenalty.filler") == 0) {
	
		if (setParameter(m_fInsertionPenaltyFiller,strValue) == false) {
			return false;
		}	
	} 
	// decoder output
	else if (strcmp(strParameter,"output.bestSinglePath") == 0) {
	
		if (setParameter(m_bOutputBestSinglePath,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"output.bestSinglePath.outputFillers") == 0) {
	
		if (setParameter(m_bOutputBestSinglePathOutputFillers,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"output.nBest") == 0) {
	
		if (setParameter(m_bOutputNBest,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"output.nBest.output.folder") == 0) {
	
		if (setParameter(m_strOutputNBestOutputFolder,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"output.graph") == 0) {
	
		if (setParameter(m_bOutputGraph,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"output.graph.outputFolder") == 0) {
	
		if (setParameter(m_strOutputGraphOutputFolder,strValue) == false) {
			return false;
		}	
	} 
	// confidence estimation
	else if (strcmp(strParameter,"confidenceEstimation.bestPath") == 0) {
	
		if (setParameter(m_bConfidenceEstimationBestPath,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"confidenceEstimation.method") == 0) {
	
		if (setParameter(m_strConfidenceEstimationMethod,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"confidenceEstimation.graphPosteriors.insertionPenalty.standard") == 0) {
	
		if (setParameter(m_fConfidenceEstimationGraphPosteriorsInsertionPenaltyStandard,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"confidenceEstimation.graphPosteriors.insertionPenalty.filler") == 0) {
	
		if (setParameter(m_fConfidenceEstimationGraphPosteriorsInsertionPenaltyFiller,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"confidenceEstimation.graphPosteriors.scalingFactor.languageModel") == 0) {
	
		if (setParameter(m_fConfidenceEstimationGraphPosteriorsScalingFactorLanguageModel,strValue) == false) {
			return false;
		}	
	} 
	else if (strcmp(strParameter,"confidenceEstimation.graphPosteriors.confidenceMeasure") == 0) {
	
		if (setParameter(m_strConfidenceEstimationGraphPosteriorsConfidenceMeasure,strValue) == false) {
			return false;
		}	
	} else {
		return false;
	}

	return true;
}
