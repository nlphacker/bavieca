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

#ifndef CONFIGURATIONDECODER_H
#define CONFIGURATIONDECODER_H

#include "Configuration.h"

#define BATCH_FILE_TYPE_AUDIO						"audio"
#define BATCH_FILE_TYPE_AUDIO_SEGMENTS			"audioSegment"
#define BATCH_FILE_TYPE_FEATURES					"featureVectors"
#define BATCH_FILE_TYPE_FEATURES_SEGMENTS		"featureVectorsSegment"

#define ACOUSTIC_MODELS_FORMAT_DEFAULT		"default"
#define ACOUSTIC_MODELS_FORMAT_SONIC		"Sonic"

#define NBEST_OUTPUT_FORMAT_DEFAULT			"default"
#define NBEST_OUTPUT_FORMAT_SRILM			"SRILM"

// ceptral normalization methods
#define CEPSTRAL_MEAN_NORMALIZATION						"CMN"
#define CEPSTRAL_MEAN_VARIANCE_NORMALIZATION			"CMVN"

/**
	@author root <root@localhost.localdomain>
*/
class ConfigurationDecoder : public Configuration {

	private:
	
		static int m_iParametersAll;
		static const char *m_strParametersAll[];
		static int m_iParametersRequired;
		static const char *m_strParametersRequired[];	
		
		// set the default parameter values for the non-required parameters
		void setDefaultParameterValues();
		
		// check that the values of the read parameters are correct
		bool checkParameterValues();	
		
	public:
		
		// parameters
	
   	int m_iCores;	// number of cores in the system (for open-mp based optimization)	
	
   	// system
   	string m_strSystemTask;
   	string m_strSystemMode;
   	
   	// decoding network
   	bool m_bDecodingNetworkCompactForward;
   	bool m_bDecodingNetworkCompactBackward;
   	
   	// input to the system
		string m_strInputType;
   	
   	// audio format
   	string m_strAudioFormat;		// [ raw ]
   	int m_iSamplingRate;			// sampling rate in Khz
   	int m_iSampleSize;				// bits per sample
   	
		// speech activity detection
		bool m_bSpeechActivityDetectionAutoEndPoint;
		int m_iSpeechActivityDetectionEndPointPadding;
		float m_fSpeechActivityDetectionPenaltySilenceToSpeech;
		float m_fSpeechActivityDetectionPenaltySpeechToSilence;
		float m_fSpeechActivityDetectionPenaltySilence;
		int m_iSpeechActivityDetectionMinSpeechFrames;
		int m_iSpeechActivityDetectionMinSilenceFrames;
		int m_iSpeechActivityDetectionSpeechMixtures;	
   	
   	// feature extraction
   	string m_strFileConfigurationFeatures;
   	bool m_bCepstralNormalization;	// whether to perform cepstral normalization	
   	string m_strCepstralNormalizationMethod;
   	int m_iCepstralNormalizationBufferSize;
   	string m_strFileFeatureTransform;
   
   	// acoustic models
   	string m_strAcousticModelsFile;
   	string m_strAcousticModelsFormat;
   	bool m_bAcousticModelsContextDependent;
   	FILE *m_fileAcousticModels;
   	int m_iTopMixturesCount;
   	
   	// transition probabilities
		bool m_bDurationModeling;
		float m_fDurationModelingScalingFactor;	
   	
   	// language model
   	string m_strLanguageModelFile;					// file containing the language model
   	string m_strLanguageModelFormat;					// language model format: ARPA, binary
   	string m_strLanguageModelType;					// language model type: ( ngram | CFG ) 
   	bool m_bLanguageModelLookAhead;					// whether to apply language model look-ahead
   	string m_strLanguageModelLookAheadOrder;		// language model look-ahead order (unigram/bigram/trigram)
   	int m_iLanguageModelLookAheadMaxCacheSize;	// whether to apply language model look-ahead
   	float m_fLanguageModelScalingFactor;			// language model scaling factor
		string m_strLanguageModelNGram;					// unigram | bigram | trigram | ...
		bool m_bLanguageModelCrossUtterance;			// cross utterance language model
   	
   	// lexicon
   	string m_strLexiconFile;
   	bool m_bLexiconAlternativePronunciations;
   	
   	// phonetic symbol set
		string m_strPhoneticSetFile;
   
      // log
      string m_strLogFile;
      int m_iLogVerbosityLevel; 
      
		// pruning
		int m_iPruningGlobalMaxTokens;	
		float m_fPruningGlobalLikelihood;
		
		int m_iPruningWordEndMaxTokens;
		float m_fPruningWordEndLikelihood;
		
		int m_iPruningWordIdentityMaxTokens;
		float m_fPruningWordIdentityLikelihood;
		
		int m_iPruningFanOutMaxTokens;
		float m_fPruningFanOutLikelihood;	
		
   	// insertion penalties
   	float m_fInsertionPenaltyStandard;
   	float m_fInsertionPenaltyFiller;
   	string m_strInsertionPenaltyFillerFile;  	
   	
		// speaker adaptation (online)
		// (VTLN)
   	bool m_bAdaptationVTLN;
   	bool m_bAdaptationVTLNEstimation;
		float m_fAdaptationVTLNWarpFactor;
		string m_strAdaptationVTLNWarpFactorFile;
   	
   	// (MLLR)
   	bool m_bAdaptationMLLR;
		bool m_bAdaptationMLLRAdaptCovariance;
   	string m_strAdaptationMLLRRegressionTreeFile;
   	float m_fAdaptationMLLRRegressionTreeMinOccupation;
   	
   	// decoder output
   	// best single path
		bool m_bOutputBestSinglePath;
		bool m_bOutputBestSinglePathOutputSentenceDelimiters;
		bool m_bOutputBestSinglePathOutputFillers;
		bool m_bOutputBestSinglePathOutputConfidenceValues;
		// graph
		bool m_bOutputGraph;
		string m_strOutputGraphOutputFolder;
		// n-best
		bool m_bOutputNBest;
		int m_iOutputNBestN;
		bool m_bOutputNBestAlignmentsUnique;
		bool m_bOutputNBestDiscardFillers;
		bool m_bOutputNBestSinglePronunciation;
		string m_strOutputNBestOutputFormat;
		bool m_bOutputNBestOutputAlignments;
		bool m_bOutputNBestOutputLikelihood;
		string m_strOutputNBestOutputFolder;	
		// audio
		bool m_bOutputAudio;
		string m_strOutputAudioOutputFolder;
		// feature vectors
		bool m_bOutputFeatureVectors;
		string m_strOutputFeatureVectorsOutputFolder;	
		// alignment
		bool m_bOutputAlignment;
		string m_strOutputAlignmentOutputFolder;
		
		// confidence estimation
		bool m_bConfidenceEstimationBestPath;
		string m_strConfidenceEstimationMethod;	
		float m_fConfidenceEstimationGraphPosteriorsInsertionPenaltyStandard;
		float m_fConfidenceEstimationGraphPosteriorsInsertionPenaltyFiller;
		float m_fConfidenceEstimationGraphPosteriorsScalingFactorLanguageModel;
		string m_strConfidenceEstimationGraphPosteriorsConfidenceMeasure;
   		
		// constructor
		ConfigurationDecoder(const char *strConfigurationFile);

		// destructor
		~ConfigurationDecoder();		
		
		// load the configuration
		bool load();
		
		// set the value of a parameter (if the value is acceptable)
		bool setParameterValue(const char* strParameter, const char *strValue);	

};

#endif
