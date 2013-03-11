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


#include "ConfigurationBavieca.h"

namespace Bavieca {

// define the parameters
void ConfigurationBavieca::defineParameters() {

	// feature extraction	
	defineParameter("feature.configurationFile","",PARAMETER_TYPE_FILE,false);
	defineParameter("feature.cepstralNormalization.mode","cepstral normalization mode",
		PARAMETER_TYPE_STRING,false,"utterance|stream|session");
	defineParameter("feature.cepstralNormalization.method","cepstral normalization method",
		PARAMETER_TYPE_STRING,false,"none|CMN|CMVN");
	defineParameter("feature.cepstralNormalization.bufferSize","cepstral buffer size (in bytes)",
		PARAMETER_TYPE_INTEGER,false);
	defineParameter("feature.warpFactor","warp factor for VTLN",
		PARAMETER_TYPE_FLOAT,false);
	defineParameter("feature.transformFile","list of feature transforms",
		PARAMETER_TYPE_FILE,true); 
		
	// phonetic symbol set
	defineParameter("phoneticSymbolSet.file","set of phonetic symbols",PARAMETER_TYPE_FILE,false);
	
	// acoustic models
	defineParameter("acousticModels.file","acoustic models",PARAMETER_TYPE_FILE,false);
	
	// speech activity detection
	defineParameter("sad.maxGaussianSilence","max Gaussian components for silence model",PARAMETER_TYPE_INTEGER,true);
	defineParameter("sad.maxGaussianSpeech","max Gaussian components for speech model",PARAMETER_TYPE_INTEGER,true);
	defineParameter("sad.speechPadding","speech padding",PARAMETER_TYPE_INTEGER,true);
	defineParameter("sad.speechPenalty","speech insertion penalty",PARAMETER_TYPE_FLOAT,true);
	
	// lexicon
	defineParameter("lexicon.file","pronunciation lexicon",PARAMETER_TYPE_FILE,true);
	
	// language model
	defineParameter("languageModel.file","language model",PARAMETER_TYPE_FILE,true);
	defineParameter("languageModel.format","language model format",PARAMETER_TYPE_STRING,true);
	defineParameter("languageModel.type","language model type",PARAMETER_TYPE_STRING,true);
	defineParameter("languageModel.scalingFactor","language model scaling factor",PARAMETER_TYPE_FLOAT,true);
	defineParameter("languageModel.ngram","language model n-gram",PARAMETER_TYPE_STRING,true);
	defineParameter("languageModel.crossUtterance","language model",PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");
	
	// insertion penalty
	defineParameter("insertionPenalty.standard","",PARAMETER_TYPE_FLOAT,true);
	defineParameter("insertionPenalty.filler","",PARAMETER_TYPE_FLOAT,true);
	defineParameter("insertionPenalty.filler.file","",PARAMETER_TYPE_FILE,true);

	// pruning 
	defineParameter("pruning.likelihoodBeam","likelihood beam used for pruning arcs",
		PARAMETER_TYPE_FLOAT,true);
	defineParameter("pruning.likelihoodBeamWE","likelihood beam used for pruning word-end arcs",
		PARAMETER_TYPE_FLOAT,true);
	defineParameter("pruning.likelihoodBeamTokensArc","likelihood beam used for pruning tokens within an arc",
		PARAMETER_TYPE_FLOAT,true);
	defineParameter("pruning.maxActiveArcs","maximum number of active arcs",
		PARAMETER_TYPE_INTEGER,true);
	defineParameter("pruning.maxActiveArcsWE","maximum number of active word-end arcs",
		PARAMETER_TYPE_INTEGER,true);
	defineParameter("pruning.maxActiveTokensArc","maximum number of active tokens per arc",
		PARAMETER_TYPE_INTEGER,true);
		
	// decoder output
	defineParameter("output.lattice.maxWordSequencesState",
		"maximum number of different word sequences exiting a state",
		PARAMETER_TYPE_INTEGER,true,"[2|100]","5");
}

// load the configuration parameters
void ConfigurationBavieca::load() {

	// load pairs parameter-value from the file
	ConfigurationFile::load();

	// get the data
	VParameterValue vParameterValue;
	for(MParameterValue::iterator it = m_mParameterValue.begin() ; it != m_mParameterValue.end() ; ++it) {
		ParameterValue parameterValue;
		parameterValue.strParameter = it->first;
		parameterValue.strValue = it->second;
		vParameterValue.push_back(parameterValue);
	}
	
	//ConfigurationFile::print();
	
	// make the parameter definitiones
	defineParameters();
	
	// parse the parameter-value pairs
	parse(vParameterValue);
}

// check whether SAD configuration parameters are set
bool ConfigurationBavieca::areSADParametersSet() {

	// make sure required configuration parameters are defined
	if (isParameterSet("sad.maxGaussianSilence") &&
		isParameterSet("sad.maxGaussianSpeech") &&
		isParameterSet("sad.speechPadding") &&
		isParameterSet("sad.speechPenalty")) {
		return true;
	}
	
	return false;
}

// check whether alignment configuration parameters are set
bool ConfigurationBavieca::areAlignerParametersSet() {

	// make sure required configuration parameters are defined
	if (isParameterSet("lexicon.file")) {
		return true;
	}
	
	return false;
}

// check whether decoding configuration parameters are set
bool ConfigurationBavieca::areDecodingParametersSet() {

	// make sure required configuration parameters are defined
	if (isParameterSet("languageModel.file") &&
		isParameterSet("languageModel.format") &&
		isParameterSet("languageModel.type") &&
		isParameterSet("languageModel.scalingFactor") &&
		isParameterSet("languageModel.ngram") &&
		isParameterSet("languageModel.crossUtterance") &&
		isParameterSet("insertionPenalty.standard") &&
		isParameterSet("insertionPenalty.filler") &&
		isParameterSet("insertionPenalty.filler.file") &&
		isParameterSet("pruning.likelihoodBeam") &&
		isParameterSet("pruning.likelihoodBeamWE") &&
		isParameterSet("pruning.likelihoodBeamTokensArc") &&
		isParameterSet("pruning.maxActiveArcs") &&
		isParameterSet("pruning.maxActiveArcsWE") &&
		isParameterSet("pruning.maxActiveTokensArc")) {
		return true;
	}
	
	return false;
}

};	// end-of-namespace

