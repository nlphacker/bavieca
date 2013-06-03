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


#include "ConfigurationDynamicDecoder.h"

namespace Bavieca {

// constructor
ConfigurationDynamicDecoder::ConfigurationDynamicDecoder() : 
	ParameterManager(), ConfigurationFile(NULL)
{
}

// constructor (from a file)
ConfigurationDynamicDecoder::ConfigurationDynamicDecoder(const char *strFile) : 
	ParameterManager(), ConfigurationFile(strFile)
{
}

// destructor
ConfigurationDynamicDecoder::~ConfigurationDynamicDecoder()
{
}

// define the parameters
void ConfigurationDynamicDecoder::defineParameters() {

	// input
	defineParameter("input.type","format of input data",PARAMETER_TYPE_STRING,false,"audio|features");
	
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
	
	// language model
	defineParameter("languageModel.file","language model",PARAMETER_TYPE_FILE,false);
	defineParameter("languageModel.format","language model format",PARAMETER_TYPE_STRING,false);
	defineParameter("languageModel.type","language model type",PARAMETER_TYPE_STRING,false);
	defineParameter("languageModel.scalingFactor","language model scaling factor",PARAMETER_TYPE_FLOAT,false);
	defineParameter("languageModel.crossUtterance","language model",PARAMETER_TYPE_BOOLEAN,true,"yes|no","no");
	
	// lexicon
	defineParameter("lexicon.file","pronunciation lexicon",PARAMETER_TYPE_FILE,false);
	
	// insertion penalty
	defineParameter("insertionPenalty.standard","",PARAMETER_TYPE_FLOAT,false);
	defineParameter("insertionPenalty.filler","",PARAMETER_TYPE_FLOAT,false);
	defineParameter("insertionPenalty.filler.file","",PARAMETER_TYPE_FILE,false);

	// pruning 
	defineParameter("pruning.likelihoodBeam","likelihood beam used for pruning arcs",
		PARAMETER_TYPE_FLOAT,false);
	defineParameter("pruning.likelihoodBeamWE","likelihood beam used for pruning word-end arcs",
		PARAMETER_TYPE_FLOAT,false);
	defineParameter("pruning.likelihoodBeamTokensArc","likelihood beam used for pruning tokens within an arc",
		PARAMETER_TYPE_FLOAT,false);
	defineParameter("pruning.maxActiveArcs","maximum number of active arcs",
		PARAMETER_TYPE_INTEGER,false);
	defineParameter("pruning.maxActiveArcsWE","maximum number of active word-end arcs",
		PARAMETER_TYPE_INTEGER,false);
	defineParameter("pruning.maxActiveTokensArc","maximum number of active tokens per arc",
		PARAMETER_TYPE_INTEGER,false);
	
	// decoder output
	defineParameter("output.bestSinglePath","",PARAMETER_TYPE_BOOLEAN,true,"yes|no","yes");
	defineParameter("output.lattice.folder","",PARAMETER_TYPE_FOLDER,true);
	defineParameter("output.lattice.maxWordSequencesState",
		"maximum number of different word sequences exiting a state",
		PARAMETER_TYPE_INTEGER,true,"[2|100]","5");
	defineParameter("output.audio.folder","folder to store the audio",PARAMETER_TYPE_FOLDER,true);
	defineParameter("output.features.folder","folder to store the features",PARAMETER_TYPE_FOLDER,true);
	defineParameter("output.alignment.folder","folder to store the alignment",PARAMETER_TYPE_FOLDER,true);
}

// load the configuration parameters
void ConfigurationDynamicDecoder::load() {

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

};	// end-of-namespace
