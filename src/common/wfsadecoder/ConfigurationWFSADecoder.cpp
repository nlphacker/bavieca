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


#include "ConfigurationWFSADecoder.h"

namespace Bavieca {

// constructor
ConfigurationWFSADecoder::ConfigurationWFSADecoder(const char *strFile) : ParameterManager(),ConfigurationFile(strFile)
{
	m_strFile = strFile;
}

// destructor
ConfigurationWFSADecoder::~ConfigurationWFSADecoder()
{
}

// define the parameters
void ConfigurationWFSADecoder::defineParameters() {

	// decoding network
	defineParameter("decodingNetwork.file","decoding network",PARAMETER_TYPE_FILE,false);
	
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
	
	// lexicon
	defineParameter("lexicon.file","",PARAMETER_TYPE_FILE,false);

	// pruning 
	defineParameter("pruning.maxActiveStates","maximum number of active states after pruning",
		PARAMETER_TYPE_INTEGER,false,"[100|1000000]");
	defineParameter("pruning.likelihoodBeam","likelihood beam used for pruning states",
		PARAMETER_TYPE_FLOAT,false,"[10.0|1000000.0]");
	
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
void ConfigurationWFSADecoder::load() {

	try {

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
		
	} catch (ExceptionBase &e) {
		BVC_ERROR << "unable to load the wfsa configuration file";
	}
}

};	// end-of-namespace



