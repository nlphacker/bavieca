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

#include "ConfigurationFeatures.h"
#include "FeatureExtractor.h"

namespace Bavieca {

// constructor
ConfigurationFeatures::ConfigurationFeatures() : 
	ParameterManager(), ConfigurationFile(NULL)
{
}

// constructor (from a file)
ConfigurationFeatures::ConfigurationFeatures(const char *strConfigurationFile) :
	ParameterManager(), ConfigurationFile(strConfigurationFile)
{
}

// destructor
ConfigurationFeatures::~ConfigurationFeatures()
{

}

// define the parameters
void ConfigurationFeatures::defineParameters() {

	// waveform
	defineParameter("waveform.samplingRate","audio sampling rate (in Hz)",PARAMETER_TYPE_INTEGER,false,"8000|16000","16000");
	defineParameter("waveform.sampleSize","audio sample size (in bits)",PARAMETER_TYPE_INTEGER,false,"8|16","16");
	
	// DC-removal
	defineParameter("dcRemoval","DC removal",PARAMETER_TYPE_BOOLEAN,false,NULL,"yes");
	
	// preemphasis
	defineParameter("preemphasis","preemphasis",PARAMETER_TYPE_BOOLEAN,false,NULL,"yes");
	
	// windowing
	defineParameter("window.width","window width (in ms)",PARAMETER_TYPE_INTEGER,false,NULL,"20");
	defineParameter("window.shift","windows shift (in ms)",PARAMETER_TYPE_INTEGER,false,NULL,"10");
	defineParameter("window.tapering","window tapering function",
		PARAMETER_TYPE_STRING,false,"Hann|HannModified|Hamming|none","Hamming");
		
	// type
	defineParameter("features.type","feature type",PARAMETER_TYPE_STRING,false,"mfcc|plp","mfcc");
	
	// filterbank analysis
	defineParameter("filterbank.frequency.min","minimum frequency (in Hz)",PARAMETER_TYPE_INTEGER,false,NULL,"0");
	defineParameter("filterbank.frequency.max","maximum frequency (in Hz)",PARAMETER_TYPE_INTEGER,false,NULL,"8000");
	defineParameter("filterbank.filters","number of filters",PARAMETER_TYPE_INTEGER,false,NULL,"20");
	
	// cepstral coefficients
	defineParameter("cepstralCoefficients","number of cepstral coefficients",
		PARAMETER_TYPE_INTEGER,false,NULL,"12");	
	
	// energy
	defineParameter("energy","energy",PARAMETER_TYPE_BOOLEAN,false,NULL,"yes");	
	
	// derivatives
	defineParameter("derivatives.order","derivatives order",PARAMETER_TYPE_INTEGER,true,NULL,"2");
	defineParameter("derivatives.delta","derivatives delta",PARAMETER_TYPE_INTEGER,true,NULL,"2");	
	
	// spliced features
	defineParameter("spliced.size","number of consecutive static vectors that will be spliced",
		PARAMETER_TYPE_INTEGER,true,NULL,"9");
}

// load the configuration parameters
void ConfigurationFeatures::load() {

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
		
		// make the parameter definitions
		defineParameters();
		
		// parse the parameter-value pairs
		parse(vParameterValue);
		
	} catch (std::runtime_error) {
		BVC_ERROR << "unable to load the feature configuration file";
	}
}

// return the feature dimensionality
int ConfigurationFeatures::getDimensionality() {

	int iStatic = atoi(getParameterValue("cepstralCoefficients"));
	if (strcmp(getParameterValue("energy"),"yes") == 0) {
		iStatic += 1;
	}
	
	if (isParameterSet("derivatives.order")) {
		return iStatic*(atoi(getParameterValue("derivatives.order"))+1);
	} else {
		assert(isParameterSet("spliced.size"));
		return iStatic*atoi(getParameterValue("spliced.size"));
	}	
}

};	// end-of-namespace
