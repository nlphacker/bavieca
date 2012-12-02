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


#ifndef CONFIGURATIONBAVIECA_H
#define CONFIGURATIONBAVIECA_H

#include "ConfigurationFile.h"
#include "ParameterManager.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ConfigurationBavieca : public ParameterManager, public ConfigurationFile {

	public:

		// constructor
		ConfigurationBavieca() : ParameterManager(), ConfigurationFile(NULL) {
		}

		// constructor (from a file)
		ConfigurationBavieca(const char *strFile) : ParameterManager(), ConfigurationFile(strFile) {
		}

		// destructor
		~ConfigurationBavieca() {
		}
		
		// define the parameters
		void defineParameters();
		
		// load the configuration parameters
		void load();
		
		// check whether SAD configuration parameters are set
		bool areSADParametersSet();
				
		// check whether aligner configuration parameters are set
		bool areAlignerParametersSet();
				
		// check whether decoding configuration parameters are set
		bool areDecodingParametersSet();		

};

};	// end-of-namespace

#endif
