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


#ifndef CONFIGURATIONDYNAMICDECODER_H
#define CONFIGURATIONDYNAMICDECODER_H

#include "ConfigurationFile.h"
#include "ParameterManager.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ConfigurationDynamicDecoder : public ParameterManager, public ConfigurationFile {

	private:
	
		// define the parameters
		void defineParameters();	

	public:
		
		// constructor
		ConfigurationDynamicDecoder();

		// constructor (from a file)
		ConfigurationDynamicDecoder(const char *strFile);

		// destructor
		~ConfigurationDynamicDecoder();
		
		// load the configuration parameters
		void load();
};

};	// end-of-namespace

#endif
