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


#ifndef CONFIGURATIONFILE_H
#define CONFIGURATIONFILE_H

// string lengths for both atributtes and values
#define MAX_CONFIGURATION_ATTRIBUTE_LENGTH		100
#define MAX_CONFIGURATION_VALUE_LENGTH				900

#include <stdlib.h>
#include <string.h>

using namespace std;

#include <string>
#include <map>
#include <vector>

#include "Global.h"

namespace Bavieca {

typedef map<string,string> MParameterValue;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ConfigurationFile { 

	protected:
	
		string m_strFile;
		MParameterValue m_mParameterValue;	

	public:

		// constructor
		ConfigurationFile(const char *strFile);

		// destructor
		~ConfigurationFile();
		
		// load the configuration file
		void load();
		
		// get the configuration data
		inline MParameterValue *getData() {
		
			return &m_mParameterValue;
		}
		
		// print the parameters
		void print();
		
		// return whether a string is blank
		bool isBlank(string &string) {
		
			for(unsigned int i=0 ; i<string.length() ; ++i) {
				if (isspace(string.c_str()[i]) == 0) {
					return false;
				}
			}
			
			return true;
		}

};

};	// end-of-namespace

#endif
