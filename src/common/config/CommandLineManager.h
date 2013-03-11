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


#ifndef COMMANDLINEMANAGER_H
#define COMMANDLINEMANAGER_H

#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "ParameterManager.h"

using namespace std;

#include <string>
#include <map>
#include <vector>

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class CommandLineManager : public ParameterManager {
	
	private:
	
		string m_strApplicationName;
		string m_strVersion;
		string m_strAuthor;
		string m_strDate;	
				
		// return the list of items from a string of characters of the form "item1|item2|...|itemN"
		VString *getItemList(const char *strItems);
		
		// checks whether a string of characters represents a floating point number
		bool isFloat(const char *strSource) {
		
			char *strAux = (char*)strSource;	
			/*float fValue = (float)*/strtod(strSource,&strAux);
			if ((errno != 0) ||					// conversion failed (EINVAL, ERANGE)
				(strSource == strAux) ||		// conversion failed (no characters consumed)
				(*strAux != 0))					// conversion failed (trailing data)
				{
				return false;	
			}	
			
			return true;
		}

		// checks whether a string of characters represents an integer
		bool isInteger(const char *strSource) {
		
			char *strAux = (char*)strSource;	
			/*long int lValue = */strtol(strSource,&strAux,10);
			if ((errno != 0) ||					// conversion failed (EINVAL, ERANGE)
				(strSource == strAux) ||		// conversion failed (no characters consumed)
				(*strAux != 0))					// conversion failed (trailing data)
				{
				return false;	
			}	
			
			return true;
		}
		
		// checks whether a string of characters represents an integer
		bool isBoolean(const char *strSource) {
		
			if ((strcmp(strSource,"yes") == 0) || (strcmp(strSource,"no") == 0)) {
				return true;
			}
			
			return false;
		}
		
		// return the parameter type in string format
		const char *getStrParameterType(char iType) {
		
			switch(iType) {
				case PARAMETER_TYPE_STRING: {
					return STR_PARAMETER_TYPE_STRING;
				}
				case PARAMETER_TYPE_BOOLEAN: {
					return STR_PARAMETER_TYPE_BOOLEAN;
				}
				case PARAMETER_TYPE_INTEGER: {
					return STR_PARAMETER_TYPE_INTEGER;
				}
				case PARAMETER_TYPE_FLOAT: {
					return STR_PARAMETER_TYPE_FLOAT;
				}
				case PARAMETER_TYPE_FILE: {
					return STR_PARAMETER_TYPE_FILE;
				}
				case PARAMETER_TYPE_FOLDER: {
					return STR_PARAMETER_TYPE_FOLDER;
				}
				default: {
					assert(0);
					return NULL;
				}
			}
		}
		
		// get the parameter type as a string
		void getStrType(Parameter *parameter, string &strType);

	public:
	
		// constructor
		CommandLineManager(const char *strApplicationName, const char *strVersion, const char *strAuthor, const char *strDate);

		// destructor
		~CommandLineManager();		
		
		// display the application usage including parameter sintax
		void displayApplicationUsage();
		
		// parse parameters read from the command line
		bool parseParameters(int argc, char *argv[]);
		
		// print parameters
		void printParameters();	
};

};	// end-of-namespace

#endif
