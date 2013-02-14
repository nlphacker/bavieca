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


#ifndef PARAMETERMANAGER_H
#define PARAMETERMANAGER_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

#include "LogMessage.h"

using namespace std;

#include <string>
#include <map>
#include <vector>

namespace Bavieca {

#define MAX_PARAMETER_NAME_LENGTH		100

// parameter types
#define PARAMETER_TYPE_STRING				0
#define PARAMETER_TYPE_BOOLEAN			1
#define PARAMETER_TYPE_INTEGER			2
#define PARAMETER_TYPE_FLOAT				3
#define PARAMETER_TYPE_FILE				4
#define PARAMETER_TYPE_FOLDER				5

// parameter types (string format)
#define STR_PARAMETER_TYPE_STRING				"string"
#define STR_PARAMETER_TYPE_BOOLEAN				"yes|no"
#define STR_PARAMETER_TYPE_INTEGER				"int"
#define STR_PARAMETER_TYPE_FLOAT					"float"
#define STR_PARAMETER_TYPE_FILE					"file"
#define STR_PARAMETER_TYPE_FOLDER				"folder"

typedef vector<string> VString;

typedef struct {
	string strMinimumValue;				// minimum value the parameter can receive
	string strMaximumValue;				// maximum value the parameter can receive
} AllowedRange;

typedef struct {
	int iOrder;
	string strName;
	string strDescription;
	string strValue;
	char iType;
	string strDefaultValue;
	bool bOptional;
	VString *vAllowedValues;
	AllowedRange *allowedRange;
} Parameter;

typedef map<string,Parameter*> MParameter;

// incompatibility between two parameters
typedef struct {
	VString *vParameters;
	bool bOneRequired;			// whether one of the parameters has to be defined
} Incompatibility;

typedef vector<Incompatibility> VIncompatibility;

typedef struct {
	string strParameter;
	string strValue;	
} ParameterValue;

typedef vector<ParameterValue> VParameterValue;

// dependenty of one parameter respect to another parameter(s) and their respective value(s), taken as an OR, so the said parameter can only be defined if at least one of those parameters(values) is set
typedef struct {
	string strParameter;
	VParameterValue vParameterValue;			// set of pairs parameter-value that need to be set in order for strParameter to be set
} Dependency;

typedef vector<Dependency> VDependency;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ParameterManager {

	protected:
	
		MParameter m_mParameter;
		VIncompatibility m_vIncompatibility;
		VDependency m_vDependency;
		
		// return the list of items from a string of characters of the form "item1|item2|...|itemN"
		VString *getItemList(const char *strItems);
		
		// return the range of values from the given string
		AllowedRange *getRange(const char *strRange);
		
		// checks whether a string of characters represents a floating point number
		bool isFloat(const char *strSource);
		
		// checks whether a string of characters represents an integer
		bool isInteger(const char *strSource);
		
		// checks whether a string of characters represents an integer
		bool isBoolean(const char *strSource);
		
		// return the parameter type in string format
		const char *getStrParameterType(char iType);

	public:

		// constructor
		ParameterManager();

		// destructor
		~ParameterManager();
		
		// define a command line parameter
		bool defineParameter(const char *strName, const char *strDescription, char iType, 
			bool bOptional = false, const char *strValueConstraints = NULL, const char *strDefaultValue = NULL);
			
		// parse the given parameter-value types
		void parse(VParameterValue &vParameterValue);	
			
		// set a parameter value
		void setParameterValue(const char *strName, const char *strValue);
		
		// return a parameter value
		const char *getParameterValue(const char *strName);	
		
		// return whether a parameter value is set
		bool isParameterSet(const char *strParameter);	
		
		// print the parameters defined
		void print();
		
		// set allowed values for the given parameter
		//bool setParameterAllowedValues(const char *strParameter, const char *strAllowedValues);
		
		// set the allowed range for a parameter value 
		//bool setParameterAllowedRange(const char *strParameter, const char *strMinimumValue, const char *strMaximumValue);
		
		// define an incompatibility between parameters (if one is set, any of the others can be set)
		bool defineIncompatibility(const char *strParameters, bool bOneRequired);
		
		// define a dependency of one parameter respect to the other parameters and their respective values
		// the parameters and values are expressed as follows: parameter1(value1)|...|parameterN(valueN)
		bool defineDependency(const char *strParameter, const char *strParametersAndValues);
		
		// convert a string to a bool
		inline static bool str2bool(const char* str) {	
			if (strcmp(str,"yes") == 0) {
				return true;
			} else {
				assert(strcmp(str,"no") == 0);
				return false;
			}
		}
		
		// return a parameter in string format
		inline const char *getStrParameterValue(const char *strParameter) {
			
			MParameter::iterator it = m_mParameter.find(strParameter);
			assert(it != m_mParameter.end());
		
			return it->second->strValue.c_str();
		}
		
		// return a parameter in integer format
		inline int getIntParameterValue(const char *strParameter) {
			
			MParameter::iterator it = m_mParameter.find(strParameter);
			assert(it != m_mParameter.end());
			assert(it->second->iType == PARAMETER_TYPE_INTEGER);
		
			return atoi(it->second->strValue.c_str());
		}

		// return a parameter in floating point format
		inline float getFloatParameterValue(const char *strParameter) {
			
			MParameter::iterator it = m_mParameter.find(strParameter);
			assert(it != m_mParameter.end());
			assert(it->second->iType == PARAMETER_TYPE_FLOAT);
		
			return (float)atof(it->second->strValue.c_str());
		}
		
		// return a parameter in bool format
		inline bool getBoolParameterValue(const char *strParameter) {
		
			MParameter::iterator it = m_mParameter.find(strParameter);
			assert(it != m_mParameter.end());
			assert(it->second->iType == PARAMETER_TYPE_BOOLEAN);
			
			return str2bool(it->second->strValue.c_str());	
		}
		
		// print a dependency
		void print(Dependency *dependency);
};

};	// end-of-namespace

#endif
