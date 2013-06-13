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


#include "ParameterManager.h"

namespace Bavieca {

// constructor
ParameterManager::ParameterManager()
{
}

// destructor
ParameterManager::~ParameterManager()
{
	for(MParameter::iterator it = m_mParameter.begin() ; it != m_mParameter.end() ; ++it) {
		if (it->second->vAllowedValues != NULL) {
			delete it->second->vAllowedValues;
		}
		if (it->second->allowedRange != NULL) {
			delete it->second->allowedRange;
		}
		delete it->second;
	}
	for(VIncompatibility::iterator it = m_vIncompatibility.begin() ; it != m_vIncompatibility.end() ; ++it) {
		delete (*it).vParameters;
	}
}

// define a command line parameter
bool ParameterManager::defineParameter(const char *strName, const char *strDescription, char iType, 
	bool bOptional, const char *strValueConstraints, const char *strDefaultValue) {

	// check valid type
	if ((iType != PARAMETER_TYPE_STRING) && (iType != PARAMETER_TYPE_BOOLEAN) && 
		(iType != PARAMETER_TYPE_INTEGER) && (iType != PARAMETER_TYPE_FLOAT) && 
		(iType != PARAMETER_TYPE_FILE) && (iType != PARAMETER_TYPE_FOLDER)) {
		return false;
	}
	
	// check that the parameter does not already exist
	if (m_mParameter.find(strName) != m_mParameter.end()) {
		return false;
	}	
	
	// check length
	if (strlen(strName) > MAX_PARAMETER_NAME_LENGTH) {
		return false;
	}
	
	Parameter *parameter = new Parameter;
	parameter->iOrder = (int)m_mParameter.size();
	parameter->strName = strName;
	parameter->strDescription = strDescription;
	parameter->strValue = "";
	parameter->iType = iType;
	parameter->bOptional = bOptional;
	if (strDefaultValue != NULL) {
		parameter->strDefaultValue = strDefaultValue;
	} else {
		strDefaultValue = "";
	}
	parameter->vAllowedValues = NULL;
	parameter->allowedRange = NULL;	
	// there can be two types of contraints: allowed values or a range
	if (strValueConstraints != NULL) {
		// is it a range?
		if (strValueConstraints[0] == '[') {
			parameter->allowedRange = getRange(strValueConstraints);
			if (parameter->allowedRange == NULL) {
				return false;
			}
		} 
		// it is a set of allowed values
		else {
			parameter->vAllowedValues = getItemList(strValueConstraints);
			if (parameter->vAllowedValues == NULL) {
				return false;
			}
		}	
	}	
	
	m_mParameter.insert(MParameter::value_type(strName,parameter));

	return true;
}

// parse the given parameter-value types
void ParameterManager::parse(VParameterValue &vParameterValue) {

	// (1) checks
	if (m_mParameter.empty()) {
		BVC_ERROR << "no parameter definition found";
	}

	// (2) read command line parameters 
	for(VParameterValue::iterator it = vParameterValue.begin() ; it != vParameterValue.end() ; ++it) {	
		// check if it is a valid parameter
		MParameter::iterator jt = m_mParameter.find((*it).strParameter);
		if (jt == m_mParameter.end()) {
			BVC_ERROR << "unrecognized parameter: " << (*it).strParameter;
		}
		// check that the parameter is assigned a value just once
		if (jt->second->strValue.compare("") != 0) {
			BVC_ERROR << "parameter " << (*it).strParameter << " defined multiple times";
		}
		// keep the parameter value
		jt->second->strValue = (*it).strValue;
	}

	// (3) check that all the required command-line parameters receive a value
	for(MParameter::iterator it = m_mParameter.begin() ; it != m_mParameter.end() ; ++it) {
		// check value
		if (it->second->strValue.compare("") == 0) {
			if (it->second->bOptional == false) {
				BVC_ERROR << "undefined parameter " << it->first << " is required";
			} else {
				it->second->strValue = it->second->strDefaultValue;
			}
		} else {
			// check type
			switch(it->second->iType) {
				case PARAMETER_TYPE_STRING: {
					break;
				}
				case PARAMETER_TYPE_BOOLEAN: {
					if (!isBoolean(it->second->strValue.c_str())) {
						BVC_ERROR << "parameter type mismatch: \"" << it->second->strName << " = " << it->second->strValue << "\"";
					}
					break;
				}
				case PARAMETER_TYPE_INTEGER: {
					if (!isInteger(it->second->strValue.c_str())) {
						BVC_ERROR << "parameter type mismatch: \"" << it->second->strName << " = " << it->second->strValue << "\"";
					}
					break;
				}
				case PARAMETER_TYPE_FLOAT: {
					if (!isFloat(it->second->strValue.c_str())) {
						BVC_ERROR << "parameter type mismatch: \"" << it->second->strName << " = " << it->second->strValue << "\"";
					}
					break;
				}
				case PARAMETER_TYPE_FILE: {
					break;
				}
				case PARAMETER_TYPE_FOLDER: {
					break;
				}
			}
			// check allowed values
			if (it->second->vAllowedValues != NULL) {
				bool bMatch = false;
				for(VString::iterator jt = it->second->vAllowedValues->begin() ; jt != it->second->vAllowedValues->end() ; ++jt) {
					if (it->second->strValue.compare(*jt) == 0) {
						bMatch = true;
						break;
					}	
				}
				if (bMatch == false) {
					BVC_ERROR << "parameter " << it->second->strName << " receives an invalid value: \"" 
						<< it->second->strValue << "\"";
				}
			}
			// check allowed range		
			if (it->second->allowedRange != NULL) {
				if (it->second->iType == PARAMETER_TYPE_INTEGER) {
					int iValue = atoi(it->second->strValue.c_str());
					// if a minimum value is defined make the corresponding check
					if (it->second->allowedRange->strMinimumValue.compare("") != 0) {
						int iMinimumValue = atoi(it->second->allowedRange->strMinimumValue.c_str());
						if (iValue < iMinimumValue) {
							BVC_ERROR << "parameter " << it->second->strName << " receives a value out of range: " << it->second->strValue.c_str();
						}
					}
					// if a maximum value is defined make the corresponding check
					if (it->second->allowedRange->strMaximumValue.compare("") != 0) {
						int iMaximumValue = atoi(it->second->allowedRange->strMaximumValue.c_str());
						if (iValue > iMaximumValue) {
							BVC_ERROR << "parameter " << it->second->strName << " receives a value out of range: " << it->second->strValue.c_str();
						}
					}
				} else if (it->second->iType == PARAMETER_TYPE_FLOAT) {
					float fValue = (float)atof(it->second->strValue.c_str());
					// if a minimum value is defined make the corresponding check
					if (it->second->allowedRange->strMinimumValue.compare("") != 0) {
						float fMinimumValue = (float)atof(it->second->allowedRange->strMinimumValue.c_str());
						if (fValue < fMinimumValue) {
							BVC_ERROR << "parameter " << it->second->strName << " receives a value out of range";
						}
					}
					// if a minimum value is defined make the corresponding check
					if (it->second->allowedRange->strMaximumValue.compare("") != 0) {
						float fMaximumValue = (float)atof(it->second->allowedRange->strMaximumValue.c_str());
						if (fValue > fMaximumValue) {
							BVC_ERROR << "parameter " << it->second->strName << " receives a value out of range";
						}
					}
				} else {
					assert(0);
				}
			}
		}
	}
	
	// (5) check incompatibilities
	for(VIncompatibility::iterator it = m_vIncompatibility.begin() ; it != m_vIncompatibility.end() ; ++it) {
	
		// if one parameter is set, none of the others can be set
		bool bOneDefined = true;
		for(VString::iterator jt = (*it).vParameters->begin() ; jt != (*it).vParameters->end() ; ++jt) {
			if (m_mParameter[*jt]->strValue.compare("") != 0) {
				VString::iterator kt = jt;
				++kt;
				for( ; kt != (*it).vParameters->end() ; ++kt) {
					if (m_mParameter[*kt]->strValue.compare("") != 0) {
						BVC_ERROR << "parameters " << *jt << " and " << *kt << " cannot be simultaneously defined";
					}
				}
			}	
		}
		if (((*it).bOneRequired) && (bOneDefined == false)) {
			string str = "one of these parameters needs to be defined:";
			for(VString::iterator jt = (*it).vParameters->begin() ; jt != (*it).vParameters->end() ; ++jt) {
				str += " " + *jt;
			}
			BVC_ERROR << str;
		}	
	}
	
	// (6) check dependencies
	for(VDependency::iterator it = m_vDependency.begin() ; it != m_vDependency.end() ; ++it) {	
		// find the parameter
		if (m_mParameter[(*it).strParameter]->strValue.compare("") == 0) {
			continue;
		}
		bool bDependencyMet = false;
		for(VParameterValue::iterator kt = (*it).vParameterValue.begin() ; kt != (*it).vParameterValue.end() ; ++kt ) {
			if (m_mParameter[(*kt).strParameter]->strValue.compare((*kt).strValue) != 0) {
				bDependencyMet = true;
				break;
			}
		}
		if (bDependencyMet == false) {
			stringstream oss;	
			oss << "parameter: " << (*it).strParameter << " cannot be defined unless one of the following parameters is defined: " << endl;
			for(VParameterValue::iterator kt = (*it).vParameterValue.begin() ; kt != (*it).vParameterValue.end() ; ++kt ) {
				oss << (*kt).strParameter.c_str() << "(" << (*kt).strValue.c_str() << ")" << endl;
			}
			BVC_ERROR << oss;
		}
	}
}
		
// set a parameter value
void ParameterManager::setParameterValue(const char *strName, const char *strValue) {

	MParameter::iterator it = m_mParameter.find(strName);
	if (it == m_mParameter.end()) {
		BVC_ERROR << "parameter \"" << strName << "\" does not exist";
	}	
	it->second->strValue = strValue;
}

// return a parameter value
const char *ParameterManager::getParameterValue(const char *strName) {

	MParameter::iterator it = m_mParameter.find(strName);
	if (it == m_mParameter.end()) {
		return NULL;
	}	
	if (isParameterSet(strName) == false) {
		return NULL;
	}
	return it->second->strValue.c_str();
}

// return whether a parameter value is set
bool ParameterManager::isParameterSet(const char *strParameter) {

	MParameter::iterator it = m_mParameter.find(strParameter);
	if (it == m_mParameter.end()) {
		return false;
	}

	return (it->second->strValue.compare("") != 0);
}

// print parameters
void ParameterManager::print() {

	cout << " ------ parameters ---------------------------------------------------------\n";
	for(MParameter::iterator it = m_mParameter.begin() ; it != m_mParameter.end() ; ++it) {
		printf("%10s = %s\n",it->first.c_str(),it->second->strValue.c_str());
	}
	cout << " ----------------------------------------------------------------------------\n";
}

// set allowed values for the given parameter
/*bool ParameterManager::setParameterAllowedValues(const char *strParameter, const char *strAllowedValues) {

	// (1) find the parameter
	MParameter::iterator it = m_mParameter.find(strParameter);
	if (it == m_mParameter.end()) {
		return false;
	}
	it->second->vAllowedValues = getItemList(strAllowedValues);
	if ((it->second->vAllowedValues == NULL) || (it->second->vAllowedValues->size() < 2)) {
		return false;
	}

	return true;
}

// set the allowed range for a parameter value 
bool ParameterManager::setParameterAllowedRange(const char *strParameter, const char *strMinimumValue, const char *strMaximumValue) {

	// (1) find the parameter
	MParameter::iterator it = m_mParameter.find(strParameter);
	if (it == m_mParameter.end()) {
		return false;
	}
	
	// (2) make sure the parameter is of type integer or float
	if ((it->second->iType != PARAMETER_TYPE_INTEGER) && (it->second->iType != PARAMETER_TYPE_FLOAT)) {
		return false;
	}
		
	// (3) keep the allowed range
	assert(it->second->allowedRange == NULL);
	it->second->allowedRange = new AllowedRange;
	it->second->allowedRange->strMinimumValue = strMinimumValue;
	it->second->allowedRange->strMaximumValue = strMaximumValue;	

	return true;
}
*/

// return the list of items from a string of characters of the form "item1|item2|...|itemN"
VString *ParameterManager::getItemList(const char *strItems) {

	VString *vString = new VString;
	
	// parse the string of allowed values
	unsigned int iLength = (unsigned int)strlen(strItems);
	char strValue[1000];
	int iCharacters = 0;
	for(unsigned int i=0 ; i<iLength ; ++i) {
		if (isspace(strItems[i]) != 0) {
			return NULL;
		}
		if (strItems[i] == '|') {
			if (iCharacters <= 0) {
				return NULL;
			}
			strValue[iCharacters] = 0;
			vString->push_back(strValue);
			iCharacters = 0;
		} else {
			strValue[iCharacters++] = strItems[i];
		}
	}
	if (iCharacters > 0) {
		strValue[iCharacters] = 0;
		vString->push_back(strValue);
	}	

	return vString;
}

// return the range of values from the given string
AllowedRange *ParameterManager::getRange(const char *strRange) {
	
	unsigned int iLength = (unsigned int)strlen(strRange);
	if ((strRange[0] != '[') || (strRange[iLength-1] != ']')) {
		return NULL;
	}
	const char *strSeparator = strchr(strRange,'|');
	if (strSeparator == NULL) {
		return NULL;
	}
	int iLenMin = (int)(strSeparator-strRange-1);
	int iLenMax = iLength-iLenMin-3;
	if ((iLenMin < 1) || (iLenMax < 1)) {
		return NULL;
	}
	
	char *strMin = new char[iLenMin+1];
	char *strMax = new char[iLenMax+1];
	strncpy(strMin,strRange+1,iLenMin);
	strMin[iLenMin] = 0;
	strncpy(strMax,strSeparator+1,iLenMax);
	strMax[iLenMax] = 0;	
	AllowedRange *allowedRange = new AllowedRange;
	allowedRange->strMinimumValue = strMin;
	allowedRange->strMaximumValue = strMax;
	delete [] strMin;
	delete [] strMax;
	
	return allowedRange;
}

// define an incompatibility between parameters (if one is set, any of the others can be set)
bool ParameterManager::defineIncompatibility(const char *strParameters, bool bOneRequired) {

	Incompatibility incompatibility;

	// get the list of incompatible parameters
	incompatibility.vParameters = getItemList(strParameters);
	if ((incompatibility.vParameters == NULL) || (incompatibility.vParameters->size() < 2)) {
		return false;
	}

	// make sure all the parameters exist
	for(VString::iterator it = incompatibility.vParameters->begin() ; it != incompatibility.vParameters->end() ; ++it) {
		if (m_mParameter.find(*it) == m_mParameter.end()) {
			return false;
		}
	}
	
	incompatibility.bOneRequired = bOneRequired;
	m_vIncompatibility.push_back(incompatibility);

	return true;
}

// define a dependency of one parameter respect to the other parameters and their respective values
// the parameters and values are expressed as follows: parameter1(value1)|...|parameterN(valueN)
bool ParameterManager::defineDependency(const char *strParameter, const char *strParametersAndValues) {

	// (1) make sure the parameter exists
	MParameter::iterator it = m_mParameter.find(strParameter);
	if (it == m_mParameter.end()) {
		return false;
	}
	
	Dependency dependency;
	dependency.strParameter = strParameter;
	
	// (2) parse the parameter(value) items
	ParameterValue parameterValue;
	VString *vString = getItemList(strParametersAndValues);
	for(VString::iterator it = vString->begin() ; it != vString->end() ; ++it) {
		size_t iOpen = (*it).find("(");
		size_t iClose = (*it).rfind(")");	
		// just a parameter
		if (iOpen == string::npos) {
			if (iClose != string::npos) {
				return false;
			}
			parameterValue.strParameter = *it;
			parameterValue.strValue = "";
		} 
		// a parameter with its corresponding value
		else {
			if (((*it).rfind("(") != iOpen) || ((*it).find(")") != iClose) || (iClose-iOpen < 2) || (iClose+1 != (*it).length())) {
				return false;
			}
			parameterValue.strParameter = (*it).substr(0,iOpen);
			parameterValue.strValue = (*it).substr(iOpen+1,iClose-(iOpen+1));
			//printf("%s %s\n",parameterValue.strParameter.c_str(),parameterValue.strValue.c_str());
		}
		// check that the parameter actually exists
		MParameter::iterator jt = m_mParameter.find(parameterValue.strParameter);
		if (jt == m_mParameter.end()) {
			return false;
		}
		// check that the parameter is different from the dependent parameter
		if (parameterValue.strParameter.compare(strParameter) == 0) {
			return false;
		}
		dependency.vParameterValue.push_back(parameterValue);	
	}
	delete vString;

	m_vDependency.push_back(dependency);	

	return true;
}

// checks whether a string of characters represents a floating point number
bool ParameterManager::isFloat(const char *strSource) {

	errno = 0;
	char *strAux = (char*)strSource;	
	/*float fValue = */strtof(strSource,&strAux);
	if ((errno != 0) ||					// conversion failed (EINVAL, ERANGE)
		(strSource == strAux) ||		// conversion failed (no characters consumed)
		(*strAux != 0))					// conversion failed (trailing data)
		{
		return false;	
	}	
	
	return true;
}

// checks whether a string of characters represents an integer
bool ParameterManager::isInteger(const char *strSource) {

	errno = 0;
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
bool ParameterManager::isBoolean(const char *strSource) {

	if ((strcmp(strSource,"yes") == 0) || (strcmp(strSource,"no") == 0)) {
		return true;
	}
	
	return false;
}

// return the parameter type in string format
const char *ParameterManager::getStrParameterType(char iType) {

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

// print a dependency
void ParameterManager::print(Dependency *dependency) {

	printf("parameter: %s\n",dependency->strParameter.c_str());
	for(VParameterValue::iterator it = dependency->vParameterValue.begin() ; it != dependency->vParameterValue.end() ; ++it) {
		printf("-> parameter: %20s value: %20s\n",(*it).strParameter.c_str(),(*it).strValue.c_str());	
	}
}

};	// end-of-namespace
