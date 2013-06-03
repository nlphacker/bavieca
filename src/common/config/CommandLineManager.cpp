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

#include "CommandLineManager.h"

namespace Bavieca {

CommandLineManager::CommandLineManager(const char *strApplicationName, const char *strVersion, 
	const char *strAuthor, const char *strDate) : ParameterManager()
{
	m_strApplicationName = strApplicationName;
	m_strVersion = strVersion;
	m_strAuthor = strAuthor;
	m_strDate = strDate;	
}

CommandLineManager::~CommandLineManager()
{
}

// display the application usage including parameter sintax
void CommandLineManager::displayApplicationUsage() {

	cout << endl;
	cout << " " << m_strApplicationName.c_str() << " (version: " << m_strVersion.c_str() 
		<< ", author: " << m_strAuthor.c_str() << ")" << endl;
	cout << endl;
	cout << " usage: " << m_strApplicationName.c_str() << " [parameters]" << endl;
	cout << " parameters:" << endl;
	string strType;
	
	bool bOptionalParameters = false;
	// compute max lenght of string types
	unsigned int iMaxLength = 0;
	for(MParameter::iterator it = m_mParameter.begin() ; it != m_mParameter.end() ; ++it) {
		if (it->second->bOptional) {
			bOptionalParameters = true;
		}	
		getStrType(it->second,strType);
		if (strType.length() > iMaxLength) {
			iMaxLength = (int)strType.length();
		}
	}
	
	for(unsigned int iPosition = 0 ; iPosition < m_mParameter.size() ; ++iPosition) {
		for(MParameter::iterator it = m_mParameter.begin() ; it != m_mParameter.end() ; ++it) {
			if (it->second->iOrder == (int)iPosition) {
				getStrType(it->second,strType);
				printf("  %-12s %-*s %-s",it->second->strName.c_str(),max(iMaxLength,24u),strType.c_str(),
					it->second->strDescription.c_str());
				if (it->second->strDefaultValue.compare("") == 0) {
					cout << endl;
				} else {
					cout << " (default: " << it->second->strDefaultValue.c_str() << " )" << endl;
				}
				break;
			}
		}
	}
	if (bOptionalParameters) {
		cout << endl << " (*optional)" << endl;
	}
	cout << endl;
}

// get the parameter type as a string
void CommandLineManager::getStrType(Parameter *parameter, string &strType) {

	strType = "[";
	if (parameter->vAllowedValues == NULL) {
		strType += getStrParameterType(parameter->iType);
	} else {
		int iItems = 0;
		for(VString::iterator jt = parameter->vAllowedValues->begin() ; jt != parameter->vAllowedValues->end() ; ++jt, ++iItems) {
			if (iItems > 0) {
				strType += "|";
			}
			strType += *jt;
		}
	}
	if (parameter->bOptional) {
		strType += "]*";
	} else {
		strType += "]";
	}
}
		
// parse parameters read from the command line
bool CommandLineManager::parseParameters(int argc, char *argv[]) {

	// (1) show usage if necessary
	if (argc == 1) {
		displayApplicationUsage();
		return false;
	}
	
	// (2) read command-line parameters
	VParameterValue vParameterValue;
	for(int i=1 ; i < argc ; i+=2) {
	
		ParameterValue parameterValue;
		parameterValue.strParameter = argv[i];
		if (i+1 >= argc) {
			BVC_ERROR << "parameter \"" << argv[i] << "\" has no value";
		}	
		parameterValue.strValue = argv[i+1];
		vParameterValue.push_back(parameterValue);
	}
	
	try {
		parse(vParameterValue);
	} catch (std::runtime_error) {
		return false;
	}
	
	return true;
}

};	// end-of-namespace
