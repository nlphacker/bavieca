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

#include "ConfigurationFile.h"
#include "FileInput.h"
#include "LogMessage.h"
#include "IOBase.h"

namespace Bavieca {

// constructor
ConfigurationFile::ConfigurationFile(const char *strFile)
{
	m_strFile = strFile;
} 

// destructor
ConfigurationFile::~ConfigurationFile()
{
}

// load the configuration file
void ConfigurationFile::load() {

	try {
	
		FileInput file(m_strFile.c_str(),false);
		file.open();	
		
		string strLine;
		int iLine = 0;
		while(std::getline(file.getStream(),strLine)) {
			++iLine;
			// skip comments and blank lines
			if (isBlank(strLine) || (strLine.c_str()[0] == '#')) {
				continue;
			}		
			std::stringstream s(strLine);
			string strAttribute,strValue,strEqual;
			IOBase::readString(s,strAttribute);
			IOBase::readString(s,strEqual);
			getline(s,strValue);	
			// remove heading and traling whitespaces
			const size_t iBeg = strValue.find_first_not_of(" \n\r");
			const size_t iEnd = strValue.find_last_not_of(" \n\r");
			if ((iBeg == string::npos) || (iEnd == string::npos)) {
				BVC_ERROR << "unable to retrieve parameter name and value from the configuration file, line: " << iLine;
			}	
			strValue = strValue.substr(iBeg,iEnd-iBeg+1);     //do this within a function called IOBase::readStringWS	
			// check that the parameter is not already defined
			if (m_mParameterValue.find(strAttribute) != m_mParameterValue.end()) {
				BVC_ERROR << "parameter " << strAttribute << " defined multiple times in the configuration file, line: " << iLine;
			}
			m_mParameterValue.insert(MParameterValue::value_type(strAttribute,strValue));	
		}
		
		file.close();
	} catch (std::runtime_error) {
		BVC_ERROR << "unable to load the configuration file: " << m_strFile;
	}
}

// print the parameters
void ConfigurationFile::print() {

	for(MParameterValue::iterator it = m_mParameterValue.begin() ; it != m_mParameterValue.end() ; ++it) {
		BVC_VERB << it->first.c_str() << " -> " << it->second.c_str();
	}
}

};	// end-of-namespace
