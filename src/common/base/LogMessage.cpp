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

#include "LogMessage.h"

namespace Bavieca {

LogMessage::LogMessage(const char *strType, const char *strFile, const char *strFunction, int iLine)
{
	m_strType = strType;
	if ((m_strType.compare("Error") == 0) || (m_strType.compare("Warning") == 0)) {
		m_stream << strFile << " " << strFunction << " " << iLine << " ";
	}
}

LogMessage::~LogMessage()
{
	if (m_strType.compare("Error") == 0) {
		std::cerr << "Error: " << m_stream.str() << endl;
		throw std::runtime_error(m_stream.str());
	} else if (m_strType.compare("Warning") == 0) {
		std::cerr << "Warning: " << m_stream.str() << endl;
	} else if (m_strType.compare("Information") == 0) {
		std::cout << m_stream.str() << endl;
	} else if (m_strType.compare("Verbose") == 0) {
#ifdef BVC_VERBOSE_ENABLED
		std::cout << "[verb] " << m_stream.str() << endl;
#endif
	} else {
		std::cerr << m_stream.str()  << endl;
	}
}

};	// end-of-namespace


