/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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

#include "LogMessage.h"

LogMessage::LogMessage(const char *strType, const char *strFile, const char *strFunction, int iLine)
{
	m_strType = strType;
	m_stream << strFile << " " << strFunction << " " << iLine << " ";
}


LogMessage::~LogMessage()
{
	if (m_strType.compare("Error") == 0) {
		std::cerr << "Error: " << m_stream.str() << "\n";
		throw std::runtime_error(m_stream.str());
	} else if (m_strType.compare("Warning") == 0) {
		std::cerr << "Warning: " << m_stream.str() << "\n";
	} else {
		std::cerr << m_stream.str()  << "\n";
	}
}


