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


#include "FileOutput.h"
#include "LogMessage.h"

namespace Bavieca {

// constructor
FileOutput::FileOutput(const char *strFile, bool bBinary)
{
	m_strFile = strFile;
	m_bBinary = bBinary;
}

// destructor
FileOutput::~FileOutput()
{
	
}

// open the file
void FileOutput::open() {
	
	if (m_bBinary) {
		m_os.open(m_strFile.c_str(),ios::binary);
	} else {
		m_os.open(m_strFile.c_str());
	}
	if (!m_os.is_open()) {
		BVC_ERROR << "unable to open the file: " << m_strFile.c_str();
	}
}

// open the file
void FileOutput::close() {

	m_os.close();
}

};	// end-of-namespace


