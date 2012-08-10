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

#include "FileInput.h"

// constructor
FileInput::FileInput(const char *strFile, bool bBinary)
{
	m_strFile = strFile;
	m_bBinary = bBinary;
}

// destructor
FileInput::~FileInput()
{
	
}

// open the file
void FileInput::open() {
	
	m_is.open(m_strFile.c_str(),ios::binary);
	if (!m_is.is_open()) {
		//EXCEPTION("unable to open the file: " << m_strFile.c_str());
		EXCEPTION("unable to open the file: ");
	}
}

// open the file
void FileInput::close() {

	m_is.close();
}

// return the file size in bytes
int FileInput::size() {

	m_is.seekg(0, ios::end);
	int iBytes = m_is.tellg();
	m_is.seekg(0, ios::beg);	
	
	return iBytes;
}




