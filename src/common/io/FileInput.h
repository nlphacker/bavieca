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


#ifndef FILEINPUT_H
#define FILEINPUT_H

using namespace std;

#include <fstream>
#include <iostream>

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FileInput {

	private:
	
		bool m_bBinary;			// whether binary or text mode
		ifstream m_is;
		string m_strFile;

	public:

		// constructor
		FileInput(const char *strFile, bool bBinary);

		// destructor
		~FileInput();
		
		// open the file
		void open();
		
		// close the file
		void close();
		
		// return the file size in bytes
		long size();
		
		// return the stream
		istream &getStream() {
		
			return m_is;
		}
};

};	// end-of-namespace

#endif
