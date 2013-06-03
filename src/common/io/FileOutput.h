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


#ifndef FILEOUTPUT_H
#define FILEOUTPUT_H

using namespace std;

#include <fstream>
#include <iostream>

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FileOutput {

	private:
	
		bool m_bBinary;			// whether binary or text mode
		ofstream m_os;
		string m_strFile;	

	public:

		// constructor
		FileOutput(const char *strFile, bool bBinary);

		// destructor
		~FileOutput();
		
		// open the file
		void open();
		
		// close the file
		void close();		
		
		// return the stream
		ostream &getStream() {
		
			return m_os;
		}
};

};	// end-of-namespace

#endif
