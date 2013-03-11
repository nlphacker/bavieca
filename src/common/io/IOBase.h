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


#ifndef IOBASE_H
#define IOBASE_H

#include "FileInput.h"
#include "LogMessage.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class IOBase {

	private:

	public:

		// read bool value
		static void read(istream &is, bool *bData, bool bBinary = true);
		
		// read unsigned char value
		static void read(istream &is, char *iData, bool bBinary = true);
		
		// read unsigned char value
		static void read(istream &is, unsigned char *iData, bool bBinary = true);
		
		// read int value
		static void read(istream &is, int *iData, bool bBinary = true);

		// read unsigned int value
		static void read(istream &is, unsigned int *iData, bool bBinary = true);

		// read float value
		static void read(istream &is, float *fData, bool bBinary = true);
		
		// read double value
		static void read(istream &is, double *fData, bool bBinary = true);
		
		// read a string
		static void readString(istream &is, char **str);
		
		// read a string
		static void readString(istream &is, string &str);	
		
		// read bytes
		static void readBytes(istream &is, char *bData, int iBytes);
		
		// read whitespaces (text mode)
		static void readWhiteSpaces(istream &is);
		
		// -------------------------------------------------------------------------------------------
		
		// write bool value
		static void write(ostream &os, bool bData, bool bBinary = true);
	
		// write unsigned char value
		static void write(ostream &os, unsigned char iData, bool bBinary = true);
	
		// write int value
		static void write(ostream &os, int iData, bool bBinary = true);
	
		// write unsigned int value
		static void write(ostream &os, unsigned int iData, bool bBinary = true);
	
		// write float value
		static void write(ostream &os, float fData, bool bBinary = true);
	
		// write double value
		static void write(ostream &os, double dData, bool bBinary = true);
	
		// write a string
		static void writeString(ostream &os, const char *str, int iLength);
		
		// write a string
		static void writeString(ostream &os, string &str);
		
		// write a string
		static void writeString(ostream &os, ostringstream &oss);
	
		// write double value
		static void writeBytes(ostream &os, char *bData, int iBytes);
		
};

};	// end-of-namespace

#endif
