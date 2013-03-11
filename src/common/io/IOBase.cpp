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


#include "IOBase.h"
#include "LogMessage.h"

namespace Bavieca {

// read bool value
void IOBase::read(istream &is, bool *bData, bool bBinary) {

	if (bBinary) {
		is.read(reinterpret_cast<char*>(bData),sizeof(bool));	
	} else {
		is >> *bData;	
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read unsigned char value
void IOBase::read(istream &is, char *iData, bool bBinary) {

	if (bBinary) {
		is.read(iData,sizeof(char));	
	} else {
		is >> *iData;	
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read unsigned char value
void IOBase::read(istream &is, unsigned char *iData, bool bBinary) {

	if (bBinary) {
		is.read(reinterpret_cast<char*>(iData),sizeof(unsigned char));	
	} else {
		is >> *iData;	
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read int value
void IOBase::read(istream &is, int *iData, bool bBinary) {

	if (bBinary) {
		is.read(reinterpret_cast<char*>(iData),sizeof(int));	
	} else {
		is >> *iData;
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read unsigned int value
void IOBase::read(istream &is, unsigned int *iData, bool bBinary) {

	if (bBinary) {
		is.read(reinterpret_cast<char*>(iData),sizeof(unsigned int));	
	} else {
		is >> *iData;
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read float value
void IOBase::read(istream &is, float *fData, bool bBinary) {

	if (bBinary) {
		is.read(reinterpret_cast<char*>(fData),sizeof(float));	
	} else {
		is >> *fData;	
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read double value
void IOBase::read(istream &is, double *dData, bool bBinary) {

	if (bBinary) {
		is.read(reinterpret_cast<char*>(dData),sizeof(double));	
	} else {
		is >> *dData;
	}
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read a string
void IOBase::readString(istream &is, char **str) {
	
	// read length
	int iElements = -1;
	is.read(reinterpret_cast<char*>(&iElements),sizeof(int));	
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
	
	// read the actual elements
	*str = new char[iElements+1];
	is.read(*str,sizeof(char)*iElements);	
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
	(*str)[iElements] = 0;	
}

// read a string
void IOBase::readString(istream &is, string &str) {

	is >> str;
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read bytes
void IOBase::readBytes(istream &is, char *cData, int iBytes) {

	is.read(cData,iBytes);	
	if (is.fail()) {
		BVC_ERROR << "error reading from stream at position: " << is.tellg();
	}
}

// read whitespaces (text mode)
void IOBase::readWhiteSpaces(istream &is) {
	
	char c;
	while(is.peek() == ' ') {
		is >> c; 
	}
}

// write bool value
void IOBase::write(ostream &os, bool bData, bool bBinary) {

	if (bBinary) {
		os.write(reinterpret_cast<char*>(&bData),sizeof(bool));	
	} else {
		os << bData;
	}
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write unsigned char value
void IOBase::write(ostream &os, unsigned char iData, bool bBinary) {

	if (bBinary) {
		os.write(reinterpret_cast<char*>(&iData),sizeof(unsigned char));	
	} else {
		os << iData;
	}
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write int value
void IOBase::write(ostream &os, int iData, bool bBinary) {

	if (bBinary) {
		os.write(reinterpret_cast<char*>(&iData),sizeof(int));	
	} else {
		os << iData;
	}
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write int value
void IOBase::write(ostream &os, unsigned int iData, bool bBinary) {

	if (bBinary) {	
		os.write(reinterpret_cast<char*>(&iData),sizeof(unsigned int));	
	} else {
		os << iData;
	}
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write float value
void IOBase::write(ostream &os, float fData, bool bBinary) {

	if (bBinary) {
		os.write(reinterpret_cast<char*>(&fData),sizeof(float));	
	} else {
		os << fData;
	}
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write double value
void IOBase::write(ostream &os, double dData, bool bBinary) {

	if (bBinary) {
		os.write(reinterpret_cast<char*>(&dData),sizeof(double));	
	} else {
		os << dData;
	}
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write a string
void IOBase::writeString(ostream &os, const char *str, int iLength) {

	// write length
	os.write(reinterpret_cast<char*>(&iLength),sizeof(int));	
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
	
	// write the actual elements
	os.write(str,sizeof(char)*iLength);	
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write a string
void IOBase::writeString(ostream &os, string &str) {

	os << str;
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

// write a string
void IOBase::writeString(ostream &os, ostringstream &oss) {

	string str = oss.str();
	os << str;
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}



// write double value
void IOBase::writeBytes(ostream &os, char *bData, int iBytes) {

	os.write(bData,iBytes);	
	if (os.fail()) {
		BVC_ERROR << "error writing to stream";
	}
}

};	// end-of-namespace

