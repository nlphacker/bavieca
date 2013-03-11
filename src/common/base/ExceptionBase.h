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


#ifndef EXCEPTIONBASE_H
#define EXCEPTIONBASE_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Bavieca {

// macro to throw exceptions
#define EXCEPTION(str) throw ExceptionBase(str,__FILE__,__func__,__LINE__);

//#define WARNING(str) std::cerr << "Warning: " << str << '\n'

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ExceptionBase : public std::runtime_error {

	private:
	
		std::string m_strMessage;

	public:
		
		// constructor
		ExceptionBase(const std::string &strArg, const char *strFile, const char *strFunction, int iLine) : 
			std::runtime_error(strArg) {
			
			std::ostringstream o;
			o << strFile << ":" << iLine << ": " << strArg;
			m_strMessage = o.str();
		}

		// destructor
		~ExceptionBase() throw() {	
		}
		
		// return the exception message
		const char *what() const throw() {
			return m_strMessage.c_str();
		}
};

};	// end-of-namespace

#endif
