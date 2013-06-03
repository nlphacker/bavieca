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


#ifndef LOGMESSAGE_H
#define LOGMESSAGE_H

#include "Global.h"

#include <iostream>

using namespace std;

#include <sstream>

namespace Bavieca {

#ifdef _MSC_VER
#define __func__ __FUNCTION__
#endif

#define BVC_MESSAGE(type) LogMessage(type,__func__,__FILE__,__LINE__).getStream()
#define BVC_ERROR BVC_MESSAGE("Error")								// error
#define BVC_WARNING BVC_MESSAGE("Warning")						// warning 
#define BVC_INFORMATION BVC_MESSAGE("Information")				// high-level information
#define BVC_VERB BVC_MESSAGE("Verbose")							// detailed information

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class LogMessage {

	private:
	
		ostringstream m_stream;
		string m_strType;

	public:

		LogMessage(const char *strType, const char *strFile, const char *strFunction, int iLine);

		~LogMessage();
		
		ostringstream &getStream() {
		
			return m_stream;
		}
};

};	// end-of-namespace

#endif
