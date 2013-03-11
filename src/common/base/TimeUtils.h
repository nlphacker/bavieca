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


#ifndef TIME_H
#define TIME_H

using namespace std;

#include <string>

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class TimeUtils {

	public:
	
		// return the current time measured in milliseconds
		static double getTimeMilliseconds();
		
		// return the date and time in a string
		static void getDateTime(string &strDateTime);
		
		// convert hundreths of a second
		static void convertHundredths(double dHundredths, int &iHours, int &iMinutes, int &iSeconds) {
		
			return convertMilliseconds(dHundredths*10,iHours,iMinutes,iSeconds);
		}		
		
		// convert milliseconds
		static void convertMilliseconds(double dMilliseconds, int &iHours, int &iMinutes, int &iSeconds) {
			
			iHours = (int)(dMilliseconds/(1000*3600));
			iMinutes = (int)(dMilliseconds/(1000*60)-iHours*60);
			iSeconds = (int)(dMilliseconds/(1000)-iHours*3600-iMinutes*60);		
		}
};

};	// end-of-namespace

#endif
