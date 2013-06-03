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


#include "TimeUtils.h"

#include <string.h>

#if defined __linux__ || defined __APPLE__ || __MINGW32__
#include <sys/time.h>
#elif _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif

namespace Bavieca {

// return the current time measured in milliseconds
double TimeUtils::getTimeMilliseconds() {

#if defined __linux__ || defined __APPLE__ || __MINGW32__
	struct timeval tv;
	struct timezone tz;
	struct tm *tm;
	gettimeofday(&tv, &tz);
	tm=localtime(&tv.tv_sec);
	//printf(" %d:%02d:%02d %d \n", tm->tm_hour, tm->tm_min, tm->tm_sec, tv.tv_usec);
	return (tv.tv_sec*1000 + tv.tv_usec/1000);
#elif _MSC_VER
	FILETIME ft;
	unsigned __int64 tmpres = 0;
    GetSystemTimeAsFileTime(&ft);
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
	return tmpres/10000.0;
	//SYSTEMTIME t;
	//GetSystemTime(&t);
	//double d = t.wMonth* +((float)t.wMilliseconds)+((float)t.wMilliseconds)+((float)t.wMilliseconds)+((float)t.wMilliseconds);  
	//return 
#endif

}

// return the date and time in a string
void TimeUtils::getDateTime(string &strDateTime) {

#if defined __linux__ || defined __APPLE__ || __MINGW32__
	char strAux[100];

	time_t timeAux;
	struct tm * tmTime;

	time(&timeAux);
	tmTime = localtime (&timeAux);

	// get the date and time in string format
	asctime_r(tmTime,strAux);	
	// remove the end of line
	while(strAux[strlen(strAux)-1] == '\n') {
		strAux[strlen(strAux)-1] = 0;
	}
	
	strDateTime = strAux;
#elif _MSC_VER
	char strAux[100];
	SYSTEMTIME t;
	//GetSystemTime(&t);    //(this is UTC time)
	GetLocalTime(&t);
	sprintf(strAux, "%04d/%02d/%02d %02d:%02d:%02d:%04d", t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond, t.wMilliseconds);
	strDateTime = strAux;
#endif
}

};	// end-of-namespace



