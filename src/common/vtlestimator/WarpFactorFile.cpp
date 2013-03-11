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


#include "WarpFactorFile.h"

namespace Bavieca {

// constructor
WarpFactorFile::WarpFactorFile(const char *strFile) {

	m_strFile = strFile;
}

// destructor
WarpFactorFile::~WarpFactorFile() {

	for(VWarpFactorSpeaker::iterator it = m_vWarpFactorSpeaker.begin() ; it != m_vWarpFactorSpeaker.end() ; ++it) {
		delete *it;	
	}
}

// load warp factors from the file
bool WarpFactorFile::load() {

	// open the file
	FILE *file = fopen(m_strFile.c_str(),"r");
	if (file == NULL) {
		return false;
	}

	char strSpeakerId[1024+1];
	float fWarpFactor;

	// read the file line by line
	while(fscanf(file,"%s %f\n",strSpeakerId,&fWarpFactor) == 2) {
		WarpFactorSpeaker *warpFactorSpeaker = new WarpFactorSpeaker();
		warpFactorSpeaker->strSpeakerId = strSpeakerId;
		warpFactorSpeaker->fWarpFactor = fWarpFactor;
		m_vWarpFactorSpeaker.push_back(warpFactorSpeaker);
	}	
	
	// close the file
	if (fclose(file) == EOF) {
		return false;
	}

	return true;
}

};	// end-of-namespace

