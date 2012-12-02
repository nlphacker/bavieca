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


#include "FeatureTransformList.h"
#include "Global.h"

namespace Bavieca {

// constructor
FeatureTransformList::FeatureTransformList(const char *strFile)
{
	m_strFile = strFile;
}

// destructor
FeatureTransformList::~FeatureTransformList()
{
	for(VTransform::iterator it = m_vTransform.begin() ; it != m_vTransform.end() ; ++it) {
		delete *it;	
	}
}

// load the feature transforms
bool FeatureTransformList::load() {

	// open the file
	FILE *file = fopen(m_strFile.c_str(),"rt");
	if (file == NULL) {
		printf("Error: unable to open the list of transforms: %s\n",m_strFile.c_str());
		return false;
	}
	
	// load the transforms one by one
	char strFileTransform[1024+1];
	char strLine[1024+1];
	while(fgets(strLine,1024,file) != NULL) {
	
		sscanf(strLine,"%s",strFileTransform);
		Transform *transform = new Transform();
		transform->load(strFileTransform);
		m_vTransform.push_back(transform);
	}
	
	// close the file
	if (fclose(file) == EOF) {
		printf("Error: unable to close the list of transforms: %s\n",m_strFile.c_str());
		return false;
	}

	return true;
}

};	// end-of-namespace
