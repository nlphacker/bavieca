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


#ifndef FEATURETRANSFORMLIST_H
#define FEATURETRANSFORMLIST_H

#include "Transform.h"

using namespace std;

#include <string>

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FeatureTransformList {

	private:
	
		string m_strFile;
		VTransform m_vTransform; 

	public:
		
		// constructor
		FeatureTransformList(const char *strFile);

		// destructor
		~FeatureTransformList();
		
		// load the feature transforms
		bool load();	
		
		// return the transforms
		inline VTransform *getTransforms() {
		
			return &m_vTransform;
		}

};

};	// end-of-namespace

#endif
