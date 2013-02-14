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


#ifndef CLIPPINGESTIMATOR_H
#define CLIPPINGESTIMATOR_H

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class ClippingEstimator {

	public:
	
		// return the number of clipped samples
		static int estimateClipping8(const char* cSamples, int iSamples);
		
		// return the number of clipped samples
		static int estimateClipping16(const short* sSamples, int iSamples);
};

};	// end-of-namespace

#endif
