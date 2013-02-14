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


#include "ClippingEstimator.h"
#include <limits.h>

namespace Bavieca {

// return the number of clipped samples (8 bits sample)
int ClippingEstimator::estimateClipping8(const char* cSamples, int iSamples) {

	int iSamplesClipped = 0;
	int iMax = 0;
	int iMin = 0;
	
	for(int i=0 ; i<iSamples ; ++i) {
		if (cSamples[i] == CHAR_MAX) {
			++iMax;
			if (iMax == 2) {
				iSamplesClipped+=2;
			} else {
				iSamplesClipped++;
			}
			iMin = 0;	
		}
		else if (cSamples[i] == CHAR_MIN) {
			++iMin;
			if (iMin == 2) {
				iSamplesClipped+=2;
			} else {
				iSamplesClipped++;
			}	
			iMax = 0;
		} else {
			iMax = 0;
			iMin = 0; 
		}
	}

	return iSamplesClipped;
}

// return the number of clipped samples (16 bits sample)
int ClippingEstimator::estimateClipping16(const short* sSamples, int iSamples) {

	int iSamplesClipped = 0;
	int iSamplesClippedMax = 0;
	int iSamplesClippedMin = 0;
	int iMax = 0;
	int iMin = 0;
	
	for(int i=0 ; i<iSamples ; ++i) {
		if (sSamples[i] == SHRT_MAX) {
			++iMax;
			if (iMax == 2) {
				iSamplesClipped += 2;
			} else if (iMax > 2) {
				iSamplesClippedMax++;
				iSamplesClipped++;
			}
			iMin = 0;	
		}
		else if (sSamples[i] == SHRT_MIN) {
			++iMin;
			if (iMin == 2) {
				iSamplesClipped += 2;
			} else if (iMin > 2) {
				iSamplesClippedMin++;
				iSamplesClipped++;
			}	
			iMax = 0;
		} else {
			iMax = 0;
			iMin = 0; 
		}
	}
	
	return iSamplesClipped;
}

};	// end-of-namespace
