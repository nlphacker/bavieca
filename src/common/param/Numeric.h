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


#ifndef NUMERIC_H
#define NUMERIC_H

using namespace std;

#include <algorithm>
#include <math.h>

namespace Bavieca {

const double PI_NUMBER = 2.0*acos(0.0);

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class Numeric {

	private:
	
		// fast fourier transform (auxiliar function)
		static void fft_aux(double dData[], unsigned long nn);	

	public:

		// addition of two log-values (natural logarithm)
		static inline float logAddition(float fLnA, float fLnB) {
		
			if (fabs(fLnA-fLnB) >= log(8388607.0)) { // 2^23-1 = 8388607
				return max(fLnA,fLnB);   // necessary to prevent dShiftedLnA = fLnA - fLnB from having a too big value
			}
			else {
				double dShiftedLnA = fLnA-fLnB; 
				double dShiftedSum = exp(dShiftedLnA)+1.0;  
				float fShiftedSumLn = (float)log(dShiftedSum); 
				float fUnshiftedSumLn = fShiftedSumLn + fLnB; 
				return fUnshiftedSumLn;
			}
		}

		// addition of two log-values (natural logarithm)
		static inline double logAddition(double dLnA, double dLnB) {
		
			if (fabs(dLnA-dLnB) >= log(8388607.0)) { // 2^23-1 = 8388607
				return max(dLnA,dLnB);   // necessary to prevent dShiftedLnA = dLnA - dLnB from having a too big value
			}
			else {
				double dShiftedLnA = dLnA-dLnB; 
				double dShiftedSum = exp(dShiftedLnA)+1.0;  
				double dShiftedSumLn = log(dShiftedSum); 
				double dUnshiftedSumLn = dShiftedSumLn + dLnB; 
				return dUnshiftedSumLn;
			}
		}
		
		// fast fourier transform
		static void fft(double dData[], unsigned long n);
		
};

};	// end-of-namespace

#endif
