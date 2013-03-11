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


#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <limits.h>

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class Waveform {

	private:
	
		// data type conversion
		static void convert(float *fV, short *sV, int iLength);
		
	public:
		
		static void frequencyWarp(float *fft_mag, int length, float warp_factor);
		
		// computes the log energy of a speech frame
		static float computeLogEnergy(float *fFrame, int iWindowLength);
		
		// applies pre emphasis to the speech signal
		static void preEmphasize(float *speech, int numsamples, float a);
		
		// remove the DC component form the waveform
		static void removeDC(short *rawdata, int numsamples);	

};

};	// end-of-namespace

#endif
