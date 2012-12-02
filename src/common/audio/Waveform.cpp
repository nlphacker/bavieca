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


#include <math.h>
#include <stdlib.h>

#include "Waveform.h"

namespace Bavieca {

// computes the log energy of a speech frame
float Waveform::computeLogEnergy(float *fFrame, int iWindowLength){

	float fEnergy = 0.0; 
	for(int i = 0 ; i < iWindowLength ; i++) {
		fEnergy += fFrame[i] * fFrame[i];
	}

	fEnergy /= ((float)iWindowLength);
	if (fEnergy > 1) {
		fEnergy = log10(fEnergy);
	} else {
		fEnergy = 0.0;
	}
	 
   return fEnergy;
}

// applies pre emphasis to the speech signal
void Waveform::preEmphasize(float *speech, int iNumSamples, float a) {

	for(int i=iNumSamples-1 ; i>=1 ; i--) {
		speech[i] -= a*speech[i-1];
	}
	speech[0] *= 1-a;
}


// remove the DC component form the waveform
void Waveform::removeDC(short *sRawData, int iNumSamples) {

	int i;	

	// compute max input
	float fMaxInput = 0;
	for(i=0;i<iNumSamples;i++){
		if (abs(sRawData[i]) > fMaxInput) {
			fMaxInput = (float)(abs(sRawData[i]));
		}
	}
	// DC
	float fDC = 0.0;
	for(i=0;i<iNumSamples;i++) {
		fDC += (float)sRawData[i];
	}
	fDC /= (float)iNumSamples;
	// enhanced data
	float *fEnhanced = new float[iNumSamples];
	for(i=0;i<iNumSamples;i++) {
		fEnhanced[i] = (float)sRawData[i]-fDC;
	}
	// compute max output
	float fMaxOutput = 0.0;
	for(i=0;i<iNumSamples;i++) {
		if (fabs(fEnhanced[i]) > fMaxOutput) {
			fMaxOutput = fabs(fEnhanced[i]);
		}
	}
	// compute gain
	float fGain = fMaxInput/fMaxOutput;
	for(i=0;i<iNumSamples;i++) {
		fEnhanced[i] *= fGain;
	}
	convert(fEnhanced,sRawData,iNumSamples);
	delete [] fEnhanced;
}

// data type conversion
void Waveform::convert(float *fV, short *sV, int iLength) {

	for(int i = 0 ; i < iLength ; ++i) {
		if (fV[i] > ((float)SHRT_MAX)) {
			sV[i] = SHRT_MAX;
		} else if (fV[i] < ((float)SHRT_MIN)) {
			sV[i] = SHRT_MIN;
		} else {
			sV[i] = (short)fV[i];
		}
	}
}

};	// end-of-namespace
