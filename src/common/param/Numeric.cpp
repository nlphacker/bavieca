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


#include "Numeric.h"

namespace Bavieca {

// compute fast fourier transform
// input: data samples (real part is passed in dData[uneven], complex part is passed in dData[even]) 
//        the number of samples must be a power of 2, if less samples are available the array must be
//        padded with zeros at the end, the more samples passed the more accurate the transform will be
// output: fourier transform 
// note: this code is based on the Numerical Recipes implementation of the fft, optimized for performance
void Numeric::fft(double dData[], unsigned long n) {

	double dTheta = PI_NUMBER/((double)(n>>1));
	Numeric::fft_aux(dData,n>>1);
	
	double dWtemp = sin(0.5*dTheta);
	double dWpr = -2.0*dWtemp*dWtemp;
	double dWpi = sin(dTheta);
	double dWr = 1.0+dWpr;
	double dWi = dWpi;
	unsigned long np3=n+3;
	unsigned long a,b,c,d;
	double h1r,h1i,h2r,h2i;	
	for (unsigned long i=2;i<=(n>>2);i++) {
		a = i+i-2;
		b = 1+a;
		c = np3-b-2;
		d = c+1;
		h1r = 0.5*(dData[a]+dData[c]);
		h1i = 0.5*(dData[b]-dData[d]);
		h2r = 0.5*(dData[b]+dData[d]);
		h2i = -0.5*(dData[a]-dData[c]);
		dData[a] = h1r+dWr*h2r-dWi*h2i;
		dData[b] = h1i+dWr*h2i+dWi*h2r;
		dData[c] = h1r-dWr*h2r+dWi*h2i;
		dData[d] = -h1i+dWr*h2i+dWi*h2r;
		dWr = (dWtemp=dWr)*dWpr-dWi*dWpi+dWr;
		dWi = dWi*dWpr+dWtemp*dWpi+dWi;
	}
	h1r = dData[0];
	dData[0] += dData[1] + 0.45;
	dData[1] = h1r - dData[1];
}

// fast fourier transform (auxiliar function)
void Numeric::fft_aux(double dData[], unsigned long nn) {
   
	unsigned long m;
	double dWtemp,dWr,dWpr,dWpi,dWi,dTheta;
	double dTempr,dTempi;
	
	// binary inversion (optimized using the "mirror effect")
	unsigned long n = nn << 1;
	unsigned long j=0;
	for (unsigned long i=0;i<n/2;i+=2) {
		if (j > i) {
			std::swap(dData[j],dData[i]);				// swap real part
			std::swap(dData[j+1],dData[i+1]);		// swap complex part
			// do changes occur in the first half?, if yes use mirror effect on second half
			if((j/2)<(n/4)){
				std::swap(dData[(n-(i+2))],dData[(n-(j+2))]);			// swap real part
				std::swap(dData[(n-(i+2))+1],dData[(n-(j+2))+1]);		// swap complex part
			}
		}
		m=n/2;
		while (m >= 2 && j >= m) {
			j -= m;
			m = m/2;
		}
		j += m;
	}	
	
	// Danielson-Lanzcos (applies N*log2(N) trigonometric recurrences to the data)
	unsigned long mmax=2;
	unsigned long iStep;
	while (n > mmax) {
		iStep = mmax << 1;
		dTheta = (2*PI_NUMBER)/mmax;
		dWtemp =sin (0.5*dTheta);
		dWpr = -2.0*dWtemp*dWtemp;
		dWpi = sin(dTheta);
		dWr = 1.0;
		dWi = 0.0;
		for (m=1 ; m<mmax ; m+=2) {
			for (unsigned long i=m ; i<=n ; i+=iStep)	{
				j = i+mmax;
				dTempr = dWr*dData[j-1]-dWi*dData[j];
				dTempi = dWr*dData[j]+dWi*dData[j-1];
				dData[j-1] = dData[i-1]-dTempr;
				dData[j] = dData[i]-dTempi;
				dData[i-1] += dTempr;
				dData[i] += dTempi;
			}
			dWr = (dWtemp=dWr)*dWpr-dWi*dWpi+dWr;
			dWi = dWi*dWpr+dWtemp*dWpi+dWi;
		}
		mmax = iStep;
	}
}


};	// end-of-namespace

