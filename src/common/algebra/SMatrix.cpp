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


#include "SMatrix.h"

#include <math.h>
#include <stdlib.h>

namespace Bavieca {

// compute matrix inverse (returns determinant)
template<>
float SMatrix<float>::invert() {

	assert(0);
	printf("code commented!!!");
	exit(-1);
	/*integer *iPivot = new integer[m_iDim];
	integer iInfo;
	integer iRows = m_iDim;

	// decomposition
	ssptrf_(const_cast<char*>("U"),&iRows,m_rData,iPivot,&iInfo);
	assert(iInfo >= 0);
	
	// matrix is singular (determinant == 0), there exist no inverse
	if (iInfo > 0) {
		delete [] iPivot;
		return 0.0;
	}
	
	// to understand this see LAPACK documentation: http://www.netlib.org/lapack/	
	float fDet = 1.0;
	for (int i = 0; i < this->m_iDim ; i++) {
		if (iPivot[i] > 0) {
			fDet *= (*this)(i, i);
		} 
		// 2x2 block
		else {  
			i++;
			fDet *= (*this)(i,i)*(*this)(i-1,i-1) - (*this)(i,i-1)*(*this)(i,i-1);		// block determinant
		}
	} 
	
	// computes the inverse
	float *fPWork = new float[m_iDim];
	ssptri_(const_cast<char *>("U"), &iRows,m_rData,iPivot,fPWork,&iInfo);
	assert(iInfo >= 0);
	
	delete [] fPWork;
	delete [] iPivot;
	
	// matrix is singular (determinant == 0), there exist no inverse
	if (iInfo != 0) {	
		return 0.0;
	}

	return fDet;*/
	return 0.0;
}

// print the matrix
template<typename Real>
void SMatrix<Real>::print() {

	cout << "dim: " << this->m_iDim << endl;
	for(int i=0 ; i < this->m_iDim ; ++i) {
		for(int j=0 ; j < this->m_iDim ; ++j) {
			cout << (*this)(i,j) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

template void SMatrix<float>::print();
template void SMatrix<double>::print();

};	// end-of-namespace

