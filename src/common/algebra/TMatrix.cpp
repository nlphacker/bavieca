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


#include "LogMessage.h"
#include "TMatrix.h"

#include <math.h>

namespace Bavieca {

// Cholesky decomposition of a symmetric matrix
template<typename Real>
void TMatrix<Real>::choleskyDecomposition(SMatrix<Real> &m) {

	Real* rS = m.getData(); 
	Real* rT = this->m_rData; 
	Real* ti = this->m_rData;
	for (int i=0; i< (int)this->m_iDim; i++) {
		Real* tj = rT; 
		Real sum; 
		int k;
		for (int j=0; j<i; j++) {
			Real* tk = ti; 
			sum = 0.0; 
			k = j;
			while (k--) { 
				sum += *tj++ * *tk++; 
			}
			*tk = (*rS++ - sum)/ *tj++;
		}
		sum = 0.0; k = i;
		while (k--) { 
			sum += (*ti)*(*ti);
			ti++;
		}
		Real rD = *rS++ - sum;
		if (rD<=0.0)  {
			BVC_ERROR << "unable to perform Cholesky decomposition, matrix is not positive definite";
		}
		*ti++ = sqrt(rD);
	}
}

// print the matrix
template<typename Real>
void TMatrix<Real>::print() {

	cout << "dim: " << this->m_iDim << endl;
	for(unsigned int i=0 ; i < this->m_iDim ; ++i) {
		for(unsigned int j=0 ; j < this->m_iDim ; ++j) {
			if (i >= j) {
				cout << (*this)(i,j) << " ";
			} else {
				cout << "0.0000000";
			}
		}
		cout << endl;
	}
	cout << endl;
}

template void TMatrix<float>::choleskyDecomposition(SMatrix<float> &m);
template void TMatrix<float>::print();
template void TMatrix<double>::print();

};	// end-of-namespace


