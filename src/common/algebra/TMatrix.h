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


#ifndef TMATRIX_H
#define TMATRIX_H

#include "PackedMatrix.h"
#include "SMatrix.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real> 
class TMatrix : public PackedMatrix<Real> {

	public:

		// contructor
		TMatrix() : PackedMatrix<Real>() {
		}

		// contructor
		TMatrix(int iDim) : PackedMatrix<Real>(iDim) {
		}

		// copy contructor
		TMatrix(TMatrix<Real> &m) : PackedMatrix<Real>(m) {
		}

		// destructor
		~TMatrix() {
		}		
		
		// Cholesky decomposition of a symmetric matrix
		void choleskyDecomposition(SMatrix<Real> &m);
		
		// matrix determinant
		double determinant() {
		
			double dDet = (*this)(0,0);
			for(int i=1 ; i < this->m_iDim ; ++i) {
				dDet *= (*this)(i,i);
			}
			
			return dDet;
		}
		
		// print the matrix
		void print();
};

};	// end-of-namespace

#endif
