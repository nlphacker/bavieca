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


#ifndef SMATRIX_H
#define SMATRIX_H

#include "PackedMatrix.h"

#include "f2c.h"
#include "clapack.h"
#include "cblas.h"

#undef min
#undef max
#undef abs

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real>
class SMatrix : public PackedMatrix<Real> {

	public:

		// contructor
		SMatrix() : PackedMatrix<Real>() {
		}

		// contructor
		SMatrix(int iDim) : PackedMatrix<Real>(iDim) {
		}

		// copy contructor
		SMatrix(SMatrix<Real> &m) : PackedMatrix<Real>(m) {
		}

		// copy contructor
		template<typename Real2>
		SMatrix(SMatrix<Real2> &m) {
		
			this->m_iDim = m.getDim();
			this->allocate();
			this->copy(m);	
		}

		// destructor
		~SMatrix() {
		}		
		
		// return a matrix element (stored by rows as a lower triangular matrix)
		Real& operator()(int iRow, int iCol) {
		
			assert((iRow >= 0) && (iRow < this->m_iDim));
			assert((iCol >= 0) && (iCol < this->m_iDim));
			if (iRow < iCol) {
				std::swap(iRow,iCol);
			}
			
			return this->m_rData[(iRow*(iRow+1))/2+iCol];
		}	
		
		// return a matrix element (stored by rows as a lower triangular matrix)
		Real operator()(int iRow, int iCol) const {
		
			assert((iRow >= 0) && (iRow < this->m_iDim));
			assert((iCol >= 0) && (iCol < this->m_iDim));
			if (iRow < iCol) {
				std::swap(iRow,iCol);
			}
			
			return this->m_rData[(iRow*(iRow+1))/2+iCol];
		}
		
		// compute matrix inverse (returns determinant)
		float invert();	
		
		// print matrix
		void print();

};

};	// end-of-namespace

#endif
