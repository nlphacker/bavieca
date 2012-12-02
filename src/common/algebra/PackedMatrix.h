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


#ifndef PACKEDMATRIX_H
#define PACKEDMATRIX_H

#include <stdio.h>

using namespace std;

#include "Vector.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real> 
class PackedMatrix {

	protected:
	
		int m_iDim;			// dimensionality of the square matrix
		Real *m_rData;		// data			
		
		PackedMatrix() {}
		
	public:

		// constructor		
		PackedMatrix(int iDim) : m_iDim(iDim), m_rData(NULL) {
			allocate();
			zero();
		}
			
		// copy constructor
		PackedMatrix(PackedMatrix<Real> &m) {
		
			m_iDim = m.m_iDim;
			allocate();
			copy(m);
		}

		// copy constructor
		template<typename Real2>
		PackedMatrix(PackedMatrix<Real2> &m) {
		
			m_iDim = m.m_iDim;
			allocate();
			copy(m);
		}

		// destructor
		~PackedMatrix() {
			if (m_rData) {
				deallocate();
			}
		}
		
		// copy data from another matrix
		void copy(const PackedMatrix<Real> &m) {
			
			int iElements = getElements();
			for(int i=0; i < iElements ; ++i) {
				m_rData[i] = m.m_rData[i];
			}
		}
		
		// copy data from another matrix
		template<typename Real2>
		void copy(const PackedMatrix<Real2> &m) {
			
			int iElements = getElements();
			for(int i=0; i < iElements ; ++i) {
				m_rData[i] = (Real)(m.getData()[i]);
			}
		}
		
		// return the packed size
		int getElements() {
		
			return (m_iDim*(m_iDim+1))/2;
		}
		
		// return the dimensionality
		int getDim() const {
		
			return m_iDim;
		}
		
		// return the data
		Real *getData() const {
		
			return m_rData;
		}
		
		// allocate memory
		void allocate() {
		
			assert(m_iDim > 0);
			m_rData = new Real[getElements()];
		}
		
		// deallocate memory
		void deallocate() {
		
			assert(m_rData);
			delete [] m_rData;
			m_rData = NULL;	
		}		
		
		// resize the matrix (also used to allocate space for the first time)
		void resize(int iDim) {
		
			if (iDim != m_iDim) {
			
				if (m_rData) {
					deallocate();
				}
				
				m_iDim = iDim;
				allocate();
			}	
			
			resize(m_iDim);
		}
		
		// return a matrix element (stored by rows as a lower triangular matrix)
		Real& operator()(int iRow, int iCol) {
		
			assert((iRow >= 0) && (iRow < m_iDim));
			assert((iCol >= 0) && (iCol < m_iDim));
			assert(iRow >= iCol);
			
			return m_rData[(iRow*(iRow+1))/2+iCol];
		}	
		
		// return a matrix element (stored by rows as a lower triangular matrix)
		Real operator()(int iRow, int iCol) const {
		
			assert((iRow >= 0) && (iRow < m_iDim));
			assert((iCol >= 0) && (iCol < m_iDim));
			assert(iRow >= iCol);
			
			return m_rData[(iRow*(iRow+1))/2+iCol];
		}	
		
		// set to zero
		void zero() {
		
			int iElements = getElements();
			for(int i=0; i < iElements ; ++i) {
				m_rData[i] = 0.0;
			}
		}
		
		// check wether all elements are zero
		bool isZero(Real rEpsilon = 0.00001) {
		
			int iElements = getElements();
			for(int i=0; i < iElements ; ++i) {
				if (fabs(m_rData[i]) > rEpsilon) {
					return false;
				}
			}
			
			return true;
		}
		
		// get the matrix diagonal
		void getDiagonal(VectorBase<Real> &vDiagonal) {
		
			assert(vDiagonal.getDim() == m_iDim);
		
			for(int i=0 ; i < vDiagonal.getDim() ; ++i) {
				vDiagonal(i) = (*this)(i,i);
			}
		}		
		
		// add two packed matrices
		void add(const PackedMatrix<Real> &m) {
		
			assert(this->m_iDim == m.m_iDim);
		
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] += m.m_rData[i];
			}
		}
		
		// add two packed matrices
		template<typename Real2>
		void add(const PackedMatrix<Real2> &m) {
		
			assert(m_iDim == m.getDim());
		
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] += m.m_rData[i];
			}
		}
		
		// add two packed matrices
		void add(Real rAlpha, const PackedMatrix<Real> &m) {
		
			assert(m_iDim == m.getDim());
			
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] += rAlpha*m.getData()[i];
			}
		}
	
		// add two packed matrices
		template<typename Real2>
		void add(Real rAlpha, const PackedMatrix<Real2> &m) {
		
			assert(m_iDim == m.getDim());
			
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] += (Real)(rAlpha*m.getData()[i]);		
			}
		}
	
		// scale elements by the given factor
		void scale(Real rAlpha) {
		
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] *= rAlpha;
			}
		}
		
		// adds the given squared vector
		void addSquare(Real rAlpha, const VectorBase<Real> &v) {
		
			assert(m_iDim == v.getDim());
			
			Real *r = m_rData;
			for(int i=0 ; i<m_iDim ; ++i) {
				for(int j=0 ; j<=i ; ++j, ++r) {
					*r += rAlpha*v(i)*v(j); 
				}	
			}
		}
		
		// adds the given squared vector
		template<typename Real2>
		void addSquare(Real rAlpha, const VectorBase<Real2> &v) {
		
			assert(m_iDim == v.getDim());
			
			Real *r = m_rData;
			for(int i=0 ; i<m_iDim ; ++i) {
				for(int j=0 ; j<=i ; ++j, ++r) {
					*r += rAlpha*v(i)*v(j); 
				}	
			}
		}		
		
		// multipl by the given constant
		void mul(Real r) {
		
			for(int i=0 ; i < getElements() ; ++i) {
				m_rData[i] *= r;
			}
		}	
		
		// add two packed matrices
		void mul(Real rAlpha, const PackedMatrix<Real> &m) {
		
			assert(m_iDim == m.getDim());
			
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] = rAlpha*m.m_rData[i];		
			}
		}
		
		// add two packed matrices
		template<typename Real2>
		void mul(Real rAlpha, const PackedMatrix<Real2> &m) {
		
			assert(m_iDim == m.getDim());
			
			int iElements = getElements();
			for(int i=0 ; i < iElements ; ++i) {
				m_rData[i] = (Real)(rAlpha*m.getData()[i]);		
			}
		}
		
		// write the matrix
		void write(ostream &os);
		
		// write the matrix data
		void writeData(ostream &os);

		// read the matrix data
		void readData(istream &is);

};

};	// end-of-namespace

#endif
