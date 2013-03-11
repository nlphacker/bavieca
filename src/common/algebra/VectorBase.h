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


#ifndef VECTORBASE_H
#define VECTORBASE_H

#include <assert.h>
#include <string.h>
#include <iostream>
using namespace std;

#undef min
#undef max
#undef abs

namespace Bavieca {

template<typename Real>
class MatrixBase;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real>
class VectorBase {

	protected: 
	
		int m_iDim;			// dimensionality
		Real *m_rData;		// data
		
		VectorBase() {
		}
		
		VectorBase(int iDim) {
			m_iDim = iDim;
			m_rData = NULL;
		}

	public:
	
		// return the dimensionality
		int getDim() const {
		
			return m_iDim;
		}
		
		// return the data
		Real *getData() {
		
			return m_rData;
		}
		
		// returns the given element as a r-value
		Real& operator() (int i) {
		
			assert((i>=0) && (i < m_iDim));
			return m_rData[i];
		}	
				
		// returns the given element as a r-value
		Real operator() (int i) const {
		
			assert((i>=0) && (i < m_iDim));
			return m_rData[i];
		}	
		
		// invert vector elements
		void invertElements() {
			
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] = (Real)(1.0/m_rData[i]);
			}
		}
		
		// return the minimum value
		Real min() {
		
			Real rMin = m_rData[0];
			for(int i=1 ; i<m_iDim ; ++i) {
				if (m_rData[i] < rMin) {
					rMin = m_rData[i]; 
				}
			}
			
			return rMin;
		}	
				
		// return the maximum value
		Real max() {
		
			Real rMax = m_rData[0];
			for(int i=1 ; i<m_iDim ; ++i) {
				if (m_rData[i] > rMax) {
					rMax = m_rData[i]; 
				}
			}
			
			return rMax;
		}	
				
		// add another vector
		void add(VectorBase<Real> &v) {
			
			add(1.0,v);
		}
		
		// add another vector
		template<typename Real2>
		void add(VectorBase<Real2> &v) {
			
			add(1.0,v);
		}
		
		// add another vector multiplied by a constant
		void add(Real r, const VectorBase<Real> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] += r*v(i);
			}
		}
		
		// add another vector multiplied by a constant
		template<typename Real2>
		void add(Real r, const VectorBase<Real2> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] += (Real)(r*v(i));
			}
		}
		
		// add another vector (squaring elements) multiplied by a constant
		void addSquare(Real r, const VectorBase<Real> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] += r*v(i)*v(i);
			}
		}
		
		// add another vector (squaring elements) multiplied by a constant
		template<typename Real2>
		void addSquare(Real r, const VectorBase<Real2> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] += r*v(i)*v(i);
			}
		}
		
		// set elements to zero
		void zero() {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
			}
		}
		
		// scale elements by the given factor
		void scale(Real rFactor) {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= rFactor;
			}
		}
		
		// copy from memory
		void copy(Real *rData, int iDim);
		
		// copy from memory
		template<typename Real2>
		void copy(Real2 *rData, int iDim);
		
		// copy from another vector (same length and type)
		void copy(const VectorBase<Real> &v) {
		
			assert(m_iDim == v.getDim());
			memcpy(m_rData,v.m_rData,m_iDim*sizeof(Real));	
		}
		
		// copy from another vector (same length and different type)
		template<typename Real2>
		void copy(const VectorBase<Real2> &v) {
		
			assert(m_iDim == v.getDim());
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = (Real)v(i);
			}		
		}
		
		// apply the square root
		void sqrt();
		
		// multiply by a scalar
		void mul(Real r) {
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= r;
			}
		}
		
		// set to a vector multiplied by a constant
		void mul(Real r, const VectorBase<Real> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] = r*v(i);
			}
		}
		
		// set to a vector multiplied by a constant
		template<typename Real2>
		void mul(Real r, const VectorBase<Real2> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] = (Real)(r*v(i));
			}
		}		
		
		// multiply by a matrix diagonal
		void mulDiagonal(MatrixBase<Real> &m) {
			
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= m(i,i);
			}
		}
		
		// multiply vector elements
		void mulElements(VectorBase<Real> &v) {
			
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= v(i);
			}	
		}
		
		// multiply vector elements
		template<typename Real2>
		void mulElements(VectorBase<Real2> &v) {
			
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= (Real)v(i);
			}	
		}
		
		// multiply the vector elmenents by a matrix diagonal
		template<typename Real2>
		void mulDiagonal(MatrixBase<Real2> &m) {
			
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= (Real)m(i,i);
			}
		}
		
		// return scalar from multiplying the two vectors
		Real mul(const VectorBase<Real> &v) {
		
			assert(m_iDim == v.getDim());
		
			Real r = 0;
			for(int i=0 ; i < m_iDim ; ++i) {
				r += m_rData[i] * v(i);
			}
				
			return r;
		}
		
		// assign the result of multiplying a vector by a matrix
		void mul(const VectorBase<Real> &v, MatrixBase<Real> &m) {
		
			assert(this != &v);					// otherwise allocation of temporal space would be needed
			assert(m_iDim == v.m_iDim);
			assert(m_iDim == m.getCols());
			
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
				for(int j=0 ; j < m.getRows() ; ++j) {
					m_rData[i] += v(j)*m(j,i);
				}
			}
		}
		
		// assign the result of multiplying a vector by a matrix
		template<typename Real2>
		void mul(const VectorBase<Real> &v, MatrixBase<Real2> &m) {
		
			assert(this != &v);	// otherwise allocation of temporal space would be needed
			assert(m_iDim == m.getRows());
			assert(m.getCols() == v.getDim());
				
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
				for(int j=0 ; j < m.getCols() ; ++j) {
					m_rData[i] += m(i,j)*v(j);
				}
			}	
		}
		
		// assign the result of multiplying a matrix by a vector
		void mul(MatrixBase<Real> &m, const VectorBase<Real> &v) {
		
			assert(this != &v);	// otherwise allocation of temporal space would be needed
			assert(m_iDim == m.getRows());
			assert(m.getCols() == v.getDim());
				
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
				for(int j=0 ; j < m.getCols() ; ++j) {
					m_rData[i] += m(i,j)*v(j);
				}
			}
		}	
		
		// assign the result of multiplying a matrix by a vector
		template<typename Real2>
		void mul(MatrixBase<Real2> &m, const VectorBase<Real> &v) {
		
			assert(this != &v);	// otherwise allocation of temporal space would be needed
			assert(m_iDim == m.getRows());
			assert(m.getCols() == v.getDim());
				
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
				for(int j=0 ; j < m.getCols() ; ++j) {
					m_rData[i] += m(i,j)*v(j);
				}
			}	
		}	
		
		// assign the result of multiplying a matrix by a vector
		template<typename Real2>
		void mul(MatrixBase<Real> &m, const VectorBase<Real2> &v) {
		
			assert(this != &v);	// otherwise allocation of temporal space would be needed
			assert(m_iDim == m.getRows());
			assert(m.getCols() == v.getDim());
				
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
				for(int j=0 ; j < m.getCols() ; ++j) {
					m_rData[i] += m(i,j)*v(j);
				}
			}	
		}	
		
		// assign the result of multiplying a matrix by a vector
		template<typename Real2>
		void mul(MatrixBase<Real2> &m, const VectorBase<Real2> &v) {
		
			assert(this != &v);	// otherwise allocation of temporal space would be needed
			assert(m_iDim == m.getRows());
			assert(m.getCols() == v.getDim());
				
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
				for(int j=0 ; j < m.getCols() ; ++j) {
					m_rData[i] += m(i,j)*v(j);
				}
			}	
		}	
		
		// multiply vector elements
		void squareElements() {
			
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= m_rData[i];
			}	
		}		
		
		// divide the vector by another vector
		void divide(VectorBase<Real> &v);
		
		// return the sum of all elements
		Real sum() {
		
			Real r = 0.0;
			for(int i=0 ; i<m_iDim ; ++i) {
				r += m_rData[i]; 
			}
			
			return r;
		}
		
		// return the sum of all log-elements 
		Real sumLog() {
		
			Real r = 0.0;
			for(int i=0 ; i<m_iDim ; ++i) {
				r += log(m_rData[i]); 
			}
			
			return r;
		}
		
		// floor the vector elements to the given values
		void floor(VectorBase<Real> &vFloor) {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				if (m_rData[i] < vFloor(i)) {
					m_rData[i] = vFloor(i);
				}
			}
		}
		
		// floor the vector elements to the given values
		template<typename Real2>
		void floor(VectorBase<Real2> &vFloor) {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				if (m_rData[i] < vFloor(i)) {
					m_rData[i] = (Real)vFloor(i);
				}
			}
		}
		
		// shift the vector elements to the right
		void shiftRight(int iOffset);	
		
		// return whether all elements are positive
		bool isPositive() {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				if (m_rData[i] <= 0) {
					return false;
				}
			}
			
			return true;
		}
		
		// print the vector
		void print();
		
		// write the vector
		void write(ostream &os);
		
		// write the vector data
		void writeData(ostream &os);
		
		// read the vector data
		void readData(istream &is);
};

};	// end-of-namespace

#endif
