/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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

#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <assert.h>

#include "SMatrix.h"
#include "TMatrix.h"

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real>
class Matrix {

	private:
	
		int m_iRows;
		int m_iCols;
		Real *m_rData;			
		
		// allocate memory
		void allocate() {
		
			assert(m_iCols > 0);
			assert(m_iRows > 0);
			m_rData = new Real[m_iCols*m_iRows];
		}
		
		// deallocate memory
		void deallocate() {
		
			assert(m_rData);
			delete [] m_rData;
		}

	public:

		// constructor
		Matrix() : m_iRows(0), m_iCols(0), m_rData(NULL) {	
		}

		// constructor for square matrices
		Matrix(int iDim) : m_iRows(iDim), m_iCols(iDim) {
			allocate();
			zero();
		}

		// constructor
		Matrix(int iRows, int iCols) : m_iRows(iRows), m_iCols(iCols) {
			allocate();
			zero();
		}
		
		// return a matrix element
		Real& operator()(int iRow, int iCol) {
			
			return m_rData[iRow*m_iCols+iCol];
		}

		// constructor from a symmetric matrix
		Matrix(SMatrix<Real> &sm) : m_iRows(sm.getDim()), m_iCols(sm.getDim()) {
		
			allocate();
			
			for(int i=0 ; i < m_iRows ; i++) {
				for(int j=0 ; j < i ; j++) {
					(*this)(i,j) = (*this)(j,i) = sm(i,j);
     			}
				(*this)(i,i) = sm(i,i);
			}	
		}

		// constructor from a triangular matrix
		Matrix(TMatrix<Real> &tm) : m_iRows(tm.getDim()), m_iCols(tm.getDim()) {
		
			allocate();
			
			for(int i=0 ; i < m_iRows ; i++) {
				for(int j=0 ; j < i ; j++) {
					(*this)(i,j) = (*this)(j,i) = tm(i,j);
     			}
				(*this)(i,i) = tm(i,i);
			}	
		}

		// destructor
		~Matrix() {
		
			if (m_rData) {
				deallocate();
			}
		}
		
		// return the number of rows
		int getRows() {
		
			return m_iRows;
		}
		
		// return the number of columns
		int getCols() {
		
			return m_iCols;
		}
		
		// return the data
		Real *getData() {
		
			return m_rData;
		}
		
		// resize the matrix
		void resize(int iRows, int iCols) {
			
			assert((iRows > 0) && (iCols > 0));
		
			if ((m_iRows != iRows) || (m_iCols != iCols)) {
			
				if (m_rData) {
					deallocate();
				}
				
				m_iRows = iRows;
				m_iCols = iCols;
				allocate();
			}
			
			zero();
		}
		
		// return the number of elements
		int getElements() {
		
			return m_iCols*m_iRows;
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
		
		// invert the matrix
		void invert() {
		
			int *pivot = new int[num_rows_];
  			int result = clapack_dgetrf(CblasColMajor, num_rows_, num_cols_, data_,
                              stride_, pivot);
		
		
			assert(0);
		}
		
		// invert matrix elements
		void invertElements() {
			
			int iElements = getElements();
			for(int i=0; i < iElements ; ++i) {
				m_rData[i] = 1.0/m_rData;
			}
		}
		
		// multiplication
		void add() {
		
			assert(0);
		}
		
		// singular value decomposition
		void svd(Vector<Real> vEigenvalues, Matrix<Real> mU, Matrix<Real> mV) {
		
		}
};

#endif
