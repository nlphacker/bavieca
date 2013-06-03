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


#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <assert.h>

#include "MatrixBase.h"
#include "SMatrix.h"
#include "TMatrix.h"
#include "LogMessage.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<class Real>
class Matrix : public MatrixBase<Real> {

	public:
		
		// allocate memory
		void allocate() {
		
			assert(this->m_iCols > 0);
			assert(this->m_iRows > 0);
		#if defined __AVX__ || defined __SSE__
			// round up the dimensionality to the nearest multiple	
			unsigned int iReminder = this->m_iCols%(ALIGN_BOUNDARY/sizeof(Real));
			this->m_iStride = (iReminder == 0) ? this->m_iCols : this->m_iCols + (ALIGN_BOUNDARY/sizeof(Real))-iReminder;
			assert(this->m_iStride % (ALIGN_BOUNDARY/sizeof(Real)) == 0);
			if (posix_memalign((void**)&this->m_rData,ALIGN_BOUNDARY,this->getSize()*sizeof(Real))) {
				BVC_ERROR << "memory allocation error, unable to allocate aligned memory using posix_memalign";
			}
		#else	
			this->m_iStride = this->m_iCols;
			this->m_rData = new Real[this->getSize()];
		#endif
		}
		
		// deallocate memory
		void deallocate() {
		
			assert(this->m_rData);
		#if defined __AVX__ || defined __SSE__
			free(this->m_rData);
		#else
			delete [] this->m_rData;
		#endif
			this->m_rData = NULL;
		}

	public:

		// constructor
		Matrix() : MatrixBase<Real>(0,0) {	
		}
		
		// constructor for square matrices
		Matrix(unsigned int iDim) : MatrixBase<Real>(iDim,iDim) {
			allocate();
			this->zero();
		}
		
		// constructor
		Matrix(unsigned int iRows, unsigned int iCols) : MatrixBase<Real>(iRows,iCols) {
			allocate();
			this->zero();
		}
		
		// constructor from another matrix
		Matrix(Matrix<Real> &m) : MatrixBase<Real>(m.getRows(),m.getCols()) {
			allocate();
			memcpy(this->m_rData,m.getData(),this->getSize()*sizeof(Real));
		}		
		
		// constructor from another matrix
		template<typename Real2>
		Matrix(Matrix<Real2> &m) : MatrixBase<Real>(m.getRows(),m.getCols()) {		
			allocate(); 
			for(unsigned int i=0 ; i < m.getSize() ; ++i) {
				this->m_rData[i] = (Real2)m.getData()[i];
			}
		}
		
		// constructor from another matrix
		Matrix(MatrixBase<Real> &m) : MatrixBase<Real>(m.getRows(),m.getCols()) {
			allocate();
			memcpy(this->m_rData,m.getData(),this->getSize()*sizeof(Real));
		}		
		
		// constructor from a symmetric matrix
		Matrix(SMatrix<Real> &sm) : MatrixBase<Real>(sm.getDim(),sm.getDim()) {
		
			allocate();
			
			for(unsigned int i=0 ; i < MatrixBase<Real>::m_iRows ; i++) {
				for(unsigned int j=0 ; j < i ; j++) {
					(*this)(i,j) = (*this)(j,i) = sm(i,j);
     			}
				(*this)(i,i) = sm(i,i);
			}	
		}

		// constructor from a triangular matrix
		Matrix(TMatrix<Real> &tm) : MatrixBase<Real>(tm.getDim(),tm.getDim()) {
		
			allocate();
			
			for(unsigned int i=0 ; i < MatrixBase<Real>::m_iRows ; i++) {
				for(unsigned int j=0 ; j < i ; j++) {
					(*this)(i,j) = (*this)(j,i) = tm(i,j);
     			}
				(*this)(i,i) = tm(i,i);
			}	
		}

		// destructor
		~Matrix() {
		
			if (MatrixBase<Real>::m_rData) {
				deallocate();
			}
		}
		
		// resize the matrix
		void resize(unsigned int iRows, unsigned int iCols) {
			
			assert((iRows > 0) && (iCols > 0));
		
			if ((this->m_iRows != iRows) || (this->m_iCols != iCols)) {
			
				if (this->m_rData) {
					deallocate();
				}
				
				this->m_iRows = iRows;
				this->m_iCols = iCols;
				allocate();
			}
			
			this->zero();
		}
		
		// swap the contents of two matrices
		void swap(Matrix<Real> &m) {
		
			assert((this->m_iRows == m.getRows()) && (this->m_iCols == m.getCols()));
			std::swap(this->m_rData,m.m_rData);
		}
		
		// read the matrix
		static Matrix<Real> *read(istream &is);
		
		// read the matrix data
		void readData(istream &is);
		
};

};	// end-of-namespace

#endif
