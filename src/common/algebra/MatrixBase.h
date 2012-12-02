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


#ifndef MATRIXBASE_H
#define MATRIXBASE_H

using namespace std;

// blas/lapack
#include "cblas.h"
//#include "f2c.h"
//#include "clapack.h"

#undef min
#undef max
#undef abs

#include "VectorBase.h"
#include "VectorStatic.h"
#include "Global.h"

#include <string.h>

namespace Bavieca {

// whether the Matrix is transposed
typedef enum {
	yes = CblasTrans,
	no = CblasNoTrans
} MatrixTransposed;

#define MATRIX_TRANSPOSED				'T'
#define MATRIX_NO_TRANSPOSED			'N'

template<typename Real>
class MatrixStatic;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real>
class MatrixBase {

	protected:
	
		int m_iRows;				// number of rows in the matrix
		int m_iCols;				// number of columns in the matrix
		int m_iStride;				// number of elements per row (stride >= # columns)
		Real *m_rData;				// matrix data
	
		MatrixBase(int iRows, int iCols) {
		
			m_iRows = iRows;
			m_iCols = iCols;
			m_rData = NULL;
		}
		
	public:

		// return the number of rows
		int getRows() {
		
			return m_iRows;
		}
		
		// return the number of columns
		int getCols() {
		
			return m_iCols;
		}
		
		// return the stride
		int getStride() {
		
			return m_iStride;
		}
		
		// return the data
		Real *getData() {
		
			return m_rData;
		}
		
		// return a matrix element as a r-value
		Real& operator()(int iRow, int iCol) {
			
			return m_rData[iRow*m_iStride+iCol];
		}
		
		// return a matrix element as a l-value
		Real operator()(int iRow, int iCol) const {
			
			return m_rData[iRow*m_iStride+iCol];
		}
		
		// return the number of elements in the matrix
		int getElements() {
		
			return m_iCols*m_iRows;
		}
		
		// return the matrix size in terms of data allocated
		int getSize() {
		
			return m_iStride*m_iRows;
		}	
		
		// return a matrix row
		VectorStatic<Real> getRow(int iRow);
		
		// return a matrix row
		Real *getRowData(int iRow);
		
		// copy data from another matrix
		void copy(MatrixBase<Real> &m) {
		
			assert((m_iRows == m.getRows()) && (m_iCols == m.getCols()));
			memcpy(m_rData,m.getData(),getSize()*sizeof(Real));
		}
		
		// copy a vector to the given row
		void copyRow(VectorBase<Real> &v, int iRow);
		
		// copy a vector to the given column
		void copyCol(VectorBase<Real> &v, int iColumn);		
		
		// copy the given vector to each row
		void copyRows(VectorBase<Real> &v);
		
		// copy the given vector to each column
		void copyCols(VectorBase<Real> &v);	
		
		// return a submatrix of the given matrix
		MatrixStatic<Real> getSubMatrix(int iRowStart, int iRows, int iColStart, int iCols) const;
		
		// set to zero
		void zero();
		
		// check wether all elements are zero
		bool isZero(Real rEpsilon);
		
		// test two matrices for equality
		bool equal(MatrixBase<Real> &m, Real rEpsilon = 0.0);
		
		void setDiagonal(Real r) {
		
			assert(isSquare());
			for(int i=0; i<m_iCols; ++i) {
				(*this)(i,i) = r;
			}
		}
		
		// set the matrix to the identity matrix
		void setIdentity() {
		
			assert(isSquare());
			zero();
			setDiagonal(1.0);
		}
		
		// get the matrix diagonal
		void getDiagonal(VectorBase<Real> &vDiagonal) {
		
			assert(isSquare());
			assert(vDiagonal.getDim() == m_iRows);
		
			for(int i=0 ; i < vDiagonal.getDim() ; ++i) {
				vDiagonal(i) = (*this)(i,i);
			}
		}
		
		// add a matrix
		void add(Real r, MatrixBase<Real> &matrix);
		
		// add a vector multiplication: *this += r * vec1 * vec2^T 
		void addVecMul(Real r, VectorBase<Real> &v1, VectorBase<Real> &v2);	
		
		// return whether the matrix is square
		bool isSquare() {
		
			return (m_iRows == m_iCols);
		}
		
		// multiply each element by the given constant
		void mul(Real r);	
		
		// multiply each row by the given vector
		void mulRows(VectorBase<Real> &v);	
		
		// multiply each column by the given vector
		void mulCols(VectorBase<Real> &v);	
		
		// add the product of two matrices
		void addMul(MatrixBase<Real> &m1, MatrixTransposed trans1, 
			MatrixBase<Real> &m2, MatrixTransposed trans2);
			
		// add the product of three matrices
		void addMul(MatrixBase<Real> &mA, MatrixTransposed transA, 
			MatrixBase<Real> &mB, MatrixTransposed transB, 
			MatrixBase<Real> &mC, MatrixTransposed transC);
		
		// transpose the matrix 
		void transpose() {
		
			for(int i=0;i<m_iRows;++i) {
				for(int j=0;j<i;++j) {
					Real r = (*this)(i,j);
					(*this)(i,j) = (*this)(j,i);
					(*this)(j,i) = r;
				}
			}	
		}
		
		// compute the matrix of cofactors of the given matrix
		void matrixCofactor();
		
		// compute matrix inverse
		Real invert();
	
		// invert matrix elements
		void invertElements();
		
		// multiplication
		void add();
		
		// singular value decomposition
		void svd(VectorBase<Real> vEigenvalues, MatrixBase<Real> mU, MatrixBase<Real> mV);
		
		// print the matrix data
		void print();	

		// write the matrix data
		void write(ostream &os);
		
		// write the matrix data
		void writeData(ostream &os);	
};

};	// end-of-namespace

#endif
