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


#include "MatrixBase.h"
#include "IOBase.h"
#include "Matrix.h"

#include <stdlib.h>

namespace Bavieca {

// copy a vector to a given row
template<typename Real>
void MatrixBase<Real>::copyRow(VectorBase<Real> &v, int iRow) {

	assert(m_iCols == v.getDim());
	memcpy(m_rData+iRow*m_iCols,v.getData(),m_iCols*sizeof(Real));
}

// copy a vector to a given column
template<typename Real>
void MatrixBase<Real>::copyCol(VectorBase<Real> &v, int iCol) {

	assert(m_iRows == v.getDim());
	for(int i=0; i < m_iRows ; ++i) {
		(*this)(i,iCol) = v(i);
	}
}

// copy the given vector to each row
template<typename Real>
void MatrixBase<Real>::copyRows(VectorBase<Real> &v) {

	assert(m_iCols == v.getDim());
	for(int i=0; i < m_iRows ; ++i) {
		memcpy(m_rData+i*m_iCols,v.getData(),m_iCols*sizeof(Real));
	}
}

// copy the given vector to each column
template<typename Real>
void MatrixBase<Real>::copyCols(VectorBase<Real> &v) {

	assert(m_iRows == v.getDim());
	for(int j=0; j < m_iRows ; ++j) {
		for(int i=0; i < m_iCols ; ++i) {
			(*this)(j,i) = v(j);
		}
	}
}

// return a matrix row
template<typename Real>
VectorStatic<Real> MatrixBase<Real>::getRow(int iRow) {

	return VectorStatic<Real>(m_rData+iRow*m_iCols,m_iCols);
}

// return a matrix row
template<typename Real>
Real *MatrixBase<Real>::getRowData(int iRow) {

	return m_rData+iRow*m_iCols;
}

// return a submatrix of the given matrix
template<typename Real>
MatrixStatic<Real> MatrixBase<Real>::getSubMatrix(int iRowStart, int iRows, int iColStart, int iCols) const {

	return MatrixStatic<Real>(*this,iRowStart,iRows,iColStart,iCols);
}

// set to zero
// needs to be done row by row because of the stride
template<typename Real>
void MatrixBase<Real>::zero() {

	for(int i=0; i < m_iRows; ++i) {
		for(int j=0; j < m_iCols ; ++j) {
			m_rData[i*m_iStride+j] = 0.0;
		}
	}
}

// check wether all elements are zero
template<typename Real>
bool MatrixBase<Real>::isZero(Real rEpsilon = 0.00001) {

	for(int i=0; i < m_iRows; ++i) {
		for(int j=0; j < m_iCols ; ++j) {
			if (fabs(m_rData[i*m_iStride+j]) > rEpsilon) {
				return false;
			}
		}
	}
	
	return true;
}

// test two matrices for equality
template<typename Real>
bool MatrixBase<Real>::equal(MatrixBase<Real> &m, Real rEpsilon) {

	for(int i=0; i < m_iRows; ++i) {
		for(int j=0; j < m_iCols ; ++j) {
			if (fabs((*this)(i,j)-m(i,j)) > rEpsilon) {
				return false;
			}
		}
	}

	return true;
}

// add a matrix
template<typename Real>
void MatrixBase<Real>::add(Real r, MatrixBase<Real> &matrix) {

	assert(m_iRows == matrix.getRows());
	assert(m_iCols == matrix.getCols());
	for(int i=0; i < m_iRows ; ++i) {
		for(int j=0; j < m_iCols ; ++j) {
			(*this)(i,j) += r*matrix(i,j);
		}
	}
}

// add a vector multiplication: *this += r * vec1 * vec2^T 
template<typename Real>
void MatrixBase<Real>::addVecMul(Real r, VectorBase<Real> &v1, VectorBase<Real> &v2) {

	assert(isSquare());	
	assert(m_iRows == v1.getDim());
	assert(m_iRows == v2.getDim());	
	for(int i=0; i < m_iRows ; ++i) {
		for(int j=0; j < m_iCols ; ++j) {
			(*this)(i,j) += r*v1(i)*v2(j);
		}
	}
}


// multiply each element by the given constant
template<typename Real>
void MatrixBase<Real>::mul(Real r) {

	for(int i=0; i < m_iRows ; ++i) {
		for(int j=0; j < m_iCols ; ++j) {
			(*this)(i,j) *= r;
		}
	}
}

// multiply each row by the given vector
template<typename Real>
void MatrixBase<Real>::mulRows(VectorBase<Real> &v) {

	assert(v.getDim() == m_iCols);
	for(int i=0 ; i < m_iRows ; ++i) {
		for(int j=0 ; j < m_iCols ; ++j) {
			(*this)(i,j) *= v(j);
		}	
	}
}

// multiply each column by the given vector
template<typename Real>
void MatrixBase<Real>::mulCols(VectorBase<Real> &v) {

	assert(v.getDim() == m_iCols);
	for(int i=0 ; i < m_iRows ; ++i) {
		for(int j=0 ; j < m_iCols ; ++j) {
			(*this)(i,j) *= v(i);
		}	
	}
}

// add the product of two matrices
template<>
void MatrixBase<float>::addMul(MatrixBase<float> &m1, MatrixTransposed trans1, 
	MatrixBase<float> &m2, MatrixTransposed trans2) {
	
	int iM1Cols = (trans1 == no) ? m1.getCols() : m1.getRows();
	int iStride1 = m1.getStride();
	int iStride2 = m2.getStride();
	int iStrideSelf = m_iStride;
	
	cblas_sgemm(CblasRowMajor,static_cast<CBLAS_TRANSPOSE>(trans1),
				static_cast<CBLAS_TRANSPOSE>(trans2),
            getRows(),getCols(),iM1Cols,
            1.0,m1.m_rData,iStride1,
            m2.m_rData,iStride2,
            0.0,m_rData,iStrideSelf);
	
	// the code below only works for "column major" storage format
	/*integer iRows = getRows();
	integer iCols = getCols();
	integer iM1Cols = (trans1 == no) ? m1.getCols() : m1.getRows();
	integer iStride1 = m1.getStride();
	integer iStride2 = m2.getStride();
	integer iStrideSelf = m_iStride;
	
	float fOne = 1.0;
	float fZero = 0;
	char cTrans1 = (trans1 == yes) ? MATRIX_TRANSPOSED : MATRIX_NO_TRANSPOSED;
	char cTrans2 = (trans2 == yes) ? MATRIX_TRANSPOSED : MATRIX_NO_TRANSPOSED;
	
	sgemm_(&cTrans1,&cTrans2,&iRows,&iCols,&iM1Cols,&fOne,m1.m_rData,&iStride1,
		m2.m_rData,&iStride2,&fZero,m_rData,&iStrideSelf);*/
}

// add the product of two matrices
template<>
void MatrixBase<double>::addMul(MatrixBase<double> &m1, MatrixTransposed trans1, 
	MatrixBase<double> &m2, MatrixTransposed trans2) {
	
	int iM1Cols = (trans1 == no) ? m1.getCols() : m1.getRows();
	int iStride1 = m1.getStride();
	int iStride2 = m2.getStride();
	int iStrideSelf = m_iStride;
	
	cblas_dgemm(CblasRowMajor,static_cast<CBLAS_TRANSPOSE>(trans1),
				static_cast<CBLAS_TRANSPOSE>(trans2),
            getRows(),getCols(),iM1Cols,
            1.0,m1.m_rData,iStride1,
            m2.m_rData,iStride2,
            0.0,m_rData,iStrideSelf);
}

// add the product of three matrices
template<typename Real>
void MatrixBase<Real>::addMul(MatrixBase<Real> &mA, MatrixTransposed transA, 
	MatrixBase<Real> &mB, MatrixTransposed transB, 
	MatrixBase<Real> &mC, MatrixTransposed transC) {
	
	int iARows = (transA == no) ? mA.getRows() : mA.getCols();
	int iBCols = (transB == no) ? mB.getCols() : mB.getRows();
	int iCCols = (transC == no) ? mC.getCols() : mC.getRows();
	
	assert(m_iRows == iARows);
	assert(m_iCols == iCCols);

	Matrix<Real> mAB(iARows,iBCols);
	mAB.addMul(mA,transA,mB,transB);
	addMul(mAB,no,mC,transC);
}

// invert the matrix
template <>
float MatrixBase<float>::invert() {

	assert(isSquare());
	
	integer iLDA = m_iRows;
	integer *iPivot = new integer[m_iRows];
	integer iInfo;
	integer iRows = m_iRows;
	integer iColumns = m_iCols;
	
	// LU decomposition (factorizes a matrix as the product of a lower triangular matrix and an upper triangular matrix)
	sgetrf_(&iRows,&iColumns,m_rData,&iLDA,iPivot,&iInfo);	
	
	// matrix is singular (determinant = 0), there exist no inverse
	if (iInfo > 0) {
		return 0.0;
	} 
	
	bool bPositive = true;
	float fDet = 1.0;
	for(int i=0 ; i<m_iCols ; ++i) {
		fDet *= (*this)(i,i);
		if (iPivot[i] != (i+1)) {
			bPositive = !bPositive;
		}
	}
	
	// inverse computation based on the LU decomposition
	integer iLWork = iRows;
	float *fPWork = new float[iLWork];
	sgetri_(&iRows,m_rData,&iLDA,iPivot,fPWork,&iLWork,&iInfo);	
	
	delete [] iPivot;
	delete [] fPWork;
	
	// return the determinant
	return bPositive ? fDet : -fDet;
}

// invert the matrix
template <>
double MatrixBase<double>::invert() {

	assert(isSquare());
	
	integer iLDA = m_iRows;
	integer *iPivot = new integer[m_iRows];
	integer iInfo;
	integer iRows = m_iRows;
	integer iColumns = m_iCols;
	
	// LU decomposition (factorizes a matrix as the product of a lower triangular matrix and an upper triangular matrix)
	dgetrf_(&iRows,&iColumns,m_rData,&iLDA,iPivot,&iInfo);	
	
	// matrix is singular (determinant = 0), there exist no inverse
	if (iInfo > 0) {
		return 0.0;
	} 
	
	bool bPositive = true;
	double dDet = 1.0;
	for(int i=0 ; i<m_iCols ; ++i) {
		dDet *= (*this)(i,i);
		if (iPivot[i] != (i+1)) {
			bPositive = !bPositive;
		}
	}
	
	// inverse computation based on the LU decomposition
	integer iLWork = iRows;
	double *dPWork = new double[iLWork];
	dgetri_(&iRows,m_rData,&iLDA,iPivot,dPWork,&iLWork,&iInfo);	
	
	delete [] iPivot;
	delete [] dPWork;
	
	// return the determinant
	return bPositive ? dDet : -dDet;
}


// compute the matrix of cofactors
template<typename Real>
void MatrixBase<Real>::matrixCofactor() {
	
	Real rDet = invert();
	mul(rDet);
	transpose();
}
		
// invert matrix elements
template<typename Real>
void MatrixBase<Real>::invertElements() {
	
	for(int i=0 ; i < m_iRows ; ++i) {
		for(int j=0 ; j < m_iCols ; ++j) {
			m_rData(i,j) = 1.0/m_rData(i,j);
		}
	}
}

// multiplication
template<typename Real>
void MatrixBase<Real>::add() {

	assert(0);
}

// singular value decomposition
template<typename Real>
void MatrixBase<Real>::svd(VectorBase<Real> vEigenvalues, MatrixBase<Real> mU, MatrixBase<Real> mV) {

}	

// print the vector
template<typename Real>
void MatrixBase<Real>::print() {
	
	cout << m_iRows << " x " << m_iCols << endl;
	for(int i=0 ; i < m_iRows ; ++i) {
		for(int j=0 ; j < m_iCols ; ++j) {
			cout << (*this)(i,j) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// write the matrix
template<typename Real>
void MatrixBase<Real>::write(ostream &os) {

	IOBase::write(os,m_iRows);
	IOBase::write(os,m_iCols);
	writeData(os);
}

// write the matrix data
template<typename Real>
void MatrixBase<Real>::writeData(ostream &os) {

	// write the data row by row since the # of columns can be different than the stride
	for(int i=0 ; i<m_iRows ; ++i) {
		IOBase::writeBytes(os,reinterpret_cast<char*>(getRowData(i)),m_iCols*sizeof(Real));
	}
}

template void MatrixBase<float>::copyRows(VectorBase<float> &v);
template void MatrixBase<float>::copyCols(VectorBase<float> &v);
template VectorStatic<float> MatrixBase<float>::getRow(int iRow);
template VectorStatic<double> MatrixBase<double>::getRow(int iRow);
template float *MatrixBase<float>::getRowData(int iRow);
template double *MatrixBase<double>::getRowData(int iRow);
template void MatrixBase<float>::zero();
template void MatrixBase<double>::zero();
template bool MatrixBase<float>::equal(MatrixBase<float> &m, float rEpsilon);
template void MatrixBase<float>::matrixCofactor();
template void MatrixBase<double>::matrixCofactor();
template void MatrixBase<float>::add(float r, MatrixBase<float> &matrix);
template void MatrixBase<double>::add(double r, MatrixBase<double> &matrix);
template void MatrixBase<float>::addMul(MatrixBase<float> &mA, MatrixTransposed transA, 
	MatrixBase<float> &mB, MatrixTransposed transB, 
	MatrixBase<float> &mC, MatrixTransposed transC);
template void MatrixBase<float>::addVecMul(float r, VectorBase<float> &v1, VectorBase<float> &v2);
template void MatrixBase<double>::addVecMul(double r, VectorBase<double> &v1, VectorBase<double> &v2);
template void MatrixBase<float>::mul(float r);
template void MatrixBase<double>::mul(double r);
template void MatrixBase<float>::mulRows(VectorBase<float> &v);
template void MatrixBase<double>::mulRows(VectorBase<double> &v);
template void MatrixBase<float>::mulCols(VectorBase<float> &v);
template void MatrixBase<double>::mulCols(VectorBase<double> &v);
template void MatrixBase<float>::print();
template void MatrixBase<double>::print();
template void MatrixBase<float>::write(ostream &os);
template void MatrixBase<double>::write(ostream &os);
template void MatrixBase<float>::writeData(ostream &os);
template void MatrixBase<double>::writeData(ostream &os);

};	// end-of-namespace


