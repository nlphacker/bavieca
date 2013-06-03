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


#include "Matrix.h"
#include "IOBase.h"

namespace Bavieca {

// read the matrix
template<typename Real>
Matrix<Real> *Matrix<Real>::read(istream &is) {

	int iRows,iCols;
	IOBase::read(is,&iRows);
	IOBase::read(is,&iCols);
	Matrix<Real> *matrix = new Matrix<Real>(iRows,iCols);
	matrix->readData(is);
	
	return matrix;
}

// read the matrix data
template<typename Real> void Matrix<Real>::readData(istream &is) {

	// read the matrix row by row to preserve memory alignment
	for(unsigned int i=0 ; i < this->m_iRows ; ++i) {
		IOBase::readBytes(is,reinterpret_cast<char*>(this->m_rData+(i*this->m_iStride)),this->m_iCols*sizeof(Real));
	}
}

// read the matrix
template Matrix<float> *Matrix<float>::read(istream &is);
template Matrix<double> *Matrix<double>::read(istream &is);
template void Matrix<float>::readData(istream &is);
template void Matrix<double>::readData(istream &is);

};	// end-of-namespace

