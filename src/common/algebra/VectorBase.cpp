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


#include "IOBase.h"
#include "VectorBase.h"
#include "MatrixBase.h"

#include <math.h>

namespace Bavieca {

// copy from memory
template<typename Real>
void VectorBase<Real>::copy(Real *rData, unsigned int iDim) {

	assert(m_iDim == iDim);
	memcpy(m_rData,rData,m_iDim*sizeof(Real));
}

// apply the square root
template<typename Real>
void VectorBase<Real>::sqrt() {

	for(unsigned int i=0 ; i < m_iDim ; ++i) {
		m_rData[i] = ::sqrt(m_rData[i]);
	}
}

// shift the vector elements to the right
template<typename Real>
void VectorBase<Real>::shiftRight(unsigned int iOffset) {
	
	for(unsigned int i=m_iDim-1 ; i >= iOffset ; --i) {
		m_rData[i] = m_rData[i-iOffset];
	}
}

// print the vector
template<typename Real>
void VectorBase<Real>::print() {
	
	for(unsigned int i=0 ; i < m_iDim ; ++i) {
		cout << m_rData[i] << " ";
	}
	cout << endl;
}

// write the vector
template<typename Real>
void VectorBase<Real>::write(ostream &os) {

	IOBase::write(os,m_iDim);
	writeData(os);
}

// write the vector data
template<typename Real>
void VectorBase<Real>::writeData(ostream &os) {

	IOBase::writeBytes(os,reinterpret_cast<char*>(m_rData),m_iDim*sizeof(Real));
}

// read the vector data
template<typename Real>
void VectorBase<Real>::readData(istream &is) {

	IOBase::readBytes(is,reinterpret_cast<char*>(m_rData),m_iDim*sizeof(Real));
}

template void VectorBase<float>::copy(float *rData, unsigned int iDim);
template void VectorBase<double>::copy(double *rData, unsigned int iDim);
template void VectorBase<float>::sqrt();
template void VectorBase<double>::sqrt();
template void VectorBase<float>::shiftRight(unsigned int iOffset);
template void VectorBase<double>::shiftRight(unsigned int iOffset);
template void VectorBase<float>::print();
template void VectorBase<double>::print();
template void VectorBase<float>::write(ostream &os);
template void VectorBase<double>::write(ostream &os);
template void VectorBase<float>::writeData(ostream &os);
template void VectorBase<double>::writeData(ostream &os);
template void VectorBase<float>::readData(istream &is);
template void VectorBase<double>::readData(istream &is);

};	// end-of-namespace

