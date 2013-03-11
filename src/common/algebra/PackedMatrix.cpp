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


#include "PackedMatrix.h"
#include "IOBase.h"

namespace Bavieca {

// write the matrix
template<typename Real>
void PackedMatrix<Real>::write(ostream &os) {

	IOBase::write(os,m_iDim);
	writeData(os);
}

// write the matrix data
template<typename Real>
void PackedMatrix<Real>::writeData(ostream &os) {

	IOBase::writeBytes(os,reinterpret_cast<char*>(m_rData),getElements()*sizeof(Real));
}

// read the matrix data
template<typename Real>
void PackedMatrix<Real>::readData(istream &is) {

	IOBase::readBytes(is,reinterpret_cast<char*>(m_rData),getElements()*sizeof(Real));
}

template void PackedMatrix<float>::write(ostream &os);
template void PackedMatrix<double>::write(ostream &os);
template void PackedMatrix<float>::writeData(ostream &os);
template void PackedMatrix<double>::writeData(ostream &os);
template void PackedMatrix<float>::readData(istream &is);
template void PackedMatrix<double>::readData(istream &is);

};	// end-of-namespace
