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


#include "MatrixStatic.h"

namespace Bavieca {

// constructor
template<typename Real>
MatrixStatic<Real>::MatrixStatic(Real *rData, unsigned int iDim) : MatrixBase<Real>(iDim,iDim) {
	MatrixStatic(rData,iDim,iDim);
}

// constructor
template<typename Real>
MatrixStatic<Real>::MatrixStatic(Real *rData, unsigned int iRows, unsigned int iCols) : MatrixBase<Real>(iRows,iCols) {
	this->m_rData = rData;	
#if defined __AVX__ || defined __SSE__
	// check memory alignment
	assert(is_aligned(rData,ALIGN_BOUNDARY));
	// round up the dimensionality to the nearest multiple	
	unsigned int iReminder = this->m_iCols%(ALIGN_BOUNDARY/sizeof(Real));
	this->m_iStride = (iReminder == 0) ? this->m_iCols : this->m_iCols + (ALIGN_BOUNDARY/sizeof(Real))-iReminder;
	assert(this->m_iStride % (ALIGN_BOUNDARY/sizeof(Real)) == 0);
#else 
	this->m_iStride = this->m_iCols;
#endif	
}


template MatrixStatic<float>::MatrixStatic(float *rData, unsigned int iDim);
template MatrixStatic<double>::MatrixStatic(double *rData, unsigned int iDim);
template MatrixStatic<float>::MatrixStatic(float *rData, unsigned int iRows, unsigned int iCols);
template MatrixStatic<double>::MatrixStatic(double *rData, unsigned int iRows, unsigned int iCols);

};	// end-of-namespace

