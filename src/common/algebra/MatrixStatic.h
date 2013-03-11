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


#ifndef MATRIXSTATIC_H
#define MATRIXSTATIC_H

#include "MatrixBase.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real>
class MatrixStatic : public MatrixBase<Real> {

	public:
	
		// constructor from another matrix
		MatrixStatic(MatrixBase<Real> &m, int iRowStart, int iRows, int iColStart, int iCols) : MatrixBase<Real>(iRows,iCols){
		
			assert((iRowStart >= 0) && (iRowStart+iRows <= m.getRows()));
			assert((iColStart >= 0) && (iColStart+iCols <= m.getCols()));
			assert((iRows >= 0) && (iCols >= 0));
		
			this->m_rData = m.getData()+iRowStart*m.getStride()+iColStart;
			this->m_iStride = m.getStride();
		}

		// constructor
		MatrixStatic(Real *rData, int iRows, int iCols);

		// constructor (square matrices)
		MatrixStatic(Real *rData, int iDim);

		// destructor
		~MatrixStatic() {}

};

};	// end-of-namespace

#endif
