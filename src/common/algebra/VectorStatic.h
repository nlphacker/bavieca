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


#ifndef VECTORSTATIC_H
#define VECTORSTATIC_H

#include "VectorBase.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<class Real>
class VectorStatic : public VectorBase<Real> {

	public:

		// constructor
		VectorStatic(Real *rData, int iDim) : VectorBase<Real>(iDim) {
			this->m_rData = rData;
		}

		// destructor
		~VectorStatic() {}

};

};	// end-of-namespace

#endif
