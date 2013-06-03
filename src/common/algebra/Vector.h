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


#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <assert.h>

#include "LogMessage.h"
#include "VectorBase.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<class Real>
class Vector : public VectorBase<Real> {

	public:

		// constructor
		Vector() : VectorBase<Real>::m_iDim(0), VectorBase<Real>::m_rData(NULL) {
		}
		
		// constructor
		Vector(int iDim) : VectorBase<Real>(iDim) {
		
			allocate();
			this->zero();
		}
		
		// constructor
		Vector(const Vector<Real> &v) : VectorBase<Real>(v.getDim()) {
		
			allocate();
			this->copy(v);
		}	
		
		// constructor
		Vector(const VectorBase<Real> &v) : VectorBase<Real>(v.getDim()) {
			
			allocate();
			this->copy(v);
		}	
		
		// constructor
		template<typename Real2>
		Vector(const VectorBase<Real2> &v) : VectorBase<Real>(v.getDim()) {
		
			allocate();
			VectorBase<Real>::copy(v);
		}		
		
		// destructor
		~Vector() {
			if (this->m_rData) {
				deallocate();	
			}
		}
			
		// allocate memory
		void allocate() {
		
			assert(this->m_iDim > 0);
		#if defined __AVX__ || defined __SSE__
			if (posix_memalign((void**)&this->m_rData,ALIGN_BOUNDARY,this->m_iDim*sizeof(Real)) !=0) {
				BVC_ERROR << "memory allocation error, unable to allocate aligned memory using posix_memalign";
			}	
		#else
			this->m_rData = new Real[this->m_iDim];
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
		
		// resize the vector (also used to allocate space for the first time)
		void resize(unsigned int iDim, bool bCopy = false) {
		
			assert(iDim > 0);
		
			unsigned int iDimOrig = this->m_iDim;
			Real *rDataCopy = this->m_rData;
			if (iDim != this->m_iDim) {
				
				this->m_iDim = iDim;
				allocate();
				if (bCopy) {
					memcpy(this->m_rData,rDataCopy,iDimOrig*sizeof(Real));
				}
				delete [] rDataCopy;
			}
		}
		
		// extend the value to the beginning of the vector
		void appendFront(Real r) {
		
			resize(this->m_iDim+1,true);
			VectorBase<Real>::shiftRight(1);
			this->m_rData[0] = r;	
		}
};

};	// end-of-namespace

#endif
