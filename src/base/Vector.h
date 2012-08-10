/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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

/**
	@author daniel <dani.bolanos@gmail.com>
*/
template<typename Real>
class Vector {

	private: 
	
		int m_iDim;			// dimensionality
		Real *m_rData;		// data

	public:

		// constructor
		Vector() : m_iDim(0), m_rData(NULL) {
		}

		// constructor
		Vector(int iDim) : m_iDim(iDim), m_rData(NULL) {
			allocate();
			zero();
		}

		// destructor
		~Vector() {
			if (m_rData) {
				deallocate();
			}
		}	
		
		// return the dimensionality
		int getDim() {
		
			return m_iDim;
		}
		
		// return the data
		Real *getData() {
		
			return m_rData;
		}
		
		// allocate memory
		void allocate() {
		
			assert(m_iDim > 0);
			m_rData = new Real[m_iDim];
		}
		
		// deallocate memory
		void deallocate() {
		
			assert(m_rData);
			delete [] m_rData;
		}
		
		// resize the matrix (also used to allocate space for the first time)
		void resize(int iDim) {
		
			assert(iDim > 0);
		
			if (iDim != m_iDim) {
			
				if (m_rData) {
					deallocate();
				}
				
				m_iDim = iDim;
				allocate();	
			}
			
			zero();
		}
		
		// adds another vector
		void add(Vector<Real> &v) {
			
			assert(m_iDim == v.getDim());
		 
			for(int i=0 ; i<m_iDim ; ++i) {
				m_rData[i] += v.m_rData[i];
			}
		}
		
		// set elements to zero
		void zero() {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] = 0.0;
			}
		}
		
		// scale elements by the given factor
		void scale(Real rFactor) {
		
			for(int i=0 ; i < m_iDim ; ++i) {
				m_rData[i] *= rFactor;
			}
		}

};

#endif
