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


#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Global.h"
#include "Matrix.h"

using namespace std;

#include <vector>

#include <stdio.h>

namespace Bavieca {

// transform type
#define TRANSFORM_TYPE_LINEAR		0		// (for example: HLDA)
#define TRANSFORM_TYPE_AFFINE		1		// (for example: fMLLR)

#define TRANSFORM_TYPE_LINEAR_STR		"linear"
#define TRANSFORM_TYPE_AFFINE_STR		"affine"

class Transform;

typedef vector<Transform*> VTransform;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class Transform {

	private:
	
		Matrix<float> *m_matrix;
		int m_iType;

	public:

		// constructor
		Transform();
		
		// constructor
		Transform(int iType, Matrix<float> &matrix);

		// destructor
		~Transform();
		
		// return the type
		inline int getType() {
			
			return m_iType;
		}
		
		// return the number of rows
		inline int getRows() {
		
			return m_matrix->getRows();
		}
		
		// return the number of columns
		inline int getColumns() {
		
			return m_matrix->getCols();
		}
		
		// return the actual transform data
		inline Matrix<float> &getTransform() {
			
			return *m_matrix;
		}
				
		// load the transform from disk
		void load(const char *strFile);
		
		// store the transform to disk
		void store(const char *strFile);
		
		// print the transform information
		void print(bool bPrintData = false);
		
		// return the transform in string format
		inline const char *getStrTransformType() {
		
			if (m_iType == TRANSFORM_TYPE_LINEAR) {
				return TRANSFORM_TYPE_LINEAR_STR;
			} else {
				assert(m_iType == TRANSFORM_TYPE_AFFINE);
				return TRANSFORM_TYPE_AFFINE_STR;
			}
		}
		
		// apply the transform
		void apply(VectorBase<float> &vInput, VectorBase<float> &vOutput);
		
		// apply the transform to a series of features
		Matrix<float> *apply(MatrixBase<float> &mFeatures);
		
};

};	// end-of-namespace

#endif
