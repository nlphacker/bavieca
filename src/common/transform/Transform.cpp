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


#include "Transform.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"

namespace Bavieca {

// constructor
Transform::Transform()
{
	m_matrix = NULL;
}

// constructor
Transform::Transform(int iType, Matrix<float> &matrix)
{
	m_iType = iType;
	m_matrix = new Matrix<float>(matrix);
}

// destructor
Transform::~Transform()
{
	if (m_matrix) {
		delete m_matrix;
	}
}

// load the transform from disk
void Transform::load(const char *strFile) {

	FileInput file(strFile,true);
	file.open();
	
	IOBase::read(file.getStream(),&m_iType);
	m_matrix = Matrix<float>::read(file.getStream());
	
	file.close();
}

// store the transform to disk
void Transform::store(const char *strFile) {

	FileOutput file(strFile,true);
	file.open();
		
	IOBase::write(file.getStream(),m_iType);
	m_matrix->write(file.getStream());
		
	file.close();
}

// print the transform information
void Transform::print(bool bPrintData) {
	
	assert(m_matrix);
	cout << "- transform -----------" << endl;
	cout << "type: " << getStrTransformType();
	m_matrix->print();
}

// apply the transform
void Transform::apply(VectorBase<float> &vInput, VectorBase<float> &vOutput) {
	
	// linear transform (columns in transform == feature dimension)
	if (m_iType == TRANSFORM_TYPE_LINEAR) {
	
		vOutput.mul(*m_matrix,vInput);
	} 
	// affine transform (columns in transform == feature dimension+1, append 1 to feature vector)
	else {
		assert(m_iType == TRANSFORM_TYPE_AFFINE);
		
		// create the extended vector
		Vector<float> vInputEx(vInput);
		vInputEx.appendFront(1.0);
		// apply transform
		vOutput.mul(*m_matrix,vInputEx);
	}
}


// apply the transform to a series of features
Matrix<float> *Transform::apply(MatrixBase<float> &mFeatures) {

	Matrix<float> *mFeaturesX = new Matrix<float>(mFeatures.getRows(),m_matrix->getRows());
	for(unsigned int i=0 ; i < mFeatures.getRows() ; ++i) {
		VectorStatic<float> vFeatureVector = mFeatures.getRow(i);
		VectorStatic<float> vFeatureVectorX = mFeaturesX->getRow(i);
		apply(vFeatureVector,vFeatureVectorX);
	}
	
	return mFeaturesX;
}


};	// end-of-namespace

