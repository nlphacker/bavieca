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

#include <iomanip>
#include "FeatureFile.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"

namespace Bavieca {

// constructor
FeatureFile::FeatureFile(const char *strFile, const char iMode, const char iFormat, int iDim) {

	m_iMode = iMode;
	m_strFile = strFile;
	m_iFormat = iFormat;
	m_iDim = iDim;
	m_mFeatures = NULL;
}

// destructor
FeatureFile::~FeatureFile() {

}

// load the features
void FeatureFile::load() {

	// check that the object is in the right mode
	assert(m_iMode == MODE_READ);
	
	FileInput file(m_strFile.c_str(),true);
	file.open();
	
	// get the number of feature vectors
	int iFeatureVectors = -1;
	int iFeatureVectorSize = m_iDim*sizeof(float);
	int iBytes = file.size();
	if (m_iFormat == FORMAT_FEATURES_FILE_DEFAULT) {		
		assert(iBytes % iFeatureVectorSize == 0);	
		iFeatureVectors = iBytes/iFeatureVectorSize;
	} else {
		assert(m_iFormat == FORMAT_FEATURES_FILE_HTK);
		assert((iBytes-12) % iFeatureVectorSize == 0);	
		iFeatureVectors = (iBytes-12) / iFeatureVectorSize;
		file.getStream().seekg(12);	
	}
	
	// read the feature vectors
	m_mFeatures = new Matrix<float>(iFeatureVectors,m_iDim);
	m_mFeatures->readData(file.getStream());

	file.close();
}

// store the features 
void FeatureFile::store(MatrixBase<float> &mFeatures) {

	// check that the object is in the right mode
	assert(m_iMode == MODE_WRITE);
	
	FileOutput file(m_strFile.c_str(),true);
	file.open();
	
	mFeatures.writeData(file.getStream());	
	
	file.close();
}

// return a reference to the features
Matrix<float> *FeatureFile::getFeatureVectors() {
	
	return m_mFeatures;
}

// print the features (debugging)
void FeatureFile::print(MatrixBase<float> &mFeatures, unsigned int iDelta) {
	
	assert((mFeatures.getCols()%iDelta) == 0);

	cout << endl;
	for(unsigned int i=0 ; i < mFeatures.getRows() ; ++i) {
		for(unsigned int j=0 ; j < iDelta ; ++j) {
			for(unsigned int h=0 ; h < mFeatures.getCols()/iDelta ; ++h) {
				cout << " " << FLT(9,6) << mFeatures(i,h+j*(mFeatures.getCols()/iDelta));
			}
			cout << endl;
		}
		cout << endl;	
	}	
}

};	// end-of-namespace
