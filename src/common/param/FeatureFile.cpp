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
	
#ifdef SIMD
	
	// find the smaller multiple of 16 that is equal or above the feature vector dimensionality
	m_iDimAligned16 = m_iDim;
	int iAux = (iFeatureDimensionality*sizeof(float))%16;
	if (iAux > 0) {
		m_iDimAligned16 += (16-iAux)/sizeof(float);
	}
	
#endif
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
	int iFeatureVectorSize = m_iDim*sizeof(float);
	int iBytes = file.size();
	if (m_iFormat == FORMAT_FEATURES_FILE_DEFAULT) {		
		assert(iBytes % iFeatureVectorSize == 0);	
		m_iFeatureVectors = iBytes/iFeatureVectorSize;
	} else {
		assert(m_iFormat == FORMAT_FEATURES_FILE_HTK);
		assert((iBytes-12) % iFeatureVectorSize == 0);	
		m_iFeatureVectors = (iBytes-12) / iFeatureVectorSize;
		file.getStream().seekg(12);	
	}	
	
	// allocate memory for the feature vectors and read them
#ifdef SIMD	

	// we need the data aligned to addresses multiple of 16 bytes so they can be loaded in the sse registers more efficiently
	int iReturnValue = posix_memalign((void**)&m_fFeatureVectors,sizeof(__m128i),
		m_iFeatureVectors*m_iDimAligned16*sizeof(float));
	assert(iReturnValue == 0);
	
 	// we need to read them one by one to preserve the memory alignment
	for(int i=0 ; i<m_iFeatureVectors ; ++i) {
		int iOffset = i*m_iDimAligned16;
		IOBase::readBytes(file.getStream(),reinterpret_cast<char*>(m_fFeatureVectors+iOffset),
			m_iDim*sizeof(float));	
	}
	
#else

	// read all the feature vectors at once
	m_fFeatureVectors = new float[m_iFeatureVectors*m_iDim];
	IOBase::readBytes(file.getStream(),reinterpret_cast<char*>(m_fFeatureVectors),
		m_iFeatureVectors*m_iDim*sizeof(float));
	
#endif

	file.close();
}

// store the features 
void FeatureFile::store(float *fFeatureVectors, unsigned int iFeatureVectors) {

	// check that the object is in the right mode
	assert(m_iMode == MODE_WRITE);
	
	FileOutput file(m_strFile.c_str(),true);
	file.open();
	
#ifdef SIMD

 	// we need to write them one by one because of the memory alignment
	for(int i=0 ; i<m_iFeatureVectors ; ++i) {
		int iOffset = i*m_iDimAligned16;
		IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(m_fFeatureVectors+iOffset),
			m_iDim*sizeof(float));	
	}
	
#else

	IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(fFeatureVectors),
		iFeatureVectors*m_iDim*sizeof(float));

#endif
	
	file.close();
}

// return a reference to the features
float *FeatureFile::getFeatureVectors(unsigned int *iFeatureVectors) {

	*iFeatureVectors = m_iFeatureVectors;
	
	return m_fFeatureVectors;
}

// print the features (debugging)
void FeatureFile::print(float *fFeatures, unsigned int iFeatures, unsigned int iDim, unsigned int iDelta) {

	assert((iDim%iDelta) == 0);

	cout << endl;
	for(unsigned int i=0 ; i < iFeatures ; ++i) {
		for(unsigned int j=0 ; j < iDelta ; ++j) {
			for(unsigned int h=0 ; h < iDim/iDelta ; ++h) {
				cout << " " << FLT(9,6) << fFeatures[i*iDim+j*(iDim/iDelta)+h];
			}
			cout << endl;
		}
		cout << endl;	
	}	
}

};	// end-of-namespace
