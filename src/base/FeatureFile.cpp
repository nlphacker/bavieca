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

#include "FeatureFile.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"

// constructor
FeatureFile::FeatureFile(const char *strFile, const char iMode, const char iFormat, int iFeatureDimensionality) {

	m_iMode = iMode;
	m_strFile = strFile;
	m_iFormat = iFormat;
	m_iFeatureDimensionality = iFeatureDimensionality;
	
#ifdef SIMD
	
	// find the smaller multiple of 16 that is equal or above the feature vector dimensionality
	m_iFeatureDimensionalityAligned16 = m_iFeatureDimensionality;
	int iAux = (iFeatureDimensionality*sizeof(float))%16;
	if (iAux > 0) {
		m_iFeatureDimensionalityAligned16 += (16-iAux)/sizeof(float);
	}
	
#endif

	return;
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
	int iFeatureVectorSize = m_iFeatureDimensionality*sizeof(float);
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
		m_iFeatureVectors*m_iFeatureDimensionalityAligned16*sizeof(float));
	assert(iReturnValue == 0);
	
 	// we need to read them one by one to preserve the memory alignment
	for(int i=0 ; i<m_iFeatureVectors ; ++i) {
		int iOffset = i*m_iFeatureDimensionalityAligned16;
		IOBase::readBytes(file.getStream(),reinterpret_cast<char*>(m_fFeatureVectors+iOffset),
			m_iFeatureDimensionality*sizeof(float));	
	}
	
#else

	// read all the feature vectors at once
	m_fFeatureVectors = new float[m_iFeatureVectors*m_iFeatureDimensionality];
	IOBase::readBytes(file.getStream(),reinterpret_cast<char*>(m_fFeatureVectors),
		m_iFeatureVectors*m_iFeatureDimensionality*sizeof(float));
	
#endif

	file.close();
}

// store the features 
void FeatureFile::store(float *fFeatureVectors, int iFeatureVectors) {

	// check that the object is in the right mode
	assert(m_iMode == MODE_WRITE);
	
	FileOutput file(m_strFile.c_str(),true);
	file.open();
	
#ifdef SIMD

 	// we need to write them one by one because of the memory alignment
	for(int i=0 ; i<m_iFeatureVectors ; ++i) {
		int iOffset = i*m_iFeatureDimensionalityAligned16;
		IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(m_fFeatureVectors+iOffset),
			m_iFeatureDimensionality*sizeof(float));	
	}
	
#else

	IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(fFeatureVectors),
		m_iFeatureDimensionality*sizeof(float));

#endif
	
	file.close();
}

// return a reference to the features
float *FeatureFile::getFeatureVectors(int *iFeatureVectors) {

	*iFeatureVectors = m_iFeatureVectors;
	
	return m_fFeatureVectors;
}

// print the features (debugging)
void FeatureFile::print(float *fFeatures, int iFeatures) {

	printf("\n");
	for(int i=0 ; i < iFeatures ; ++i) {
		for(int j=0 ; j < 3 ; ++j) {
			for(int h=0 ; h < 13 ; ++h) {
				printf(" %9.6f",fFeatures[i*(13*3)+j*13+h]);
			}
			printf("\n");
		}
		printf("\n");	
	}	
	
	return;
}
