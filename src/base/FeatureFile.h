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

#ifndef FEATUREFILE_H
#define FEATUREFILE_H

#include <stdio.h>
#include <stdlib.h>

#include <xmmintrin.h>
//#include <smmintrin.h>
//#include <pmmintrin.h>
#include <malloc.h>

using namespace std;

#include <string>
#include "Global.h"

#include "FeatureExtractor.h"

#define MODE_READ			0
#define MODE_WRITE		1

// features file format
#define FORMAT_FEATURES_FILE_DEFAULT			0
#define FORMAT_FEATURES_FILE_HTK					1

/**
	@author Daniel Bolanos <bolanos@cslr.colorado.edu>
*/
class FeatureFile {

	private:
	
		string m_strFile;						// file name
	#ifdef SIMD		
		__attribute__((aligned(16))) float *m_fFeatureVectors;		// aligned to 16 bytes to speed up memory access
		int m_iFeatureDimensionalityAligned16;		// dimensionality needed to store a 16 bytes aligned feature vector
	#else
		float *m_fFeatureVectors;			// feature vectors
	#endif
		int m_iFeatureVectors;				// # feature vectors
		char m_iMode;							// mode
		char m_iFormat;						// file format (default, HTK, etc)
		int m_iFeatureDimensionality;		// dimensionality

	public:
    	
    	// constructor
		FeatureFile(const char *strFileName, const char iMode, const char iFormat = FORMAT_FEATURES_FILE_DEFAULT, int iFeatureDimensionality = FEATURE_DIMENSIONALITY_DEFAULT);

		// destructor
		~FeatureFile();
		
		// load the features
		void load();
		
		// store the features 
		void store(float *fFeatureVectors, int iFeatureVectors);
		
		// return a reference to the features
		float *getFeatureVectors(int *iFeatureVectors);		
		
		// print the features (debugging)
		static void print(float *fFeatures, int iFeatures);
};

#endif
