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


#ifndef FEATUREFILE_H
#define FEATUREFILE_H

#include <stdio.h>
#include <stdlib.h>

// support for SIMD instructions
#ifdef __SSE__
#include <xmmintrin.h>
#endif

// support for AVX (Intel Advanced Vector Extensions) instructions
#ifdef __AVX__
#include <immintrin.h>
#endif

#if defined __linux__ || defined __MINGW32__ || defined _MSC_VER
#include <malloc.h>
#elif __APPLE__
#include <sys/malloc.h>
#else
	#error "unsupported platform"
#endif

using namespace std;

#include <string>
#include "Global.h"
#include "Matrix.h"

namespace Bavieca {

#define MODE_READ		0
#define MODE_WRITE		1

// default feature dimensionality
#define FEATURE_DIMENSIONALITY_DEFAULT			39

// features file format
#define FORMAT_FEATURES_FILE_DEFAULT			0
#define FORMAT_FEATURES_FILE_HTK				1

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FeatureFile {

	private:
	
		string m_strFile;						// file name
		Matrix<float> *m_mFeatures;		// feature vectors (each row in the matrix is a feature vector)
		char m_iMode;							// mode
		char m_iFormat;						// file format (default, HTK, etc)
		int m_iDim;								// dimensionality

	public:
    	
    	// constructor
		FeatureFile(const char *strFileName, const char iMode, const char iFormat = FORMAT_FEATURES_FILE_DEFAULT, int iDim = FEATURE_DIMENSIONALITY_DEFAULT);

		// destructor
		~FeatureFile();
		
		// load the features
		void load();
		
		// store the features 
		void store(MatrixBase<float> &mFeatures);
		
		// return a reference to the features
		Matrix<float> *getFeatureVectors();
		
		// print the features (debugging)
		static void print(MatrixBase<float> &mFeatures, unsigned int iDelta);
};

};	// end-of-namespace


#endif
