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


#ifndef HMMSTATEDECODING_H
#define HMMSTATEDECODING_H

using namespace std;

#include <vector>
#include <list>

#include <stdio.h>
#include <string.h>
#include <math.h>

// support for SIMD instructions
#ifdef __SSE__
#include <xmmintrin.h>
#endif

// support for AVX (Intel Advanced Vector Extensions) instructions
#ifdef __AVX__
#include <immintrin.h>
#endif

// for gcc #include <x86intrin.h> can be used instead, it includes whatever is needed

#if defined __linux__ || defined __MINGW32__ || defined _MSC_VER
#include <malloc.h>
#elif __APPLE__
#include <sys/malloc.h>
#else 
	#error "unsupported platform"
#endif

#include "Global.h"

namespace Bavieca {

class FileInput;
class FileOutput;
class IOBase;
class PhoneSet;

// state's within-word position
#define WITHIN_WORD_POSITION_START			0
#define WITHIN_WORD_POSITION_INTERNAL		1
#define WITHIN_WORD_POSITION_END				2
#define WITHIN_WORD_POSITION_MONOPHONE		3

#define STR_WITHIN_WORD_POSITION_START			"start"
#define STR_WITHIN_WORD_POSITION_INTERNAL		"internal"
#define STR_WITHIN_WORD_POSITION_END			"end"
#define STR_WITHIN_WORD_POSITION_MONOPHONE	"monophone"

class HMMStateDecoding;

typedef vector<HMMStateDecoding*> VHMMStateDecoding;

#define DIMENSIONALITY 39

typedef struct {
	int iId;											// unique identifier (across all Gaussian component in the system)
	float fWeight;											// gaussian weight
	float fConstant;
#if defined __AVX__ || defined __SSE__
	#if defined __linux__ || defined __APPLE__ || __MINGW32__
	__attribute__((aligned(ALIGN_BOUNDARY))) float fMean[DIMENSIONALITY];			// gaussian mean
	__attribute__((aligned(ALIGN_BOUNDARY))) float fCovariance[DIMENSIONALITY];	// gaussian covariance matrix (diagonal)
	#elif _MSC_VER
	__declspec(align(ALIGN_BOUNDARY)) float fMean[DIMENSIONALITY];						// gaussian mean
	__declspec(align(ALIGN_BOUNDARY)) float fCovariance[DIMENSIONALITY];				// gaussian covariance matrix (diagonal)
	#endif
#else
	float fMean[DIMENSIONALITY];				// gaussian mean
	float fCovariance[DIMENSIONALITY];		// gaussian covariance matrix (diagonal)	
#endif	
} GaussianDecoding;

typedef vector<GaussianDecoding*> VGaussianDecoding;
typedef list<GaussianDecoding*> LGaussianDecoding;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class HMMStateDecoding {

	private:
	
		int m_iDim;
		PhoneSet *m_phoneSet;

		// state identity
		unsigned char m_iPhone;						// basephone
		unsigned char m_iState;						// HMM-state 
		unsigned char m_iPosition;					// HMM-state within-word position
		int m_iId;										// HMM-state unique identifier
		
		// gaussians
		int m_iGaussianComponents;
		GaussianDecoding *m_gaussians;			// block of memory containing all gaussians
		
		// constants
		float m_fConstant;
		
		// caching computations
		int m_iTimestamp;							// time-stamp (for caching emission probability computations)
		float m_fProbabilityCached;			// cached probability
		
		// whether the covariance was modified to accelerate the computation of emission probabilities
		bool m_bCovarianceOriginal;
		
		// computes the emission probability of the state given the feature vector ("brute force")
		float computeEmissionProbabilityBruteForce(float *fFeatures, int iTime);
		
		// computes the emission probability of the state given the feature vector 
		// uses nearest-neighbor approximation
		float computeEmissionProbabilityNearestNeighbor(float *fFeatures, int iTime);	
		
		// computes the emission probability of the state given the feature vector 
		// uses nearest-neighbor approximation
		// uses Partial Distance Elimination (PDE)
		float computeEmissionProbabilityNearestNeighborPDE(float *fFeatures, int iTime);	
		
		// computes the emission probability of the state given the feature vector 
		// uses nearest-neighbor approximation
		// uses SIMD instructions (SSE) (sse support must be enabled during compilation!)
		#ifdef __SSE__
		float computeEmissionProbabilityNearestNeighborSSE(float *fFeatures, int iTime);	
		#endif
				
		// computes the emission probability of the state given the feature vector 
		// uses nearest-neighbor approximation
		// uses Partial Distance Elimination (PDE)
		// uses SIMD instructions (SSE) (sse support must be enabled during compilation!)
		#ifdef __SSE__
		float computeEmissionProbabilityNearestNeighborPDE_SSE(float *fFeatures, int iTime);
		#endif
		
		// computes the emission probability of the state given the feature vector 
		// uses nearest-neighbor approximation
		// uses AVX (Intel® Advanced Vector Extensions) instructions (avx support must be enabled during compilation!)
		// compared to SSE, AVX offers 256 bit registers instead of 128 bit registers 
		#ifdef __AVX__	
		float computeEmissionProbabilityNearestNeighborAVX(float *fFeatures, int iTime);
		#endif	

	public:

		// constructor
		HMMStateDecoding(int iDim, PhoneSet *phoneSet, int iId);
		
		// constructor
		HMMStateDecoding(int iDim, PhoneSet *phoneSet, unsigned char iPhone, unsigned char iState, 
			unsigned char iPosition, int iId, int iGaussians, GaussianDecoding *gaussians);
		
		// default constructor
		HMMStateDecoding();
		
		// set initial parameters
		void setInitialParameters(int iDim, PhoneSet *phoneSet, int iId);

		// destructor
		~HMMStateDecoding();
		
		// initialize the estimation by precomputing constants and invariant terms
		void initialize();
		
		// load the HMM from a file
		void load(FileInput &file);
		
		// store the HMM into a file
		void store(FileOutput &file);	
		
		inline void resetTimeStamp() {
		
			m_iTimestamp = -1;
		}
				
		inline unsigned char getPhone() {
		
			return m_iPhone;
		}
		
		inline unsigned char getState() {
		
			return m_iState;
		}
		
		inline unsigned char getPosition() {
		
			return m_iPosition;
		}
		
		inline int getId() {
		
			return m_iId;
		}
		
		// return the array of gaussians
		inline GaussianDecoding* getGaussians(int &iGaussians) {
		
			iGaussians = m_iGaussianComponents;
			
			return m_gaussians;
		}	
		
		// return a Gaussian component by its index
		inline GaussianDecoding *getGaussian(int iGaussian) {
		
			assert((iGaussian >= 0) && (iGaussian < m_iGaussianComponents));
			
			return &m_gaussians[iGaussian];
		}
		
		// return the number of Gaussian components in the mixture
		inline unsigned int getMixtureSize() {
		
			return m_iGaussianComponents;
		}
		
		inline int getGaussianComponents() {
		
			return m_iGaussianComponents;
		}
		
		// creates a copy of the Gaussian
		inline static void copyGaussian(GaussianDecoding *gDest, GaussianDecoding *gSource, unsigned int iDim) {
		
			gDest->fWeight = gSource->fWeight;
			gDest->fConstant = gSource->fConstant;
			for(unsigned int i = 0 ; i < iDim ; ++i) {
				gDest->fMean[i] = gSource->fMean[i];
				gDest->fCovariance[i] = gSource->fCovariance[i];
			}
		}
		
		// comparison function to sort Gaussian components by weight
		static bool compareGaussianByWeight(const GaussianDecoding *g1, const GaussianDecoding *g2) {
		
			return (g1->fWeight > g2->fWeight);
		}
		
		// computes the emission probability of the state given the feature vector
		inline float computeEmissionProbability(float *fFeatures, int iTime) {
		
		#ifdef __AVX__
			return computeEmissionProbabilityNearestNeighborAVX(fFeatures,iTime);
		#elif __SSE__
			return computeEmissionProbabilityNearestNeighborSSE(fFeatures,iTime);
		#else	
			return computeEmissionProbabilityNearestNeighborPDE(fFeatures,iTime);	
		#endif
		}
		
		// return the best scoring gaussian for a given feature vector
		GaussianDecoding *getBestScoringGaussian(float *fFeatures, float *fScore);	
			
		// return the Gaussian likelihood for the given feature vector
		// uses nearest-neighbor approximation
		// uses Partial Distance Elimination (PDE)
		double computeGaussianProbability(int iGaussian, float *fFeatures);
		
		inline void setGaussianIds(int *iId) {
		
			for(int i=0 ; i < m_iGaussianComponents ; ++i) {
				m_gaussians[i].iId = (*iId)++;
			}	
		}		
};

};	// end-of-namespace

#endif
