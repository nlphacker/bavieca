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

#include <xmmintrin.h>
//#include <smmintrin.h>
//#include <pmmintrin.h>

#if defined __linux__ || defined _WIN32
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

#ifdef SIMD

typedef struct {
	int iId;											// unique identifier (across all Gaussian component in the system)
	float fWeight;											// gaussian weight
	float fConstant;
	__attribute__((aligned(16))) float fMean[DIMENSIONALITY];						// gaussian mean
	__attribute__((aligned(16))) float fCovariance[DIMENSIONALITY];				// gaussian covariance matrix (diagonal)	
} GaussianDecoding;

#else

typedef struct {
	int iId;											// unique identifier (across all Gaussian component in the system)
	float fWeight;									// gaussian weight
	float fConstant;								// constant needed for efficient computation of emission probabilities
	float fMean[DIMENSIONALITY];				// gaussian mean
	float fCovariance[DIMENSIONALITY];		// gaussian covariance matrix (diagonal)	
} GaussianDecoding;

#endif

typedef vector<GaussianDecoding*> VGaussianDecoding;
typedef list<GaussianDecoding*> LGaussianDecoding;

/**
	@author root <root@localhost.localdomain>
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
			//gDest->iBaseClass = gSource->iBaseClass;
			//gDest->accumulator = gSource->accumulator;
		}
		
		// comparison function to sort Gaussian components by weight
		static bool compareGaussianByWeight(const GaussianDecoding *g1, const GaussianDecoding *g2) {
		
			return (g1->fWeight > g2->fWeight);
		}
		
		// computes the emission probability of the state given the feature vector ("brute force")
		// each of the gaussians is a multivariate normal distribution		
		float computeEmissionProbabilityBruteForce(float *fFeatures, int iTime);
		
		// computes the emission probability of the state given the feature vector 
		// each of the gaussians is a multivariate normal distribution
		// uses Partial Distance Elimination (PDE)
		float computeEmissionProbabilityNearestNeighbor(float *fFeatures, int iTime);	
		
		// computes the emission probability of the state given the feature vector 
		// each of the gaussians is a multivariate normal distribution
		// uses nearest-neighbor approximation
		// uses Partial Distance Elimination (PDE)
		float computeEmissionProbabilityNearestNeighborPDE(float *fFeatures, int iTime);	
		
		// computes the emission probability of the state given the feature vector 
		// each of the gaussians is a multivariate normal distribution
		// uses nearest-neighbor approximation
		// uses SIMD instructions (to compile sse support has to be enabled!)
		float computeEmissionProbabilityNearestNeighborSIMD(float *fFeatures, int iTime);
		
		// computes the emission probability of the state given the feature vector 
		// each of the gaussians is a multivariate normal distribution
		// uses nearest-neighbor approximation
		// uses SIMD instructions (to compile sse support has to be enabled!)
		/*float computeEmissionProbabilityNearestNeighborSIMD(float *fFeatures, int iTime) {
		
			if (iTime == m_iTimestamp) {	
				return m_fProbabilityCached;
			}
		
		#if defined __linux__ || defined __APPLE__
			__attribute__((aligned(16))) float *fMean,*fCovariance,fAcc;
			__attribute__((aligned(16))) float fLogLikelihood = LOG_LIKELIHOOD_FLOOR;	
			__attribute__((aligned(16))) float tmpf[4];
		#elif _WIN32
			__declspec(align(16)) float *fMean,*fCovariance,fAcc;
			__declspec(align(16)) float fLogLikelihood = LOG_LIKELIHOOD_FLOOR;	
			__declspec(align(16)) float tmpf[4];
		#endif
		
		
		__m128 tmp;
		__m128 ans;
		__m128 obs128;
		__m128 mean128;
		__m128 cov128;	
		
			for(int iGaussian = 0 ; iGaussian < this->m_iGaussianComponents ; ++iGaussian) {
			
				fMean = m_gaussians[iGaussian].fMean;
				fCovariance = m_gaussians[iGaussian].fCovariance;	
				fAcc = m_gaussians[iGaussian].fConstant;
				
				obs128  = _mm_load_ps(fFeatures);
				mean128 = _mm_load_ps(fMean);
				cov128  = _mm_load_ps(fCovariance);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				ans     = _mm_mul_ps(tmp,cov128);
			
				obs128  = _mm_load_ps(fFeatures+4);
				mean128 = _mm_load_ps(fMean+4);
				cov128  = _mm_load_ps(fCovariance+4);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+8);
				mean128 = _mm_load_ps(fMean+8);
				cov128  = _mm_load_ps(fCovariance+8);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+12);
				mean128 = _mm_load_ps(fMean+12);
				cov128  = _mm_load_ps(fCovariance+12);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+16);
				mean128 = _mm_load_ps(fMean+16);
				cov128  = _mm_load_ps(fCovariance+16);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+20);
				mean128 = _mm_load_ps(fMean+20);
				cov128  = _mm_load_ps(fCovariance+20);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+24);
				mean128 = _mm_load_ps(fMean+24);
				cov128  = _mm_load_ps(fCovariance+24);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+28);
				mean128 = _mm_load_ps(fMean+28);
				cov128  = _mm_load_ps(fCovariance+28);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_load_ps(fFeatures+32);
				mean128 = _mm_load_ps(fMean+32);
				cov128  = _mm_load_ps(fCovariance+32);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				obs128  = _mm_set_ps(0, fFeatures[38], fFeatures[37], fFeatures[36]);
				mean128 = _mm_set_ps(0, fMean[38], fMean[37], fMean[36]);
				cov128  = _mm_set_ps(0, fCovariance[38], fCovariance[37], fCovariance[36]);
				tmp     = _mm_sub_ps(obs128, mean128);
				tmp     = _mm_mul_ps(tmp,tmp);
				tmp     = _mm_mul_ps(tmp,cov128);
				ans     = _mm_add_ps(ans, tmp);
				
				//__attribute__((aligned(16))) float f = 0.0;
				//ans = _mm_dp_ps(tmp,cov128,0xF1);		
				//_mm_store_ss(&f,ans);	
				//fAcc -= f;
		
				_mm_store_ps(tmpf, ans);
				fAcc -= (tmpf[0]+tmpf[1]+tmpf[2]+tmpf[3]);
				
				fLogLikelihood = max(fAcc,fLogLikelihood);
			}
		
			// cache the probability
			m_iTimestamp = iTime;
			m_fProbabilityCached = fLogLikelihood;
			
			return fLogLikelihood;
		}*/
		
		
		// computes the emission probability of the state given the feature vector 
		// each of the gaussians is a multivariate normal distribution
		// uses nearest-neighbor approximation
		// uses Partial Distance Elimination (PDE)
		// uses SIMD instructions (to compile sse support has to be enabled!)
		float computeEmissionProbabilityNearestNeighborPDE_SIMD(float *fFeatures, int iTime);
		
		// computes the emission probability of the state given the feature vector 
		// each of the gaussians is a multivariate normal distribution
		// uses nearest-neighbor approximation
		// uses Partial Distance Elimination (PDE)
		// uses SIMD instructions (to compile sse support has to be enabled!)
		// uses ASM code instead of 
		float computeEmissionProbabilityNearestNeighborPDE_SIMD_ASM(float *fFeatures, int iTime);
		
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
