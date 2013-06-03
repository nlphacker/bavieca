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


#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <float.h>
#include <math.h>

using namespace std;

#include <vector>
#include <list>

#include "Vector.h"
#include "Matrix.h"
#include "SMatrix.h"
#include "TMatrix.h"

namespace Bavieca {

// covariance modeling type
#define COVARIANCE_MODELLING_TYPE_DEFAULT			-1
#define COVARIANCE_MODELLING_TYPE_DIAGONAL		0
#define COVARIANCE_MODELLING_TYPE_FULL				1

// covariance modeling type
#define COVARIANCE_MODELLING_TYPE_DIAGONAL_STR		"diagonal"
#define COVARIANCE_MODELLING_TYPE_FULL_STR			"full"

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class Gaussian {

	private:
	
		int m_iDim;									// feature dimensionality
		int m_iCovarianceType;					// covariance modelling type
		int m_iTimestamp;
		float m_fProbabilityCached;		

		float m_fWeight;							// component weight
		Vector<float> *m_vMean;					// mean
		Vector<float> *m_vCovarianceDiag;	// covariance
		SMatrix<float> *m_mCovarianceFull;	// covariance
		
		// precomputed term to speed-up Gaussian evaluation
		double m_dConstant;
		double m_dLogConstant;	
		
		// auxiliar fields for mixture splitting/incrementing
		bool m_bMerged;							// whether the Gaussian was just created from merging 2 Gaussian components
		float m_fWeightMixtureIncrement;		// weight used for weight-based mixture increment		
		
		// 
		TMatrix<float> *m_mCholeskyFactor;
		
		// evaluate the Gaussian distribution
		inline float evaluateDiagonalCovariance(float *fFeatures, int iTime)  {
		
			assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL);
			
			double dAcc;
			double dDeterminant;
			float fProbability = 0.0;
			float fExponent;
			
			double dConstant = pow(PI_NUMBER*2.0,m_iDim/2.0);	
			
			dDeterminant = 1.0;
			fExponent = 0.0;
			
			for(int j = 0 ; j < m_iDim ; ++j) {
				// compute the determinant of the covariance matrix
				dDeterminant *= (*m_vCovarianceDiag)(j);
				// compute the numerator
				dAcc = fFeatures[j]-(*m_vMean)(j);
				fExponent += (float)(dAcc*dAcc*(1.0/(*m_vCovarianceDiag)(j)));
			}		
			dDeterminant = pow(dDeterminant,0.5);
			//fProbability = (float)max((double)LOG_LIKELIHOOD_FLOOR,(-0.5*fExponent)-log(fConstant*dDeterminant));	
			// we do not floor the probability since the contribution of each Gaussian component to the 
			// probability of the mixture (as computed in other methods) is not floored
			fProbability = (float)((-0.5*fExponent)-log(dConstant*dDeterminant));	
			if (finite(fProbability) == 0) {
				fProbability = LOG_LIKELIHOOD_FLOOR;
			}
		
			return fProbability;
		}
			
		inline float evaluateFullCovariance(float *fFeatures, int iTime) {
		
			assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_FULL);
			
			double dAcc;
			float fProbability = 0.0;
			double dExponent = 0.0;
			float *fMean = m_vMean->getData();
			float *fCholeskyFactorRow = m_mCholeskyFactor->getData();
			
			// for each row in the Cholesky Factor
			for(int j = 0 ; j < m_iDim ; ++j) {
				dAcc = 0.0;
				for(int k = j ; k < m_iDim ; ++k) {
					dAcc += (fFeatures[k]-fMean[k])*(*fCholeskyFactorRow);
					++fCholeskyFactorRow;
				}
				dExponent += dAcc*dAcc;
				// is it possible to stop the computation based on a very high partial exponent?
			}
			
			fProbability = (float)(exp(-0.5*dExponent)/m_dConstant);	
			fProbability = (float)std::max<double>((double)LOG_LIKELIHOOD_FLOOR,log(fProbability));	
			if (finite(fProbability) == 0) {
				fProbability = LOG_LIKELIHOOD_FLOOR;
			}
			
			return fProbability;
		}

	public:
	
		float *m_fCholeskyUpper;

		// constructor
		Gaussian(int iDim, int iCovarianceType);

		// destructor
		~Gaussian();	
		
		// return the mean
		inline Vector<float> &mean() {
		
			return *m_vMean;
		}
		
		// return the covariance
		inline Vector<float> &covarianceDiag() {
		
			return *m_vCovarianceDiag;
		}
		
		// return the covariance
		inline SMatrix<float> &covarianceFull() {
		
			return *m_mCovarianceFull;
		}		
		
		// return the weight
		inline float &weight() {
		
			return m_fWeight;
		}
		
		// set the mean and covariance (maybe different dimensionality)
		inline void setMeanCov(Vector<float> &vMean, Vector<float> &vCovariance) {
		
			assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL);
			assert(vMean.getDim() == vCovariance.getDim());
			
			// reallocate if different dimensionality?	
			if (vMean.getDim() != (unsigned int)m_iDim) {
				m_iDim = vMean.getDim();
				delete m_vMean;
				m_vMean = new Vector<float>(vMean);
				delete m_vCovarianceDiag;
				m_vCovarianceDiag = new Vector<float>(vCovariance);	
			} else {
				m_vMean->copy(vMean);
				m_vCovarianceDiag->copy(vCovariance);	
			}
		}
		
		// return the Cholesky factor
		inline TMatrix<float> &choleskyFactor() {
		
			return *m_mCholeskyFactor;
		}
		
		// return the evaluation constant
		inline double evaluationConstant() {
		
			return m_dConstant;
		}
		
		// return the logarithm of the evaluation constant
		inline double evaluationLogConstant() {
		
			return m_dLogConstant;
		}
		
		inline void resetTimeStamp() {
		
			m_iTimestamp = -1;
		}
		
		// set merged
		inline void setMerged(bool bMerged) {
		
			m_bMerged = bMerged;
		}
		
		// remove a Gaussian component
		void removeGaussianComponent(int g);
		
		// add a Gaussian component
		void addGaussianComponent(Gaussian *gaussian);	
		
		// precompute the constant used to speed-up Gaussian evaluations
		void precomputeConstant();

		// evaluate the Gaussian distribution
		inline float evaluate(float *fFeatures, int iTime) {
		
			if (m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				return evaluateDiagonalCovariance(fFeatures,iTime);
			} else {
				assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_FULL);
				return evaluateFullCovariance(fFeatures,iTime);
			}
		}				
		
		// print Gaussian information
		void print(ostream &os);
		
		// return the covariance modelling type as a string
		inline static const char *getCovarianceModellingStringFormat(int iType) {
		
			if (iType == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				return COVARIANCE_MODELLING_TYPE_DIAGONAL_STR;
			} else {
				return COVARIANCE_MODELLING_TYPE_FULL_STR;
			}			
		}
		
		// return the covariance modelling type
		inline static int getCovarianceModellingType(const char *strType) {
		
			if (strcmp(strType,COVARIANCE_MODELLING_TYPE_DIAGONAL_STR) == 0) {
				return COVARIANCE_MODELLING_TYPE_DIAGONAL;
			} else {
				assert(strcmp(strType,COVARIANCE_MODELLING_TYPE_FULL_STR) == 0);
				return COVARIANCE_MODELLING_TYPE_FULL;
			}
		}
		
		// return the weight for mixture increment
		inline float &weightMixtureIncrement() {
				
			return m_fWeightMixtureIncrement;
		}
		
};

typedef vector<Gaussian*> VGaussian;
typedef list<Gaussian*> LGaussian;

};	// end-of-namespace

#endif

