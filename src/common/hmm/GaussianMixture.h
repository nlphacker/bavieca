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


#ifndef GAUSSIANMIXTURE_H
#define GAUSSIANMIXTURE_H

#include "Gaussian.h"

namespace Bavieca {

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class GaussianMixture {

	private:
	
		int m_iDim;									// feature dimensionality
		int m_iCovarianceType;					// covariance modelling type
		int m_iTimestamp;
		float m_fProbabilityCached;
		VGaussian m_vGaussian;			// Gaussian components
		
		// mixture evaluation (single Gaussian distribution)
		inline float evaluateDiagonalCovarianceSingleGaussian(float *fFeatures, int iTime)  {
		
			if ((iTime == m_iTimestamp) && (iTime != -1)) {	
				return m_fProbabilityCached;
			}
			
			Gaussian *gaussian = m_vGaussian[0];	
			
			double dAcc;
			double dExponent = 0.0;	
			float *fMean = gaussian->mean().getData();
			float *fCovariance = gaussian->covarianceDiag().getData();
			
			// partial distance elimination
			float fThresholdPDE = (float)((gaussian->evaluationLogConstant()+LOG_LIKELIHOOD_FLOOR)/-0.5f);
			
			for(int j = 0 ; j < m_iDim ; ++j) {
				dAcc = fFeatures[j]-fMean[j];
				dExponent += (dAcc*dAcc)/fCovariance[j];	
				if (dExponent >= fThresholdPDE) {
					return LOG_LIKELIHOOD_FLOOR;
				}
			}
			// for the single Gaussian case the constant is in log mode
			float fProbability = (float)((-0.5f*dExponent)-gaussian->evaluationLogConstant());	
			fProbability = (float)max((float)LOG_LIKELIHOOD_FLOOR,fProbability);	
			if (finite(fProbability) == 0) {
				fProbability = LOG_LIKELIHOOD_FLOOR;
			}
			
			// cache the probability
			m_iTimestamp = iTime;
			m_fProbabilityCached = fProbability;
			
			return fProbability;
		}
		
		// mixture evaluation
		inline float evaluateDiagonalCovariance(float *fFeatures, int iTime)  {
		
			if ((iTime == m_iTimestamp) && (iTime != -1)) {	
				return m_fProbabilityCached;
			}
			
			assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL);
			
			double dAcc;
			float fProbability = 0.0;
			float fExponent;	
			float *fMean;
			float *fCovariance;
			
			// accumulate emission probabilities from each gaussian
			for(unsigned int i = 0 ; i < m_vGaussian.size() ; ++i) {	
			
				fExponent = 0.0;
				fMean = m_vGaussian[i]->mean().getData();
				fCovariance = m_vGaussian[i]->covarianceDiag().getData();
				
				for(int j = 0 ; j < m_iDim ; ++j) {
					dAcc = fFeatures[j]-fMean[j];
					fExponent += (float)((dAcc*dAcc)/fCovariance[j]);
				}		
				fProbability += (float)(m_vGaussian[i]->weight()*exp(-0.5*fExponent)/m_vGaussian[i]->evaluationConstant());
			}
			
			fProbability = (float)std::max<double>((double)LOG_LIKELIHOOD_FLOOR,log(fProbability));	
			if (finite(fProbability) == 0) {
				fProbability = LOG_LIKELIHOOD_FLOOR;
			}
			
			// cache the probability
			m_iTimestamp = iTime;
			m_fProbabilityCached = fProbability;
			
			return fProbability;
		}
		
		// mixture evaluation
		// note: it performs a Cholesky decomposition of the covariance matrix (which is symmetric) in order
		//       to express it as the product of two triangular matrices, which allow a more efficient computation
		inline float evaluateFullCovariance(float *fFeatures, int iTime)  {
		
			if ((iTime == m_iTimestamp) && (iTime != -1)) {	
				return m_fProbabilityCached;
			}	
			
			assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_FULL);
			
			double dAcc;
			float fProbability = 0.0;
			float fExponent;	
			float *fMean;
			float *fCholeskyFactorRow;
			
			// accumulate emission probabilities from each Gaussian component
			for(unsigned int i = 0 ; i < m_vGaussian.size() ; ++i) {	
			
				fExponent = 0.0;
				fMean = m_vGaussian[i]->mean().getData();
				//fCholeskyFactorRow = m_vGaussian[i]->choleskyFactor().getData();
				fCholeskyFactorRow = m_vGaussian[i]->m_fCholeskyUpper;
				
				// for each row in the Cholesky Factor
				for(int j = 0 ; j < m_iDim ; ++j) {
					dAcc = 0.0;
					for(int k = j ; k < m_iDim ; ++k) {
						dAcc += (fFeatures[k]-fMean[k])*(*fCholeskyFactorRow);
						++fCholeskyFactorRow;
					}
					fExponent += (float)(dAcc*dAcc);
					// is it possible to stop the computation based on a very high partial exponent?
				}
							
				fProbability += (float)(m_vGaussian[i]->weight()*exp(-0.5*fExponent)/m_vGaussian[i]->evaluationConstant());	
			}	
			
			fProbability = (float)std::max<double>((double)LOG_LIKELIHOOD_FLOOR,log(fProbability));	
			if (finite(fProbability) == 0) {
				fProbability = LOG_LIKELIHOOD_FLOOR;
			}
			
			// cache the probability
			m_iTimestamp = iTime;
			m_fProbabilityCached = fProbability;
			
			return fProbability;
		}
	

	public:

		GaussianMixture(int iDim, int iCovarianceType, int iComponents);

		~GaussianMixture();
		
		// returns the given element as a r-value
		inline Gaussian *operator() (unsigned int i) {
		
			assert(i < m_vGaussian.size());
			return m_vGaussian[i];
		}
		
		inline void resetTimeStamp() {
		
			m_iTimestamp = -1;
			for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
				(*it)->resetTimeStamp();
			}
		}
		
		inline void setMerged(bool bMerged) {
		
			for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
				(*it)->setMerged(bMerged);
			}
		}
		
		inline unsigned int getNumberComponents() {
		
			return (unsigned int)m_vGaussian.size();
		}
		
		// return the number of free parameters associated to this HMM-state
		inline unsigned int getNumberFreeParameters() {
		
			// digonal covariance
			if (m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				return (unsigned int)(m_vGaussian.size()*(2*m_iDim+1));	
			} 
			// full covariance
			else {
				assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_FULL);
				return (unsigned int)(m_vGaussian.size()*(m_iDim+((m_iDim*(m_iDim+1))/2)+1));
			}	
		}
		
		// set a new mean
		inline void setMeanAllComponents(const VectorBase<float> &vMean) {		
			
			for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
				(*it)->mean().copy(vMean);
			}
		}
		
		// set a new covariance
		inline void setCovarianceAllComponents(const VectorBase<float> &vCovarianceDiag) {
			
			for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
				(*it)->covarianceDiag().copy(vCovarianceDiag);
			}
		}
		
		// set a new covariance
		inline void setCovarianceAllComponents(SMatrix<float> &mCovarianceFull) {
			
			for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
				(*it)->covarianceFull().copy(mCovarianceFull);
			}
		}		
		
		// remove a Gaussian component
		void removeGaussianComponent(int g);
		
		// add a Gaussian component
		void addGaussianComponent(Gaussian *gaussian);
		
		// precompute evaluation constants
		void precomputeConstant();
			
		// compute the emission probability of a state given the feature vector
		inline float evaluate(float *fFeatures, int iTime)  {
			
			// diagonal covariance
			if (m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				if (m_vGaussian.size() == 1) {
					return evaluateDiagonalCovarianceSingleGaussian(fFeatures,iTime);
				} else {
					return evaluateDiagonalCovariance(fFeatures,iTime);
				}
			} 
			// full covariance
			else {
				return evaluateFullCovariance(fFeatures,iTime);
			}		
		}
};

};	// end-of-namespace

#endif
