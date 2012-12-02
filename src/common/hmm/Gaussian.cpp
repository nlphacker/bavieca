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


#include "Gaussian.h"
#include "Global.h"
#include "LogMessage.h"

#include <stdlib.h>

namespace Bavieca {

Gaussian::Gaussian(int iDim, int iCovarianceType)
{
	m_iDim = iDim;
	m_iCovarianceType = iCovarianceType;
	m_iTimestamp = -1;
	m_fProbabilityCached = -FLT_MAX;
	m_dConstant = -FLT_MAX;
	m_dLogConstant = -FLT_MAX;	
	m_bMerged = false;

	// memory allocation
	m_vMean = new Vector<float>(m_iDim);
	if (m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_vCovarianceDiag = new Vector<float>(m_iDim);
		m_mCovarianceFull = NULL;
		m_mCholeskyFactor = NULL;
	} else {
		m_vCovarianceDiag = NULL;
		m_mCovarianceFull = new SMatrix<float>(m_iDim);
		m_mCholeskyFactor = new TMatrix<float>(m_iDim);
	}
}


Gaussian::~Gaussian()
{
	delete m_vMean;
	if (m_vCovarianceDiag) {
		delete m_vCovarianceDiag;
	}
	if (m_mCovarianceFull) {
		delete m_mCovarianceFull;
	}
	if (m_mCholeskyFactor) {
		delete m_mCholeskyFactor;
	}
}

// precompute the constant used to speed-up Gaussian evaluations
void Gaussian::precomputeConstant() {

	// diagonal covariance
	if (m_iCovarianceType == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
	
		double dDeterminant;
		
		// compute the determinant
		dDeterminant = 1.0;
		for(int i = 0 ; i < m_iDim ; ++i) {
			dDeterminant *= (*m_vCovarianceDiag)(i);
		}
		// compute the constant
		m_dConstant = sqrt(dDeterminant)*pow(PI_NUMBER*2.0,m_iDim/2.0);
		// for single Gaussian: the constant is in log mode
		m_dLogConstant = log(m_dConstant);
		assert(finite(m_dConstant) != 0);
	}
	// full covariance
	else {
		assert(m_iCovarianceType == COVARIANCE_MODELLING_TYPE_FULL);
		
		// compute the Cholesky factor of the inverse covariance matrix
		SMatrix<float> mCovarianceInverse(*m_mCovarianceFull);
		mCovarianceInverse.invert();
		TMatrix<float> mCholesky(m_iDim);
		mCholesky.choleskyDecomposition(mCovarianceInverse);
		//mCholesky.print();
		
		//compute the determinant of the inverse covariance matrix from Cholesky factor
		double dDetCholesky = mCholesky.determinant();
		double dDetInverse = dDetCholesky*dDetCholesky;
		
		// compute the evaluation constant
		m_dConstant = sqrt(1.0/dDetInverse)*pow(PI_NUMBER*2.0,m_iDim/2.0);
		
		// keep the Cholesky factor as an upper-triangular matrix, it is
		// needed for likelihood computation
		m_fCholeskyUpper = new float[(m_iDim*(m_iDim+1))/2];
		int iIndex = 0;
		for(int i=0 ; i < m_iDim ; ++i) {
			for(int j=i ; j < m_iDim ; ++j) {
				m_fCholeskyUpper[iIndex++] = mCholesky(j,i); 
			}	
		}	
	}	
}


// print Gaussian information
void Gaussian::print(ostream &os) {

}

};	// end-of-namespace





