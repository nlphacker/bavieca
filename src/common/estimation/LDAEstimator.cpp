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


#include "HMMManager.h"
#include "LDAEstimator.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "Transform.h"

namespace Bavieca {

// constructor
LDAEstimator::LDAEstimator(int iDimSource, int iDimReduced, int iClasses)
{
	m_iClasses = iClasses;
	m_iDimSource = iDimSource;
	m_iDimReduced = iDimReduced;
	
	m_stats = NULL;
	m_smObsSquare = NULL;	
}

// destructor
LDAEstimator::~LDAEstimator()
{
	if (m_stats) {
		delete [] m_stats;
	}
	if (m_smObsSquare) {
		delete m_smObsSquare;
	}
}

// initialization
void LDAEstimator::initialize() {

	assert(m_iClasses > 1);

	// create an accumulator for each class
	for(int i=0 ; i < m_iClasses ; ++i) {
		m_stats = new ClassStats[m_iClasses];
	}
	m_smObsSquare = new SMatrix<double>(m_iDimSource);
}

// accumulate data
void LDAEstimator::accumulate(Vector<double> &vObs, SMatrix<double> &smObsSquare, 
	double dOcc, int iClass) {
/*
	assert(vObs.getDim() == m_iDimSource);
	assert(vObs.getDim() == smObsSquare.getDim());
	assert((iClass < m_iClasses) && (iClass >= 0));
	assert(dOcc >= 0.0);
	
	m_stats[iClass].dOcc += dOcc;
	m_stats[iClass].vObs.add(vObs);
	m_smObsSquare->add(smObsSquare);*/
}

// estimate the transform
void LDAEstimator::estimate() {

	// (1) compute global covariance (global cov = between class cov + within class cov)
/*	
	// compute global mean
	Vector<double> vMean(m_iDimSource);
	double dOcc = 0.0;
	for(int i=0 ; i < m_iClasses ; ++i) {
		dOcc += m_stats[i].dOcc;
		vMean.add(m_stats[i].vObs);
	}
	vMean.scale(1.0/dOcc);
	// compute global covariance
	SMatrix<double> smGlobalCov(*m_smObsSquare);
	smGlobalCov.scale(1.0/dOcc);
	smGlobalCov.addSquare(-1.0,vMean);	
	
	// (2) compute between-class covariance
	SMatrix<double> smBClassCov(m_iDimSource);
	Vector<double> vClassMean(m_iDimSource);
	for(int i=0 ; i < m_iClasses ; ++i) {	
		vClassMean.add(m_stats[i].vObs);
		vClassMean.scale(1.0/m_stats[i].dOcc);
		smBClassCov.addSquare(1.0,vClassMean);
	}
	smBClassCov.addSquare(-1.0,vMean);	
	
	// (3) compute within-class covariance from the already computed "covariances"
	SMatrix<double> smWClassCov(smGlobalCov);
	smWClassCov.add(-1.0,smBClassCov);
	TMatrix<double> tmWClassCovSqrt(m_iDimSource);
	//tmWClassCovSqrt.cholesky(smWClassCov);
	Matrix<double> mWClassCovSqrt(tmWClassCovSqrt);
	mWClassCovSqrt.invert();
	
	// (4) get the eigenvalues
	Matrix<double> mBClassCov(smBClassCov);
	Matrix<double> mSVD(m_iDimSource);
	//mSVD.addMMM(1.0,mWClassCovSqrt,kNoTrans,mBClassCov,kNoTrans,mWClassCovSqrt,kTrans,0.0);
	Vector<double> vEigenvalues(m_iDimSource);
	Matrix<double> mU(m_iDimSource),mV(m_iDimSource);
	mSVD.svd(vEigenvalues,mU,mV);	

	*/
}

};	// end-of-namespace


