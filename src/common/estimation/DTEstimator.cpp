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


#include "DTEstimator.h"
#include "FeatureFile.h"
#include "MLEstimator.h"

#include "ForwardBackward.h"
#include "LexiconManager.h"
#include "MLFFile.h"
#include "PhoneSet.h"

namespace Bavieca {

// constructor
DTEstimator::DTEstimator(HMMManager *hmmManager) {
	
	m_hmmManager = hmmManager;
	m_iDim = hmmManager->getFeatureDimensionality();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iDim,m_iCovarianceModeling);
	m_iHMMStates = -1;
	m_hmmStates = hmmManager->getHMMStates(&m_iHMMStates);	
}

// destructor
DTEstimator::~DTEstimator() {

	// destroy the accumulators
	Accumulator::destroy(m_mAccumulatorNum);
	Accumulator::destroy(m_mAccumulatorDen);	
}

// estimate the HMM parameters
void DTEstimator::estimateParameters(const char *strFileAccListNum, const char *strFileAccListDen, 
	float fE, const char *strISmoothingType, float fTau, bool bUpdateCovariance) {

	// (1) load numerator and denominator accumulators
	AccMetadata metadata;
	Accumulator::loadAccumulatorList(strFileAccListNum,m_mAccumulatorNum,metadata);
	Accumulator::loadAccumulatorList(strFileAccListDen,m_mAccumulatorDen,metadata);
	
	// (2) estimate the parameters
	
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		bool bDataNum = false;
		bool bDataDen = false;
		int iComponents = m_hmmStates[i]->getMixture().getNumberComponents();
		Accumulator **accumulatorsNum = new Accumulator*[iComponents];
		Accumulator **accumulatorsDen = new Accumulator*[iComponents];
		for(int g=0 ; g < iComponents ; ++g) {
			// numerator accumulator
			MAccumulatorPhysical::iterator it = m_mAccumulatorNum.find(Accumulator::getPhysicalAccumulatorKey(i,g));
			if (it != m_mAccumulatorNum.end()) {
				accumulatorsNum[g] = it->second;
				bDataNum = true;
			} else {
				accumulatorsNum[g] = NULL;
			}
			// denominator accumulator
			MAccumulatorPhysical::iterator jt = m_mAccumulatorDen.find(Accumulator::getPhysicalAccumulatorKey(i,g));
			if (jt != m_mAccumulatorDen.end()) {
				accumulatorsDen[g] = jt->second;
				bDataDen = true;
			} else {
				accumulatorsDen[g] = NULL;
			}
		}
		// if there is data to update the HMM-state
		if (bDataNum && bDataDen) {
			estimateParameters(m_hmmStates[i],accumulatorsNum,accumulatorsDen,fE,strISmoothingType,fTau,bUpdateCovariance);
		} else {
			BVC_WARNING << "no data to estimate the hmm-state";
		}
		delete [] accumulatorsNum;
		delete [] accumulatorsDen;
	}
}

// estimate the parameters of the given HMM-state
void DTEstimator::estimateParameters(HMMState *hmmState, Accumulator **accumulatorsNum, 
	Accumulator **accumulatorsDen, float fE, const char *strISmoothingType, float fTau, bool bUpdateCovariance) {
	
	// (1) update the mean and the covariance of each Gaussian component
	for(unsigned int g = 0 ; g < hmmState->getMixture().getNumberComponents() ; ++g) {
	
		Gaussian *gaussian = hmmState->getMixture()(g);
		Accumulator *accumulatorNum = accumulatorsNum[g];
		Accumulator *accumulatorDen = accumulatorsDen[g];
		assert(accumulatorNum && accumulatorDen);	
		assert(accumulatorNum->getOccupation() > 0.0);
		assert(accumulatorDen->getOccupation() > 0.0);
		
		// compute the Gaussian-specific learning-constant
		double dD;
		if (strcmp(strISmoothingType,I_SMOOTHING_PREVIOUS_ITERATION) == 0) {
			dD = computeLearningConstant(gaussian,accumulatorNum,accumulatorDen,fE,true,fTau);
		} else {
			dD = computeLearningConstant(gaussian,accumulatorNum,accumulatorDen,fE);
		}	
		
		// keep the original mean
		Vector<float> vMeanBase(gaussian->mean());
		
		// mean
		gaussian->mean().copy(accumulatorNum->getObservation());
		gaussian->mean().add(-1.0,accumulatorDen->getObservation());
		gaussian->mean().add((float)dD,vMeanBase);
		gaussian->mean().mul((float)(1.0/(accumulatorNum->getOccupation()-accumulatorDen->getOccupation()+dD)));	
		
		// covariance
		assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL);
		gaussian->covarianceDiag().addSquare(1.0,vMeanBase);
		gaussian->covarianceDiag().mul((float)dD);
		gaussian->covarianceDiag().add(accumulatorNum->getObservationSquareDiag());
		gaussian->covarianceDiag().add(-1.0,accumulatorDen->getObservationSquareDiag());	
		gaussian->covarianceDiag().addSquare((float)dD,vMeanBase);
		gaussian->covarianceDiag().mul((float)(1.0/(accumulatorNum->getOccupation()-accumulatorDen->getOccupation()+dD)));	
		gaussian->covarianceDiag().addSquare(-1.0,gaussian->mean());
	}
	
	// (2) compute the Gaussian priors
	
	// so far leave them unchanged...

}

// compute the Gaussian-specific learning-rate constant
double DTEstimator::computeLearningConstant(Gaussian *gaussian, Accumulator *accumulatorNum, 
	Accumulator *accumulatorDen, float fE, bool bISmoothingPreviousIteration, float fTau) {

	// case 1: double the Gaussian occupation
	double dCase1 = fE*accumulatorDen->getOccupation();
	if (bISmoothingPreviousIteration) {
		dCase1 += fTau;
	}
	
	// case 2: minimum value needed to make the covariance positive definite
	double dCase2 = dCase1/2.0;
	
	double dIncrement = dCase2*0.01;
	while(isPositiveDefinite(gaussian,accumulatorNum,accumulatorDen,dCase2) == false) {
		dCase2 += dIncrement;	
	}
	
	return dCase2*2.0;
}

// return wether the given value of the learning-constant makes the covariance positive definite
bool DTEstimator::isPositiveDefinite(Gaussian *gaussian, Accumulator *accumulatorNum, 
	Accumulator *accumulatorDen, double dD) {
	
	// keep the original mean
	Vector<float> vMeanBase(gaussian->mean());
	
	// mean
	Vector<float> vMean(m_iDim);
	vMean.copy(accumulatorNum->getObservation());
	vMean.add(-1.0,accumulatorDen->getObservation());
	vMean.add((float)dD,vMeanBase);
	vMean.mul((float)(1.0/(accumulatorNum->getOccupation()-accumulatorDen->getOccupation()+dD)));	

	// covariance
	assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL);
	
	for(int i = 0 ; i < m_iDim ; ++i) {
		double dCovarianceElement = ((accumulatorNum->getObservationSquareDiag()(i)
			-accumulatorDen->getObservationSquareDiag()(i))
			+dD*(gaussian->covarianceDiag()(i)+(vMeanBase(i)*vMeanBase(i))))/
			((accumulatorNum->getOccupation()-accumulatorDen->getOccupation())+dD);
		dCovarianceElement -= vMean(i)*vMean(i);
		if (dCovarianceElement <= 0.0) {
			return false;
		}
	}

	return true;
}

// covariance flooring
void DTEstimator::floorCovariances(float fScaleFactor) {
	
	MLEstimator::floorCovariances(m_hmmManager,m_mAccumulatorNum,fScaleFactor);
}

};	// end-of-namespace
