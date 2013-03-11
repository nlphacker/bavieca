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


#include "MLEstimator.h"

namespace Bavieca {

// constructor
MLEstimator::MLEstimator(HMMManager *hmmManager) {

	m_hmmManager = hmmManager;
	m_iDim = hmmManager->getFeatureDimensionality();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iDim,m_iCovarianceModeling);
	m_iHMMStates = -1;
	m_hmmStates = hmmManager->getHMMStates(&m_iHMMStates);	
}

// destructor
MLEstimator::~MLEstimator() {


}

// estimate the HMM parameters
void MLEstimator::estimateParameters(MAccumulatorPhysical &mAccumulator, 
	bool bUpdateCovariance) {
	
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		bool bData = false;
		int iComponents = m_hmmStates[i]->getMixture().getNumberComponents();
		Accumulator **accumulators = new Accumulator*[iComponents];
		for(int g=0 ; g < iComponents ; ++g) {
			MAccumulatorPhysical::iterator it = mAccumulator.find(Accumulator::getPhysicalAccumulatorKey(i,g));
			if (it != mAccumulator.end()) {
				accumulators[g] = it->second;
				bData = true;
			} else {
				accumulators[g] = NULL;
			}
		}
		// is there data to update the HMM-state
		if (bData) {
			estimateParameters(m_hmmStates[i],accumulators,bUpdateCovariance);
		}
		delete [] accumulators;
	}
}


// estimate the parameters of the given HMM-state
void MLEstimator::estimateParameters(HMMState *hmmState, Accumulator **accumulators, 
	bool bUpdateCovariance) {

	int iCovarianceModeling = hmmState->getCovarianceModelling();

	// (1) update the mean and the covariance of each Gaussian component
	for(unsigned int g = 0 ; g < hmmState->getMixture().getNumberComponents() ; ++g) {
	
		Gaussian *gaussian = hmmState->getMixture()(g);
		Accumulator *accumulator = accumulators[g];
		// is there occupation for this component?
		if (accumulator == NULL) {
			// remove the Gaussian if there are other Gaussian components
			if (hmmState->getMixture().getNumberComponents() > 1) {
				hmmState->getMixture().removeGaussianComponent(g);
				g = -1;
				continue;
			} 
			// otherwise end the parameter estimation of this HMM-state
			else {
				return;
			}
		}
		
		assert(accumulator->getOccupation() > 0.0);
		
		// mean (needs to be computed first in the case of full covariance)
		gaussian->mean().mul((float)(1.0/accumulator->getOccupation()),accumulator->getObservation());
		
		// covariance
		if (bUpdateCovariance) {
			// diagonal
			if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				gaussian->covarianceDiag().mul((float)(1.0/accumulator->getOccupation()),accumulator->getObservationSquareDiag());
				gaussian->covarianceDiag().addSquare(-1.0,gaussian->mean());
			// full
			} else {
				assert(iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
				gaussian->covarianceFull().mul((float)(1.0/accumulator->getOccupation()),accumulator->getObservationSquareFull());
				gaussian->covarianceFull().addSquare(-1.0,gaussian->mean());
			}	
		}
	}
	
	// (2) compute the component weights
	
	// multiple Gaussian components
	int iComponents = hmmState->getMixture().getNumberComponents();
	if (iComponents > 1) {
		// compute state occupation
		double dOccupationState = 0.0;
		for(int g = 0 ; g < iComponents ; ++g) {
			dOccupationState += accumulators[g]->getOccupation();
		}	
		// update weights: ratio between mixture occupation and state occupation
		for(int g = 0 ; g < iComponents ; ++g) {
			hmmState->getMixture()(g)->weight() = (float)(accumulators[g]->getOccupation()/dOccupationState);
		}	
	} 
	// single Gaussian component
	else {
		assert(iComponents == 1);	
		hmmState->getMixture()(0)->weight() = 1.0;
	}
	
	assert(iComponents > 0);
}

// compute the covariance floor
void MLEstimator::computeCovarianceFloor(HMMManager *hmmManager, MAccumulatorPhysical &mAccumulator, 
	float fScalingFactor, Vector<double> &vCovarianceFloor) {
	
	assert(hmmManager->getFeatureDimensionality() == vCovarianceFloor.getDim());
	assert((fScalingFactor > 0.0) && (fScalingFactor < 1.0));
	int iHMMStates = -1;
	HMMState **hmmStates = hmmManager->getHMMStates(&iHMMStates);

	// allocate memory for the covariance floor
	vCovarianceFloor.zero();
	double dOccupation = 0.0;
	SMatrix<double> mAcc(vCovarianceFloor.getDim());	
	
	// compute the covariance used as reference (the covariance of each mixture is weighted according to its occupation)
	for(int i=0 ; i < iHMMStates ; ++i) {	
		for(unsigned int g=0 ; g < hmmStates[i]->getMixture().getNumberComponents() ; ++g) {
			MAccumulatorPhysical::iterator jt = mAccumulator.find(Accumulator::getPhysicalAccumulatorKey(i,g));
			if (jt == mAccumulator.end()) {
				continue;
			}
			Accumulator *accumulator = jt->second;
			assert(accumulator != NULL);
			if (accumulator->getOccupation() > 0.0) {
				if (hmmManager->getCovarianceModelling() == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
					vCovarianceFloor.add(accumulator->getOccupation(),hmmStates[i]->getMixture()(g)->covarianceDiag());
				} else {
					assert(hmmManager->getCovarianceModelling() == COVARIANCE_MODELLING_TYPE_FULL);	
					mAcc.add(accumulator->getOccupation(),hmmStates[i]->getMixture()(g)->covarianceFull());
				}
				dOccupation += accumulator->getOccupation();
			}
		}
	}
	
	if (hmmManager->getCovarianceModelling() == COVARIANCE_MODELLING_TYPE_FULL) {
		mAcc.getDiagonal(vCovarianceFloor);
	}
		
	// divide by the total occupation and multiply by the scaling factor
	vCovarianceFloor.mul(fScalingFactor/dOccupation);
}

// set a floot to each of the HMMs covariances
// note: only the diagonal elements of a full covariance matrix are floored
void MLEstimator::floorCovariances(MAccumulatorPhysical &mAccumulator, float fScalingFactor) {

	assert((fScalingFactor > 0.0) && (fScalingFactor < 1.0));	
	
	Vector<double> vCovarianceFloor(m_iDim);
	computeCovarianceFloor(m_hmmManager,mAccumulator,fScalingFactor,vCovarianceFloor);		
	applyCovarianceFloor(m_hmmManager,vCovarianceFloor);
}

// set a floot to each of the HMMs covariances
// note: only the diagonal elements of a full covariance matrix are floored
void MLEstimator::floorCovariances(HMMManager *hmmManager, MAccumulatorPhysical &mAccumulator, 
	float fScalingFactor) {

	assert((fScalingFactor > 0.0) && (fScalingFactor < 1.0));	
	
	Vector<double> vCovarianceFloor(hmmManager->getFeatureDimensionality());
	computeCovarianceFloor(hmmManager,mAccumulator,fScalingFactor,vCovarianceFloor);
	applyCovarianceFloor(hmmManager,vCovarianceFloor);
}

// do the actual covariance flooring
void MLEstimator::applyCovarianceFloor(HMMManager *hmmManager, Vector<double> &vCovarianceFloor) {

	int iHMMStates = -1;
	HMMState **hmmStates = hmmManager->getHMMStates(&iHMMStates);
	for(int i=0 ; i < iHMMStates ; ++i) {	
		HMMState *hmmState = hmmStates[i];
		int iGaussianComponents = hmmState->getMixture().getNumberComponents();	
		for(int m = 0 ; m < iGaussianComponents ; ++m) {
			Gaussian *gaussian = hmmState->getMixture()(m);	
			// diagonal
			if (hmmManager->getCovarianceModelling() == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				for(int i = 0 ; i < vCovarianceFloor.getDim() ; ++i) {
					if (gaussian->covarianceDiag()(i) < vCovarianceFloor(i)) {
						gaussian->covarianceDiag()(i) = (float)vCovarianceFloor(i);
					}
				}
			}
			// full
			else {
				assert(hmmManager->getCovarianceModelling() == COVARIANCE_MODELLING_TYPE_FULL);
				for(int i = 0 ; i < vCovarianceFloor.getDim() ; ++i) {
					if (gaussian->covarianceFull()(i,i) < vCovarianceFloor(i)) {
						gaussian->covarianceFull()(i,i) = (float)vCovarianceFloor(i);
					}
				}
			}
		}	
	}
}

};	// end-of-namespace
