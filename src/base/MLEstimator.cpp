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

#include "MLEstimator.h"

// constructor
MLEstimator::MLEstimator(HMMManager *hmmManager) {

	m_hmmManager = hmmManager;
	m_iFeatureDimensionality = hmmManager->getFeatureDimensionality();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iFeatureDimensionality,m_iCovarianceModeling);
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
		int iComponents = m_hmmStates[i]->getGaussianComponents();
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
		if (bData == true) {
			estimateParameters(m_hmmStates[i],accumulators,bUpdateCovariance);
		}
		delete [] accumulators;
	}
}


// estimate the parameters of the given HMM-state
void MLEstimator::estimateParameters(HMMState *hmmState, Accumulator **accumulators, 
	bool bUpdateCovariance) {

	int iFeatureDimensionality = hmmState->getFeatureDimensionality();
	int iCovarianceModeling = hmmState->getCovarianceModelling();

	// (1) update the mean and the covariance of each Gaussian component
	for(int g = 0 ; g < hmmState->getGaussianComponents() ; ++g) {
	
		Gaussian *gaussian = hmmState->getGaussian(g);
		Accumulator *accumulator = accumulators[g];
		// is there occupation for this component?
		if (accumulator == NULL) {
			// remove the Gaussian if there are other Gaussian components
			if (hmmState->getGaussianComponents() > 1) {
				hmmState->removeGaussianComponent(g);
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
		for(int i = 0 ; i < iFeatureDimensionality ; ++i) {
			gaussian->fMean[i] = (float)(accumulator->getObservation()[i]/accumulator->getOccupation());
		}
		
		// covariance
		if (bUpdateCovariance == true) {
			// diagonal
			if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				for(int i = 0 ; i < iFeatureDimensionality ; ++i) {
					gaussian->fCovariance[i] = (accumulator->getObservationSquare()[i]/accumulator->getOccupation())-(gaussian->fMean[i]*gaussian->fMean[i]);
				}
			}
			// full
			else {
				assert(iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);	
				int iIndex = 0;
				for(int i = 0 ; i < iFeatureDimensionality ; ++i) {
					for(int j = i ; j < iFeatureDimensionality ; ++j, ++iIndex) {
						gaussian->fCovariance[iIndex] = (accumulator->getObservationSquare()[iIndex]/accumulator->getOccupation())-
							(gaussian->fMean[i]*gaussian->fMean[j]);
						if (i == j) {
							assert(gaussian->fCovariance[iIndex] > 0.0);
						}
					}
				}
			}	
		}
	}
	
	// (2) compute the component weights
	
	// multiple Gaussian components
	if (hmmState->getGaussianComponents() > 1) {
		// compute state occupation
		double dOccupationState = 0.0;
		for(int g = 0 ; g < hmmState->getGaussianComponents() ; ++g) {
			dOccupationState += accumulators[g]->getOccupation();
		}	
		// update weights: ratio between mixture occupation and state occupation
		for(int g = 0 ; g < hmmState->getGaussianComponents() ; ++g) {
			hmmState->getGaussian(g)->fWeight = (float)(accumulators[g]->getOccupation()/dOccupationState);
		}	
	} 
	// single Gaussian component
	else {
		assert(hmmState->getGaussianComponents() == 1);	
		hmmState->getGaussian(0)->fWeight = 1.0;
	}
	
	assert(hmmState->getGaussianComponents() > 0);
}

// compute the covariance floor
double *MLEstimator::computeCovarianceFloor(HMMManager *hmmManager, MAccumulatorPhysical &mAccumulator, float fScalingFactor) {

	int iHMMStates = -1;
	HMMState **hmmStates = hmmManager->getHMMStates(&iHMMStates);
	int iFeatureDimensionality = hmmManager->getFeatureDimensionality();
	int iCovarianceElements = hmmManager->getCovarianceElements();

	// allocate memory for the covariance floor
	double *dCovarianceFloor = new double[iCovarianceElements];
	for(int i=0 ; i < iCovarianceElements ; ++i) {
		dCovarianceFloor[i] = 0.0;
	}
	double dOccupation = 0.0;
	
	// compute the covariance used as reference (the covariance of each mixture is weighted according to its occupation)
	for(int i=0 ; i < iHMMStates ; ++i) {
		if (hmmStates[i]->isOccupied() == true) {
			int g = 0;
			for(VGaussian::iterator it = hmmStates[i]->m_vGaussian.begin() ; it != hmmStates[i]->m_vGaussian.end() ;
				++it, ++g) {
				// get the Gaussian accumulator
				MAccumulatorPhysical::iterator jt = mAccumulator.find(Accumulator::getPhysicalAccumulatorKey(i,g));
				if (jt == mAccumulator.end()) {
					continue;
				}
				Accumulator *accumulator = jt->second;
				assert(accumulator != NULL);
				for(int j=0 ; j < iCovarianceElements ; ++j) {
					dCovarianceFloor[j] += accumulator->getOccupation()*(*it)->fCovariance[j];
				}
				dOccupation += accumulator->getOccupation();
			}
		}
	}
		
	// divide by the total occupation and multiply by the scaling factor
	for(int i=0 ; i < iCovarianceElements ; ++i) {
		dCovarianceFloor[i] = (fScalingFactor*dCovarianceFloor[i])/dOccupation;
	}

	return dCovarianceFloor;
}

// set a floot to each of the HMMs covariances
// note: only the diagonal elements of a full covariance matrix are floored
void MLEstimator::floorCovariances(MAccumulatorPhysical &mAccumulator, float fScalingFactor) {

	assert((fScalingFactor > 0.0) && (fScalingFactor < 1.0));	
	
	// (1) compute the covariance floor	
	
	// allocate memory for the covariance floor
	float *fCovarianceFloor = new float[m_iCovarianceElements];
	double *dCovarianceFloor = new double[m_iCovarianceElements];
	for(int i=0 ; i < m_iCovarianceElements ; ++i) {
		dCovarianceFloor[i] = 0.0;
	}
	double dOccupation = 0.0;
	
	// compute the covariance used as reference (the covariance of each mixture is weighted according to its occupation)
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		if (m_hmmStates[i]->isOccupied() == true) {
			int g = 0;
			for(VGaussian::iterator it = m_hmmStates[i]->m_vGaussian.begin() ; it != m_hmmStates[i]->m_vGaussian.end() ;
				++it, ++g) {
				// get the Gaussian accumulator
				MAccumulatorPhysical::iterator jt = mAccumulator.find(Accumulator::getPhysicalAccumulatorKey(i,g));
				if (jt == mAccumulator.end()) {
					continue;
				}
				Accumulator *accumulator = jt->second;
				assert(accumulator != NULL);
				for(int j=0 ; j < m_iCovarianceElements ; ++j) {
					dCovarianceFloor[j] += accumulator->getOccupation()*(*it)->fCovariance[j];
				}
				dOccupation += accumulator->getOccupation();
			}
		}
	}
		
	// divide by the total occupation and multiply by the scaling factor
	for(int i=0 ; i < m_iCovarianceElements ; ++i) {
		fCovarianceFloor[i] = (fScalingFactor*dCovarianceFloor[i])/dOccupation;
	}

	delete [] dCovarianceFloor;
	
	// (2) do the actual covariance flooring
	for(int i=0 ; i<m_iHMMStates ; ++i) {
		floorCovariances(m_iFeatureDimensionality,m_iCovarianceModeling,m_hmmStates[i],fCovarianceFloor);
	}
}

// set a floor to the covariance of all the Gaussian components
void MLEstimator::floorCovariances(int iFeatureDimensionality, int iCovarianceModeling, 
	HMMState *hmmState, float *fCovarianceFloor) {

	int iGaussianComponents = hmmState->getGaussianComponents();	

	// (1) set a floor to the covariance of each Gaussian component (diagonal/full)
	// diagonal-covariance
	if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		for(int m = 0 ; m < iGaussianComponents ; ++m) {
			Gaussian *gaussian = hmmState->getGaussian(m);	
			for(int i = 0 ; i < iFeatureDimensionality ; ++i) {
				if (gaussian->fCovariance[i] < fCovarianceFloor[i]) {
					gaussian->fCovariance[i] = fCovarianceFloor[i];
				}
			}
		}
	}
	// full-covariance 
	// - important: only the elements in the main diagonal need to be floored
	else {
		assert(iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
		for(int m = 0 ; m < iGaussianComponents ; ++m) {
			Gaussian *gaussian = hmmState->getGaussian(m);	
			for(int i = 0 ; i < iFeatureDimensionality ; ++i) {
				if (HMMState::getElement(gaussian->fCovariance,i,i,iFeatureDimensionality) < 
					HMMState::getElement(fCovarianceFloor,i,i,iFeatureDimensionality)) {
					HMMState::setElement(gaussian->fCovariance,i,i,iFeatureDimensionality,
						HMMState::getElement(fCovarianceFloor,i,i,iFeatureDimensionality));
				}
			}
		}
	}
}
