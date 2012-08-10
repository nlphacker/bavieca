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

#include "MAPEstimator.h"

// constructor
MAPEstimator::MAPEstimator(HMMManager *hmmManager) {

	m_hmmManager = hmmManager;
	m_iFeatureDimensionality = hmmManager->getFeatureDimensionality();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iFeatureDimensionality,m_iCovarianceModeling);
	m_iHMMStates = -1;
	m_hmmStates = hmmManager->getHMMStates(&m_iHMMStates);	
}

// destructor
MAPEstimator::~MAPEstimator() {


}

// estimate the HMM parameters
bool MAPEstimator::estimateParameters(MAccumulatorPhysical &mAccumulator, float fPriorKnowledgeWeight) {
	
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
			estimateParameters(m_hmmStates[i],accumulators,fPriorKnowledgeWeight);
		}
		delete [] accumulators;
	}

	return true;
}

// estimate the parameters of the given HMM-state
bool MAPEstimator::estimateParameters(HMMState *hmmState, Accumulator **accumulators, float fPriorKnowledgeWeight) {

	int iFeatureDimensionality = hmmState->getFeatureDimensionality();
	int iCovarianceModeling = hmmState->getCovarianceModelling();	

	// (1) update the mean of each Gaussian component
	for(int g = 0 ; g < hmmState->getGaussianComponents() ; ++g) {
	
		Gaussian *gaussian = hmmState->getGaussian(g);
		Accumulator *accumulator = accumulators[g];
		// is there occupation for this component?
		if (accumulator == NULL) {
			continue;
		}
		
		// mean (needs to be computed first in the case of full covariance)
		for(int i = 0 ; i < iFeatureDimensionality ; ++i) {	
			gaussian->fMean[i] = (gaussian->fMean[i]*fPriorKnowledgeWeight+ accumulator->getObservation()[i])/(fPriorKnowledgeWeight+accumulator->getOccupation());
		}
	}	

	assert(hmmState->getGaussianComponents() > 0);
	
	return true;
}

