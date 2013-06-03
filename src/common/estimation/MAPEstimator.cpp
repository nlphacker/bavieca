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


#include "Accumulator.h"
#include "MAPEstimator.h"

namespace Bavieca {

// constructor
MAPEstimator::MAPEstimator(HMMManager *hmmManager) {

	m_hmmManager = hmmManager;
	m_iFeatureDimensionality = hmmManager->getFeatureDim();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iFeatureDimensionality,m_iCovarianceModeling);
	m_iHMMStates = -1;
	m_hmmStates = hmmManager->getHMMStates(&m_iHMMStates);	
}

// destructor
MAPEstimator::~MAPEstimator() {
}

// estimate the HMM parameters
void MAPEstimator::estimateParameters(MAccumulatorPhysical &mAccumulator, float fPriorKnowledgeWeight) {
	
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
			estimateParameters(m_hmmStates[i],accumulators,fPriorKnowledgeWeight);
		}
		delete [] accumulators;
	}
}

// estimate the parameters of the given HMM-state
void MAPEstimator::estimateParameters(HMMState *hmmState, Accumulator **accumulators, float fPriorKnowledgeWeight) {

	// (1) update the mean of each Gaussian component
	for(unsigned int g = 0 ; g < hmmState->getMixture().getNumberComponents() ; ++g) {
	
		Gaussian *gaussian = hmmState->getMixture()(g);
		Accumulator *accumulator = accumulators[g];
		// if there occupation for this component?
		if (accumulator == NULL) {
			continue;
		}
		
		// mean (needs to be computed first in the case of full covariance)
		gaussian->mean().mul(fPriorKnowledgeWeight);
		gaussian->mean().add(accumulator->getObservation());
		gaussian->mean().mul((float)(1.0f/(fPriorKnowledgeWeight+accumulator->getOccupation())));
	}	

	assert(hmmState->getMixture().getNumberComponents() > 0);
}

};	// end-of-namespace
