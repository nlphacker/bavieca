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


#ifndef MAPESTIMATOR_H
#define MAPESTIMATOR_H

#include "HMMManager.h"

namespace Bavieca {

class Accumulator;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class MAPEstimator {

	private:	
	
		HMMManager *m_hmmManager;
		int m_iFeatureDimensionality;
		int m_iCovarianceModeling;
		int m_iCovarianceElements;
		int m_iHMMStates;
		HMMState **m_hmmStates;

		// estimate the parameters of the given HMM-state
		void estimateParameters(HMMState *hmmState, Accumulator **accumulators, float fPriorKnowledgeWeight);
	
	public:

		// constructor
		MAPEstimator(HMMManager *hmmManager);

		// destructor
		~MAPEstimator();
		
		// estimate the HMM parameters
		void estimateParameters(MAccumulatorPhysical &mAccumulator, float fPriorKnowledgeWeight);

};

};	// end-of-namespace

#endif
