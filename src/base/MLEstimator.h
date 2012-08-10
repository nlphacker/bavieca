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

#ifndef MLESTIMATOR_H
#define MLESTIMATOR_H

#include "Accumulator.h"
#include "HMMManager.h"

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class MLEstimator {

	private:	
	
		HMMManager *m_hmmManager;
		int m_iFeatureDimensionality;
		int m_iCovarianceModeling;
		int m_iCovarianceElements;
		int m_iHMMStates;
		HMMState **m_hmmStates;

		// estimate the parameters of the given HMM-state
		void estimateParameters(HMMState *hmmState, Accumulator **accumulators, 
			bool bUpdateCovariance);
	
	
	public:

		// constructor
		MLEstimator(HMMManager *hmmManager);

		// destructor
		~MLEstimator();
		
		// set a floor to the covariance of all the Gaussian components
		static void floorCovariances(int iFeatureDimensionality, int iCovarianceModeling, 
			HMMState *hmmState, float *fCovarianceFloor);
		
		// estimate the HMM parameters
		void estimateParameters(MAccumulatorPhysical &mAccumulator, 
			bool bUpdateCovariance);
			
		// set a floot to each of the HMMs covariances
		// note: only the diagonal elements of a full covariance matrix are floored
		void floorCovariances(MAccumulatorPhysical &mAccumulator, float fScalingFactor);
		
		// compute the covariance floor
		static double *computeCovarianceFloor(HMMManager *hmmManager, 
			MAccumulatorPhysical &mAccumulator, float fScalingFactor);

};

#endif
