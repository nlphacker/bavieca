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


#ifndef DTESTIMATOR_H
#define DTESTIMATOR_H

using namespace std;

#include <string>

#include "HMMManager.h"

namespace Bavieca {

class ForwardBackward;
class LexiconManager;
class MLFFile;
class PhoneSet;

// I-smoothing
#define I_SMOOTHING_NONE								"none"
#define I_SMOOTHING_PREVIOUS_ITERATION				"prev"

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class DTEstimator {

	private:
	
		HMMManager *m_hmmManager;
		int m_iDim;
		int m_iCovarianceModeling;
		int m_iCovarianceElements;
		int m_iHMMStates;
		HMMState **m_hmmStates;
		
		// accumulators
		MAccumulatorPhysical m_mAccumulatorNum;
		MAccumulatorPhysical m_mAccumulatorDen;

		// estimate the parameters of the given HMM-state
		void estimateParameters(HMMState *hmmState, Accumulator **accumulatorsNum, 
			Accumulator **accumulatorsDen, float fE, const char *strISmoothingType, float fTau, bool bUpdateCovariance = true);
		
		// compute the Gaussian-specific learning-rate constant
		double computeLearningConstant(Gaussian *gaussian, Accumulator *accumulatorNum, 
			Accumulator *accumulatorDen, float fE, bool bISmoothingPreviousIteration = false, float fTau = 0.0);
			
		// return wether the given value of the learning-constant makes the covariance positive definite
		bool isPositiveDefinite(Gaussian *gaussian, Accumulator *accumulatorNum, 
			Accumulator *accumulatorDen, double dD);

	public:

		// constructor
		DTEstimator(HMMManager *hmmManager);

		// destructor
		~DTEstimator();
		
		// estimate the HMM parameters	
		void estimateParameters(const char *strFileAccListNum, const char *strFileAccListDen, 
			float fE, const char *strISmoothingType, float fTau, bool bUpdateCovariance = true);
			
		// set a floor to each of the HMMs covariances
		// note: only the diagonal elements of a full covariance matrix are floored
		void floorCovariances(float fScalingFactor);

};

};	// end-of-namespace

#endif
