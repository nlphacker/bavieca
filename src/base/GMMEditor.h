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

#ifndef GMMEDITOR_H
#define GMMEDITOR_H

using namespace std;

#include <map>

#ifdef __linux__
using namespace __gnu_cxx;
#include <ext/hash_map>
#endif

#ifdef _WIN32
#include <hash_map>
#endif

#include "HMMManager.h"

typedef struct {
	double dOccupation;
	bool bMerged;
} GaussianInfo;

// structure to keep physical accumulators
typedef map<Gaussian*,GaussianInfo> MGaussianInfo;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class GMMEditor {

	private:
	
		HMMManager *m_hmmManager;
		MAccumulatorPhysical *m_mAccumulatorPhysical;
		MGaussianInfo *m_mGaussianInfo;
		int m_iFeatureDimensionality;
		int m_iCovarianceModeling;
		int m_iCovarianceElements;
		int m_iHMMStates;
		HMMState **m_hmmStates;	

		// merge Gaussian components with not enough training data associated
		// when a Gaussian is split into two gaussians, and the parameters of the newly created 
		// Gaussian components are reestimated (iteration), it is possible that one of the gaussians 
		// has very low occupation, in that case we need to merge that Gaussian to the closest one
		bool mixtureMerge(HMMState *hmmState, float fMinimumGaussianOccupation, float fMinimumGaussianWeight, double *dCovarianceFloor);
		
		// adds a Gaussian to the mixture 
		// note: there is no need to floor the covariance of the new components, since they stay the same
		int mixtureIncrement(HMMState *hmmState, int iIncrement, int iSplittingCriterion, 
			bool bAllowMultipleSplitsPerComponent, float fMinimumGaussianOccupation, float fSplittingEpsilon);	
		
		// double the number of Gaussian components in the mixture
		// there is no need to select which Gaussian components will be split, since all of them will (except those 
		// that were just merged or dont have enough occupation)
		int mixtureDouble(HMMState *hmmState, float fMinimumGaussianOccupation, float fSplittingEpsilon);
	
	public:
    
		// contructor
		GMMEditor(HMMManager *hmmManager, MAccumulatorPhysical *mAccumulatorPhysical);

		// destructor
		~GMMEditor();
		
		// initialize the editor
		void initialize();
		
		// merge Gaussian components in the mixture
		bool mixtureMerge(float fMinimumGaussianOccupation, float fMinimumGaussianWeight, float fCovarianceFlooringRatio);	
		
		// add Gaussian components to each mixture
		int mixtureIncrement(int iIncrement, int iSplittingCriterion, bool bAllowMultipleSplitsPerComponent, 
			float fMinimumGaussianOccupation, float fSplittingEpsilon);	
			
		// double the number of Gaussian components in the mixture (if possible)
		int mixtureDouble(float fMinimumGaussianOccupation, float fSplittingEpsilon);
	
};

#endif
