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


#ifndef HMMSTATE_H
#define HMMSTATE_H

#include <math.h>
#include <float.h>
#include <string.h>

using namespace std;

#include <list>

#include "Accumulator.h"
#include "GaussianMixture.h"
#include "Global.h"

#include "Vector.h"
#include "Matrix.h"
#include "SMatrix.h"

namespace Bavieca {

class Accumulator;
class FileInput;
class FileOutput;
class IOBase;
class PhoneSet;

// state's within-word position
#define WITHIN_WORD_POSITION_START			0
#define WITHIN_WORD_POSITION_INTERNAL		1
#define WITHIN_WORD_POSITION_END				2
#define WITHIN_WORD_POSITION_MONOPHONE		3

#define STR_WITHIN_WORD_POSITION_START			"start"
#define STR_WITHIN_WORD_POSITION_INTERNAL		"internal"
#define STR_WITHIN_WORD_POSITION_END			"end"
#define STR_WITHIN_WORD_POSITION_MONOPHONE	"monophone"

// Gaussian mixture size
#define HMM_MIXTURE_SIZE_SINGLE				0
#define HMM_MIXTURE_SIZE_MULTIPLE			1

class HMMState;

typedef vector<HMMState*> VHMMState;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class HMMState {

	private:
	
		int m_iDim;										// feature dimensionality
		int m_iCovarianceModeling;					// covariance modelling	
		int m_iCovarianceElements;					// relevant elements in the covariance matrix
		PhoneSet *m_phoneSet;
		GaussianMixture *m_gaussianMixture;
			
		// state identity
		unsigned char m_iPhone;						// basephone
		unsigned char m_iState;						// HMM-state 
		unsigned char m_iPosition;					// HMM-state within-word position
		int m_iId;										// HMM-state unique identifier
		
	public:
	
		// constructor
		HMMState(int iFeatureDimensionality, int iCovarianceModeling, PhoneSet *phoneSet, int iPhone, int iState, int iPosition, int iGaussians, int iId);
		
		// constructor (to be used when loading the HMM-state from a file)
		HMMState(int iFeatureDimensionality, int iCovarianceModeling, PhoneSet *phoneSet, int iId);
		
		// destructor
		~HMMState();
		
		// compute the emission probability of a state given the feature vector
		inline float computeEmissionProbability(float *fFeatures, int iTime)  {
		
			return m_gaussianMixture->evaluate(fFeatures,iTime);
		}
		
		// compute the emission probability of a state given the feature vector
		inline float computeEmissionProbabilityGaussian(int iGaussian, float *fFeatures, int iTime)  {
		
			return (*m_gaussianMixture)(iGaussian)->evaluate(fFeatures,iTime);
		}	
		
		// store the HMM into a file
		void store(FileOutput &file);
		
		// load the HMM from a file 
		void load(FileInput &file, unsigned char iEstimationMethod);
		
		// convert a within word position to string format
		inline const char *getStrPosition(int iPosition) {
		
			switch(iPosition) {
				case WITHIN_WORD_POSITION_START: {
					return STR_WITHIN_WORD_POSITION_START;
					break;
				}
				case WITHIN_WORD_POSITION_INTERNAL: {
					return STR_WITHIN_WORD_POSITION_INTERNAL;
					break;
				}
				case WITHIN_WORD_POSITION_END: {
					return STR_WITHIN_WORD_POSITION_END;
					break;
				}
				case WITHIN_WORD_POSITION_MONOPHONE: {
					return STR_WITHIN_WORD_POSITION_MONOPHONE;
					break;
				}
				default: {
					assert(0);
					exit(-1);	
				}
			}
			
			return NULL;
		}		
		
		// return the phone 
		inline unsigned char getPhone() {
			
			return m_iPhone;
		}
		
		// return the state
		inline unsigned char getState() {
		
			return m_iState;
		}
		
		//return the position
		inline unsigned char getPosition() {
		
			return m_iPosition;
		}		
		
		// return the HMM-state unique identifier
		inline int getId() {
		
			return m_iId;
		}
		
		// return the Gaussian mixture
		inline GaussianMixture &getMixture() {
			
			return *m_gaussianMixture;
		}
		
		// return the feature dimensionality
		inline int getFeatureDimensionality() {
		
			return m_iDim;
		}
		
		// return the covartiance modelling type
		inline int getCovarianceModelling() {
		
			return m_iCovarianceModeling;
		}
		
		// set the covartiance modelling type
		inline void setCovarianceModelling(int iCovarianceModeling) {
		
			m_iCovarianceModeling = iCovarianceModeling;
			if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				m_iCovarianceElements = m_iDim;
			} else {
				assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
				m_iCovarianceElements = (m_iDim*(m_iDim+1))/2;
			}
		}
		
		// set the feature dimensionality
		inline void setFeatureDimensionality(int iFeatureDimensionality) {
		
			m_iDim = iFeatureDimensionality;
		}
		
		// return whether an HMM-state number is valid
		inline static bool isStateValid(unsigned char iState) {
			
			if (iState < NUMBER_HMM_STATES) {
				return true;
			}
			return false;
		}
		
		// return whether a within-word posistion is valid
		inline static bool isPositionValid(unsigned char iPosition) {
			
			if ((iPosition != WITHIN_WORD_POSITION_START) && 
				(iPosition != WITHIN_WORD_POSITION_INTERNAL) &&
				(iPosition != WITHIN_WORD_POSITION_END) &&
				(iPosition != WITHIN_WORD_POSITION_MONOPHONE)) {
				return false;
			}
			return true;
		}

};

};	// end-of-namespace

#endif
