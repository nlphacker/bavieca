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


#ifndef HMMMANAGER_H
#define HMMMANAGER_H

#include <limits.h>

using namespace std;

#include <list>

#include "HMMState.h"
#include "HMMStateDecoding.h"

namespace Bavieca {

class ContextDecisionTree;
class PhoneSet;
class PhoneticRulesManager;
class Transform;

// models purpose
#define HMM_PURPOSE_ESTIMATION			0
#define HMM_PURPOSE_EVALUATION			1

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class HMMManager {

	private:
	
		PhoneSet *m_phoneSet;						// phonetic symbol set
		unsigned char m_iPurpose;					// whether HMMs will be used to estimation or just evaluation
		
		// HMM-states (estimation)
		HMMState **m_hmmStates;						// monophones: [iBasePhone x iState x iPosition] 
															// triphones: array of physical clustered context dependent HMM-states
		
		// HMM-states (evaluation)
		HMMStateDecoding *m_hmmStatesDecoding;	// monophones: [iBasePhone x iState x iPosition] 
															// triphones: array of physical clustered context dependent HMM-states
		
		int m_iDim;										// feature dimensionality
		int m_iCovarianceModeling;					// covariance modeling type (diagonal/full)
		int m_iCovarianceElements;
		int m_iHMMStates;								// number of unique HMM-states
		bool m_bInitialized;							// whether the models are already initialized
		int m_iBasePhones;							// number of basephones
		unsigned char m_iContextModelingOrderHMM;					// HMMs context modeling order (within-word)
		unsigned char m_iContextModelingOrderHMMCW;				// HMMs context modeling order (cross-word)
		unsigned char m_iContextModelingOrderAccumulators;		// Accumulator's context modeling order (within-word)
		unsigned char m_iContextModelingOrderAccumulatorsCW;	// Accumulator context modeling order (cross-word)
		unsigned char m_iIdentityAux[MAX_IDENTITY_LENGTH];
		bool m_bSingleGaussian;						// whether all the HMM-states are single gaussian
		
		// accumulators
		unsigned char m_iAccumulatorType;		// whether to use logical accumulators (logical vs physical)
															// logical acc: one accumulator per context modeling unit
															// physical acc: one accumulator per Gaussian component
		
		// logical accumulators (needed for context clustering)
		MAccumulatorLogical m_mAccumulatorLogical;		// accumulators (one for each logical triphone)
		MAccumulatorPhysical m_mAccumulatorPhysical;		// accumulators (one for each Gaussian component)
		
		// phonetic rules
		PhoneticRulesManager *m_phoneticRulesManager;
		// context decision trees
		ContextDecisionTree **m_contextDecisionTrees;
		int m_iContextDecisionTrees;		
		
		// destroy
		void destroy();		
		
		// destroy logical accumulators
		void destroyAccumulators();
		
		// covariance floor
		float *m_fCovarianceFloor;
		
		// version of the hmmsystem that created the models
		int m_iSystemVersion;
		
	public:
	
		// constructor
		HMMManager(PhoneSet *phoneSet, unsigned char iPurpose);

		// destructor
		~HMMManager();
				
		// initializes the HMM's for the estimation
		void initializeEstimation(unsigned char iAccumulatorType, unsigned char iContextModelingOrderWW, 
			unsigned char iContextModelingOrderCW);
		
		// set the context modeling order of the accumulators
		void setContextModelingOrderAccumulators(unsigned char iContextModelingOrderAccumulators);	
		
		// get the covariance elements (dimensionality)
		inline static int getCovarianceElements(int iDim, int iCovarianceModelling) {
		
			if (iCovarianceModelling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				return iDim;
			} 
			else {
				assert(iCovarianceModelling == COVARIANCE_MODELLING_TYPE_FULL);
				return (iDim*(iDim+1))/2;
			}	
		}
		
		// set the feature dimensionality
		inline void setFeatureDimensionality(int iDim) {
		
			m_iDim = iDim;
		}
				
		// set the covariance modelling type
		// set the covartiance modelling type
		inline void setCovarianceModelling(int iCovarianceModeling) {
		
			m_iCovarianceModeling = iCovarianceModeling;
			m_iCovarianceElements = getCovarianceElements(m_iDim,iCovarianceModeling);
			
			for(int i=0 ; i < m_iHMMStates ; ++i) {
				m_hmmStates[i]->setCovarianceModelling(iCovarianceModeling);
			}
		}		
		
		// return the logical accumulators
		inline MAccumulatorLogical &getAccumulators() {
		
			return m_mAccumulatorLogical;
		}
		
		// return the phonetic symbol associated to an HMM-state
		inline unsigned char getPhoneticSymbolFromHMMStateDecoding(unsigned int iHMMStateDecoding) {
		
			assert(m_iPurpose == HMM_PURPOSE_EVALUATION);
		
			assert((int)iHMMStateDecoding <= m_iHMMStates);
			
			return m_hmmStatesDecoding[iHMMStateDecoding].getPhone();	
		}
		
		// get the accumulator for the given n-phone (logical accumulator)
		// note: accumulators are created on demand
		inline Accumulator *getAccumulator(unsigned char *iPhoneLeft, unsigned char iPhone, unsigned char *iPhoneRight,
			unsigned char iPosition, unsigned char iState) {
		
			Accumulator::buildIdentity(m_iIdentityAux,iPhoneLeft,iPhone,iPhoneRight,iPosition,iState,
				m_iContextModelingOrderAccumulators);
			
			MAccumulatorLogical::iterator it = m_mAccumulatorLogical.find(m_iIdentityAux);
			if (it == m_mAccumulatorLogical.end()) {
				
				Accumulator *accumulator = new Accumulator(m_iDim,m_iCovarianceModeling,m_iIdentityAux,m_iContextModelingOrderAccumulators);	m_mAccumulatorLogical.insert(MAccumulatorLogical::value_type(accumulator->getIdentity(),accumulator));
				return accumulator;
			} 
			else {	
				return it->second;
			}
		}	
		
		// get the accumulator for the given Gaussian component
		// note: accumulators are created on demand
		inline Accumulator *getAccumulator(int iHMMState, int iGaussian) {
		
			unsigned int iKey = Accumulator::getPhysicalAccumulatorKey(iHMMState,iGaussian);
			MAccumulatorPhysical::iterator it = m_mAccumulatorPhysical.find(iKey);
			if (it == m_mAccumulatorPhysical.end()) {
				Accumulator *accumulator = new Accumulator(m_iDim,m_iCovarianceModeling,iHMMState,iGaussian);
				m_mAccumulatorPhysical.insert(MAccumulatorPhysical::value_type(iKey,accumulator));
				return accumulator;
			} else {
				return it->second;
			}
		}	
		
		// get the accumulator for Gaussian components in the given HMM state
		// note: accumulators are created on demand
		inline void getAccumulators(int iHMMState, VAccumulator *vAccumulator) {
		
			HMMState *hmmState = getHMMState(iHMMState);
			for(unsigned int g=0 ; g<hmmState->getMixture().getNumberComponents() ; ++g) {
				vAccumulator->push_back(getAccumulator(iHMMState,g));
			}
		}	
		
		// return the purpose
		inline unsigned char getPurpose() {
		
 			return m_iPurpose;
		}
		
		// create models' prototype
		bool createSingleGaussianMonophoneModelsPrototype(int iDim, int iCovarianceModelling, int iComponents);
		
		// convert the set of contex-independent HMM-states to context-dependent clustered units
		void toContextDependentUnits(HMMState **hmmStates, int iHMMStates, unsigned char iContextModelingOrderWW, 
			unsigned char iContextModelingOrderCW, PhoneticRulesManager *phoneticRulesManager, 
			ContextDecisionTree **contextDecisionTrees, int iContextDecisionTrees);	
		
		// load models from file
		void load(const char *strFile);
		
		// store models into a file
		void store(const char *strFile);
		
		// return a pointer to the array of HMM-states
		HMMState **getHMMStates(int *iStates);	
		
		// return a pointer to a HMM-state given its identity
		HMMState *getHMMState(unsigned char *iPhoneLeft, unsigned char iPhone, unsigned char *iPhoneRight, 
			unsigned char iPosition, unsigned char iState);	
			
		// return a pointer to a HMM-state given its id (index)
		HMMState *getHMMState(int iId);
		
		// return the index of the HMM
		unsigned int getHMMStateIndex(unsigned char *iPhoneLeft, unsigned char iPhone, unsigned char *iPhoneRight, 
			unsigned char iPosition, unsigned char iState);
		
		// return a pointer to the array of HMM-states
		inline HMMStateDecoding *getHMMStatesDecoding(int *iHMMStates) {
		
			*iHMMStates = m_iHMMStates;
		
			return m_hmmStatesDecoding;
		}
		
		// return a pointer to a HMM-state given its identity
		HMMStateDecoding *getHMMStateDecoding(unsigned char *iPhoneLeft, unsigned char iPhone, unsigned char *iPhoneRight, unsigned char iPosition, unsigned char iState);	
		
		// return a pointer to a HMM-state given its index
		inline HMMStateDecoding *getHMMStateDecoding(int iIndex) {
			
			return &m_hmmStatesDecoding[iIndex];	
		}
		
		// return the number of free parameters
		int getNumberFreeParameters();	
		
		// return the number of Gaussian components in the system
		int getNumberGaussianComponents();
		
		// return the number of physical HMM-states
		inline int getNumberHMMStatesPhysical() {
		
			return m_iHMMStates;
		}
		
		// initialize the model parameters to the given mean and covariance
		void initializeModels(Vector<float> &vMean, SMatrix<float> &mCovariance);	
		
		// initialize the models from a list of logical accumulators
		bool initializeModels(MAccumulatorLogical &mAccumulatorLogical, int iDim, int iCovarianceModelling);
		
		// return whether the models are already initialized
		inline bool areInitialized() {
		
			return m_bInitialized;
		}

		// change the initialization state of the models
		inline void setInitialized(bool bInitialized) {
		
			m_bInitialized = bInitialized;
		}
		
		// reset the HMM accumulators used for the parameter estimation
		void resetAccumulators();
		
		// update HMM-parameters
		void updateHMMParameters(bool bUpdateCovariance, float fCovarianceFlooringRatio);	
		
		// update HMM-parameters
		void updateHMMParametersMAP(float fPriorKnowledgeWeight, bool bUpdateCovariance, int *iGaussianComponentsUpdated);	
		
		// reset the emission probability computation
		void resetHMMEmissionProbabilityComputation();	
				
		// reset the emission probability computation of the given HMM-states
		void resetHMMEmissionProbabilityComputation(VHMMState &vHMMState);
		
		// reset the emission probability computation of the given HMM-states
		void resetHMMEmissionProbabilityComputation(VHMMStateDecoding &vHMMStateDecoding);
		
		// print models information
		void printModelsInfo(bool bVerbose);
		
		// initialize HMM-states for decoding
		void initializeDecoding();
		
		inline void setSingleGaussian(bool b) {
		
			m_bSingleGaussian = b;
		}
		
		// return whether the HMMs are single gaussian
		inline bool isSingleGaussian() {
		
			return m_bSingleGaussian;
		}
		
		// return whether the accumulators are global
		inline bool areAccumulatorsLogical() {
		
			return (m_iAccumulatorType == ACCUMULATOR_TYPE_LOGICAL);
		}
		
		// precompute constants to speed-up emission probability computation
		void precomputeConstants();
		
		// return the context modeling size
		inline unsigned char getContextSizeAccumulators() {
		
			return (m_iContextModelingOrderAccumulators-1)/2;
		}
		
		// return the context modeling size (cross-word)
		inline unsigned char getContextSizeAccumulatorsCW() {
		
			return (m_iContextModelingOrderAccumulatorsCW-1)/2;
		}
		
		// return the context modeling size
		inline unsigned char getContextSizeHMM() {
		
			return (m_iContextModelingOrderHMM-1)/2;
		}
		
		// return the context modeling size
		inline unsigned char getContextSizeHMMCW() {
		
			return (m_iContextModelingOrderHMMCW-1)/2;
		}
		
		// return the context modeling order
		inline unsigned char getContextModelingOrderHMM() {
		
			return m_iContextModelingOrderHMM;
		}
		
		// return the context modeling order
		inline unsigned char getContextModelingOrderHMMCW() {
		
			return m_iContextModelingOrderHMMCW;
		}
		
		// return the feature dimensionality
		inline unsigned int getFeatureDim() {
		
			return m_iDim;
		}
		
		// return the dimensionality of a feature container (considers memory alignment)
		inline unsigned int getFeatureDimContainer() {
		#if defined __AVX__ || defined __SSE__
			// round up the dimensionality to the nearest multiple	
			unsigned int iReminder = m_iDim%(ALIGN_BOUNDARY/sizeof(float));
			unsigned int iDimContainer = (iReminder == 0) ? m_iDim : m_iDim + (ALIGN_BOUNDARY/sizeof(float))-iReminder;
			assert(iDimContainer % ALIGN_BOUNDARY == 0);
			return iDimContainer;
		#else 
			return m_iDim;
		#endif	
		}	
		
		// return the covariance modelling type
		inline int getCovarianceModelling() {
		
			return m_iCovarianceModeling;
		}
		
		// return the covariance elements
		inline int getCovarianceElements() {
		
			return m_iCovarianceElements;
		}
		
		// return the accumulator type
		inline unsigned char getAccumulatorType() {
		
			return m_iAccumulatorType;
		}
		
		// initialize the models from other models
		void initializeModels(HMMManager *hmmManagerSeed, int iDim, int iCovarianceModellingUpdate);
		
		// dump accumulators (either logical or physical) to disk
		void dumpAccumulators(const char *strFile);

};

};	// end-of-namespace

#endif
