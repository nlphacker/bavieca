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


#include "ContextDecisionTree.h"
#include "HMMManager.h"
#include "PhoneSet.h"
#include "PhoneticRulesManager.h"
#include "Transform.h"

#include <iomanip>

namespace Bavieca {

// constructor
HMMManager::HMMManager(PhoneSet *phoneSet, unsigned char iPurpose)
{
	assert((iPurpose == HMM_PURPOSE_ESTIMATION) || (iPurpose == HMM_PURPOSE_EVALUATION));

	m_phoneSet = phoneSet;
	m_iPurpose = iPurpose;
	m_bInitialized = false;
	
	// HMMs (estimation)
	m_hmmStates = NULL;
	// HMMs (evaluation)
	m_hmmStatesDecoding = NULL;	
	
	// get the number of basephones
	m_iBasePhones = m_phoneSet->size();	
	
	// phonetic rules
	m_phoneticRulesManager = NULL;
	
	// context decision trees
	m_contextDecisionTrees = NULL;
	m_iContextDecisionTrees = -1;	
	
	// covariance
	m_iCovarianceModeling = -1;
	m_iCovarianceElements = -1;
	m_fCovarianceFloor = NULL;
	
	// feature dimensionality
	m_iDim = -1;
}

// destructor
HMMManager::~HMMManager()
{
	destroy();
}

// initializes the HMM's for the estimation
void HMMManager::initializeEstimation(unsigned char iAccumulatorType, unsigned char iContextModelingOrderWW, 
	unsigned char iContextModelingOrderCW) {

	m_iAccumulatorType = iAccumulatorType;

	// physical accumulators
	if (iAccumulatorType == ACCUMULATOR_TYPE_PHYSICAL) {	
	} 
	// logical accumulators
	else {
		assert(iAccumulatorType == ACCUMULATOR_TYPE_LOGICAL);
		m_iContextModelingOrderAccumulators = iContextModelingOrderWW;	
		m_iContextModelingOrderAccumulatorsCW = iContextModelingOrderCW;	
	}
}

// empty accumulators
void HMMManager::resetAccumulators() {

	// logical accumulators
	if (m_iAccumulatorType == ACCUMULATOR_TYPE_LOGICAL) {
		for(MAccumulatorLogical::iterator it = m_mAccumulatorLogical.begin() ; it != m_mAccumulatorLogical.end() ; ++it) {
			it->second->reset();
		}
	} 
	// physical accumulators
	else {	
		for(MAccumulatorPhysical::iterator it = m_mAccumulatorPhysical.begin() ; it != m_mAccumulatorPhysical.end() ; ++it) {
			it->second->reset();
		}
	}
}	


// destroy accumulators
void HMMManager::destroyAccumulators() {

	VAccumulator vAccumulator;

	// accumulators cannot be directly removed from the map without corrupting it
	
	// (1) logical accumulators	
	for(MAccumulatorLogical::iterator it = m_mAccumulatorLogical.begin() ; it != m_mAccumulatorLogical.end() ; ++it) {
		vAccumulator.push_back(it->second);
	}
	m_mAccumulatorLogical.clear();
	
	// (2) physical accumulators
	for(MAccumulatorPhysical::iterator it = m_mAccumulatorPhysical.begin() ; it != m_mAccumulatorPhysical.end() ; ++it) {
		vAccumulator.push_back(it->second);
	}
	m_mAccumulatorPhysical.clear();
		
	for(VAccumulator::iterator it = vAccumulator.begin() ; it != vAccumulator.end() ; ++it) {
		delete *it;	
	}
	vAccumulator.clear();
}

// create models' prototype
bool HMMManager::createSingleGaussianMonophoneModelsPrototype(int iDim, int iCovarianceModelling, int iComponents) {

	m_iDim = iDim;
	m_iCovarianceModeling = iCovarianceModelling;
	if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_iCovarianceElements = m_iDim;
	} 
	else {
		assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
		m_iCovarianceElements = (m_iDim*(m_iDim+1))/2;
	}	
	
	m_iHMMStates = m_iBasePhones * NUMBER_HMM_STATES;	
	m_iContextModelingOrderHMM = HMM_CONTEXT_MODELING_MONOPHONES;
	m_iContextModelingOrderHMMCW = HMM_CONTEXT_MODELING_MONOPHONES;
	m_iAccumulatorType = ACCUMULATOR_TYPE_LOGICAL;
	m_bSingleGaussian = true;
	assert(m_iHMMStates > 0);
	
	// allocate memory
	m_hmmStates = new HMMState*[m_iHMMStates];

	// create the HMM-states
	int iId = 0;
	for(int iPhone = 0 ; iPhone < m_iBasePhones ; ++iPhone) {
		for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			// it will be a "middle position" phone
			m_hmmStates[iPhone*NUMBER_HMM_STATES+iState] = new HMMState(m_iDim,m_iCovarianceModeling,m_phoneSet,iPhone,iState,WITHIN_WORD_POSITION_INTERNAL,iComponents,iId++);
		}
	}

	return true;
}

// convert the set of contex-independent HMM-states to context-dependent clustered units
void HMMManager::toContextDependentUnits(HMMState **hmmStates, int iHMMStates, unsigned char iContextModelingOrderWW,
	unsigned char iContextModelingOrderCW, PhoneticRulesManager *phoneticRulesManager, 
	ContextDecisionTree **contextDecisionTrees, int iContextDecisionTrees) {

	assert(m_iContextModelingOrderHMM == HMM_CONTEXT_MODELING_MONOPHONES);
	//assert(m_iContextModelingOrderHMMCW == HMM_CONTEXT_MODELING_MONOPHONES);
	
	m_iContextModelingOrderHMM = iContextModelingOrderWW;
	m_iContextModelingOrderHMMCW = iContextModelingOrderCW;
	
	// global accumulators will no longer be used
	m_iAccumulatorType = ACCUMULATOR_TYPE_PHYSICAL;
	// destroy global accumulators
	destroyAccumulators();	
	
	// delete the monophones
	if (m_hmmStates != NULL) {
		for(int i=0 ; i < m_iHMMStates; ++i) {
			delete m_hmmStates[i];
		}
		delete [] m_hmmStates;
		m_hmmStates = NULL;
	}

	m_hmmStates = hmmStates;
	m_iHMMStates = iHMMStates;
	
	// phonetic rules
	m_phoneticRulesManager = phoneticRulesManager;

	// context decision trees
	m_contextDecisionTrees = contextDecisionTrees;
	m_iContextDecisionTrees = iContextDecisionTrees;
}

// destroy
void HMMManager::destroy() {

	// delete the memory used to store the HMM-states
	// (estimation)
	if (m_iPurpose == HMM_PURPOSE_ESTIMATION) {
		if (m_hmmStates != NULL) {
			for(int i=0 ; i < m_iHMMStates; ++i) {
				delete m_hmmStates[i];
			}
			delete [] m_hmmStates;
			m_hmmStates = NULL;
		}	
		
		if (m_fCovarianceFloor != NULL) {
			delete [] m_fCovarianceFloor;
		}		
		destroyAccumulators();
	} 
	// (evaluation)
	else if (m_iPurpose == HMM_PURPOSE_EVALUATION) {	
		if (m_hmmStatesDecoding != NULL) {
			delete [] m_hmmStatesDecoding;
			m_hmmStatesDecoding = NULL;
		}	
	} else {
		assert(0);
	}
	
	if (m_phoneticRulesManager != NULL) {
		delete m_phoneticRulesManager;
	}
	if (m_contextDecisionTrees != NULL) {
		for(int i=0 ; i < m_iContextDecisionTrees ; ++i) {
			delete m_contextDecisionTrees[i];
		}
		delete [] m_contextDecisionTrees;	
	}
}

// load models from file
void HMMManager::load(const char *strFile) {

	try {

		FileInput file(strFile,true);
		file.open();
		
		// load the version of the system that created the models
		char strVersion[SYSTEM_VERSION_FIELD_WIDTH+1];
		IOBase::readBytes(file.getStream(),strVersion,SYSTEM_VERSION_FIELD_WIDTH);
		m_iSystemVersion = atoi(strVersion);
		
		// load the dimensionality and the covariance modeling type
		IOBase::read(file.getStream(),&m_iDim);
		IOBase::read(file.getStream(),&m_iCovarianceModeling);
		m_iCovarianceElements = getCovarianceElements(m_iDim,m_iCovarianceModeling);
		
		// load the number of phonemes
		int iPhones = -1;
		IOBase::read(file.getStream(),&iPhones);
		assert(iPhones == (int)m_phoneSet->size());
		
		// load the context modeling type
		IOBase::read(file.getStream(),&m_iContextModelingOrderHMM);
		IOBase::read(file.getStream(),&m_iContextModelingOrderHMMCW);
		
		// load whether the models are single or multiple gaussian
		m_bSingleGaussian = false;
		unsigned char iGaussianMixtureSize = UCHAR_MAX;
		IOBase::read(file.getStream(),&iGaussianMixtureSize);
		if (iGaussianMixtureSize == HMM_MIXTURE_SIZE_SINGLE) {
			m_bSingleGaussian = true;
		} else {
			assert(iGaussianMixtureSize == HMM_MIXTURE_SIZE_MULTIPLE);
			m_bSingleGaussian = false;
		}
		
		// load the number of different HMM-states (monophones or triphones)
		IOBase::read(file.getStream(),&m_iHMMStates);
		
		// HMM-estimation
		if (m_iPurpose == HMM_PURPOSE_ESTIMATION) { 
		
			// allocate memory for the states
			m_hmmStates = new HMMState*[m_iHMMStates];
			
			// load the states
			for(int i=0 ; i < m_iHMMStates; ++i) {
				m_hmmStates[i] = new HMMState(m_iDim,m_iCovarianceModeling,m_phoneSet,i);
				m_hmmStates[i]->load(file,true);
			}
		} 
		// HMM-evaluation
		else {
			
			assert(m_iPurpose == HMM_PURPOSE_EVALUATION);
			
			// allocate memory for the states
			m_hmmStatesDecoding = new HMMStateDecoding[m_iHMMStates];
			
			// set initial parameters and load the states
			int iGaussianId = 0;
			for(int i=0 ; i < m_iHMMStates; ++i) {
				m_hmmStatesDecoding[i].setInitialParameters(m_iDim,m_phoneSet,i);
				m_hmmStatesDecoding[i].load(file);
				m_hmmStatesDecoding[i].setGaussianIds(&iGaussianId);
			}
		}
		
		// context dependent models (it is necessary to load the decision trees)
		if (m_iContextModelingOrderHMM > HMM_CONTEXT_MODELING_MONOPHONES) {
		
			// local accumulators have to be used
			m_iAccumulatorType = ACCUMULATOR_TYPE_PHYSICAL;
			
			// load the phonetic rules 
			m_phoneticRulesManager = PhoneticRulesManager::load(file,m_phoneSet);
			
			// load the # of context decision tree(s)
			IOBase::read(file.getStream(),&m_iContextDecisionTrees);
			// load the context decision tree(s)
			if (m_iContextDecisionTrees == 1) {
				m_contextDecisionTrees = new ContextDecisionTree*[1];
				m_contextDecisionTrees[0] = ContextDecisionTree::load(file,m_iDim,m_iCovarianceModeling,
					m_phoneSet,m_phoneticRulesManager,m_iContextModelingOrderHMM);
				assert(m_contextDecisionTrees[0]);
			} else {
				assert(m_iContextDecisionTrees == (int)m_phoneSet->size()*NUMBER_HMM_STATES);
				m_contextDecisionTrees = new ContextDecisionTree*[m_iContextDecisionTrees];
				for(int i=0 ; i < m_iContextDecisionTrees ; ++i) {
					m_contextDecisionTrees[i] = ContextDecisionTree::load(file,m_iDim,m_iCovarianceModeling,
						m_phoneSet,m_phoneticRulesManager,m_iContextModelingOrderHMM);
					assert(m_contextDecisionTrees[i]);
					assert(m_contextDecisionTrees[i]->checkConsistency());
				}
			}
		}
		
		// mark the HMMs as initialized
		m_bInitialized = true;
		
	} catch(ExceptionBase) {
		BVC_ERROR << "unable to load the acoustic models from file: " << strFile;
	}
}

// store models into a file
void HMMManager::store(const char *strFile) {

	try {

		FileOutput file(strFile,true);
		file.open();
		
		// store the system version
		char strVersion[SYSTEM_VERSION_FIELD_WIDTH+1];
		memset(strVersion,0,SYSTEM_VERSION_FIELD_WIDTH+1);
		assert(strlen(SYSTEM_VERSION) <= SYSTEM_VERSION_FIELD_WIDTH);
		strcpy(strVersion,SYSTEM_VERSION);
		IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(strVersion),SYSTEM_VERSION_FIELD_WIDTH);
		
		IOBase::write(file.getStream(),m_iDim);
		IOBase::write(file.getStream(),m_iCovarianceModeling);
		int iPhones = m_phoneSet->size();
		IOBase::write(file.getStream(),iPhones);
		
		IOBase::write(file.getStream(),m_iContextModelingOrderHMM);
		IOBase::write(file.getStream(),m_iContextModelingOrderHMMCW);
		
		// store whether models are single or multiple gaussian
		unsigned char iGaussianMixtureSize = HMM_MIXTURE_SIZE_SINGLE;
		if (m_bSingleGaussian == false) {
			iGaussianMixtureSize = HMM_MIXTURE_SIZE_MULTIPLE;
		}
		IOBase::write(file.getStream(),iGaussianMixtureSize);	
		IOBase::write(file.getStream(),m_iHMMStates);
		
		// store the states
		if (m_iPurpose == HMM_PURPOSE_ESTIMATION) {
			for(int i=0 ; i < m_iHMMStates; ++i) {
				HMMState *hmmState = m_hmmStates[i];
				hmmState->store(file);
			}
		} else {
			for(int i=0 ; i < m_iHMMStates; ++i) {	
				HMMStateDecoding *hmmStateDecoding = &m_hmmStatesDecoding[i];
				hmmStateDecoding->store(file);
			}
		}		
		
		// context dependent models
		if (m_iContextModelingOrderHMM > HMM_CONTEXT_MODELING_MONOPHONES) {		
			
			// store the phonetic rules 
			assert(m_phoneticRulesManager != NULL);
			m_phoneticRulesManager->store(file);
			
			// store the # of context decision tree(s)
			IOBase::write(file.getStream(),m_iContextDecisionTrees);
			
			// store the context decision tree(s)
			for(int i=0 ; i < m_iContextDecisionTrees ; ++i) {
				m_contextDecisionTrees[i]->store(file);
			}	
		}		
		
		file.close();
	} catch(ExceptionBase) {
		BVC_ERROR << "unable to store the acoustic models to file: " << strFile;
	}
}

// return a pointer to the array of HMM-states
HMMState **HMMManager::getHMMStates(int *iStates) {

	*iStates = m_iHMMStates;
	return m_hmmStates;
}

// return a pointer to a HMM-state given its identity
HMMState *HMMManager::getHMMState(unsigned char *iPhoneLeft, unsigned char iPhone, unsigned char *iPhoneRight, 
	unsigned char iPosition, unsigned char iState) {
	
	int iIndex = getHMMStateIndex(iPhoneLeft,iPhone,iPhoneRight,iPosition,iState);
	return m_hmmStates[iIndex];
}

// return a pointer to a HMM-state given its id (index)
HMMState *HMMManager::getHMMState(int iId) {

	assert(m_iPurpose == HMM_PURPOSE_ESTIMATION);
	
	if ((iId < 0) || (iId >= m_iHMMStates)) {
		return NULL;
	}
	assert(m_hmmStates[iId]->getId() == iId);
	
	return m_hmmStates[iId];
}

// return a pointer to a HMM-state given its identity
HMMStateDecoding *HMMManager::getHMMStateDecoding(unsigned char *iPhoneLeft, unsigned char iPhone, 
	unsigned char *iPhoneRight, unsigned char iPosition, unsigned char iState) {
	
	int iIndex = getHMMStateIndex(iPhoneLeft,iPhone,iPhoneRight,iPosition,iState);
	return &m_hmmStatesDecoding[iIndex];	
}

// return a pointer to a HMM-state given its identity
unsigned int HMMManager::getHMMStateIndex(unsigned char *iPhoneLeft, unsigned char iPhone, 
	unsigned char *iPhoneRight, unsigned char iPosition, unsigned char iState) {

	// monophones
	if (m_iContextModelingOrderHMM == HMM_CONTEXT_MODELING_MONOPHONES) {	
	
		return iPhone*NUMBER_HMM_STATES+iState;	
	} 
	// context dependent nphones
	else {
	
		// single tree (global clustering)
		if (m_iContextDecisionTrees == 1) {
			int iIndex = (*m_contextDecisionTrees)->getHMMIndex(iPhoneLeft,iPhone,iPhoneRight,iPosition,iState);
			assert((iIndex >= 0) && (iIndex < m_iHMMStates));
			return iIndex;
		} 
		// local clustering
		else {
			int iIndex = m_contextDecisionTrees[iPhone*NUMBER_HMM_STATES+iState]->getHMMIndex(iPhoneLeft,
				iPhone,iPhoneRight,iPosition,iState);
			assert((iIndex >= 0) && (iIndex < m_iHMMStates));
			return iIndex;
		}
	}
}


// return the number of free parameters
int HMMManager::getNumberFreeParameters() {

	int iFreeParameters = 0;

	// accumulate free parameters across HMM-states
	for(int i=0 ; i<m_iHMMStates ; ++i) {
		iFreeParameters += m_hmmStates[i]->getMixture().getNumberFreeParameters();
	}

	return iFreeParameters;
}

// return the number of free parameters
int HMMManager::getNumberGaussianComponents() {

	int iGaussianComponents = 0;

	// compute the new distributions
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		if (m_hmmStates != NULL) {
			iGaussianComponents += m_hmmStates[i]->getMixture().getNumberComponents();
		} else {
			iGaussianComponents += m_hmmStatesDecoding[i].getGaussianComponents();
		}
	}

	return iGaussianComponents;
}


// initialize the acoustic model parameters to the given mean and covariance
void HMMManager::initializeModels(Vector<float> &vMean, SMatrix<float> &mCovariance) {

	for(unsigned int iBasePhone = 0 ; iBasePhone < m_phoneSet->size() ; ++iBasePhone) {
		for(int iState = 0 ; iState < NUMBER_HMM_STATES ; ++iState) {
			HMMState *hmmState = m_hmmStates[iBasePhone*NUMBER_HMM_STATES+iState];
			if (hmmState == NULL) {
				continue;
			}
			assert(hmmState != NULL);
			hmmState->getMixture().setMeanAllComponents(vMean);
			if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				Vector<float> vCovariance(m_iDim);
				mCovariance.getDiagonal(vCovariance);
				hmmState->getMixture().setCovarianceAllComponents(vCovariance);
			} else {
				hmmState->getMixture().setCovarianceAllComponents(mCovariance);
			}
			assert(hmmState->getMixture().getNumberComponents() == 1);
			hmmState->getMixture()(0)->weight() = 1.0;
		}	
	}
}

// initialize the models from other models
void HMMManager::initializeModels(HMMManager *hmmManagerSeed, int iDim, int iCovarianceModelling) {

	assert(m_hmmStates == NULL);
	assert(m_iDim == -1);
	
	// the feature dimensionality might be different than that of the original models
	m_iDim = iDim;
	
	m_iCovarianceModeling = iCovarianceModelling;
	m_iCovarianceElements = getCovarianceElements(iDim,iCovarianceModelling);
	m_iContextModelingOrderHMM = hmmManagerSeed->getContextModelingOrderHMM();
	m_iContextModelingOrderHMMCW = hmmManagerSeed->getContextModelingOrderHMMCW();
	m_iContextModelingOrderAccumulators = hmmManagerSeed->getContextSizeAccumulators();
	m_iContextModelingOrderAccumulatorsCW = hmmManagerSeed->getContextSizeAccumulatorsCW();
	m_bSingleGaussian = hmmManagerSeed->isSingleGaussian();
	
	// accumulators should be analogous to those of the original models
	m_iAccumulatorType = hmmManagerSeed->getAccumulatorType();

	m_iHMMStates = -1;
	HMMState **hmmStates = hmmManagerSeed->getHMMStates(&m_iHMMStates);
	m_hmmStates = new HMMState*[m_iHMMStates];	
	for(int i=0 ; i<m_iHMMStates ; ++i) {
		m_hmmStates[i] = new HMMState(m_iDim,m_iCovarianceModeling,m_phoneSet,hmmStates[i]->getPhone(),
			hmmStates[i]->getState(),hmmStates[i]->getPosition(),hmmStates[i]->getMixture().getNumberComponents(),
			hmmStates[i]->getId());
		if (m_hmmStates[i]->getMixture().getNumberComponents() == 1) {
			m_hmmStates[i]->getMixture()(0)->weight() = 1.0;
		}
	}
	
	// context decision trees
	if (m_iContextModelingOrderHMM > HMM_CONTEXT_MODELING_MONOPHONES) {
		m_iContextDecisionTrees = hmmManagerSeed->m_iContextDecisionTrees;
		m_contextDecisionTrees = hmmManagerSeed->m_contextDecisionTrees;
		m_phoneticRulesManager = hmmManagerSeed->m_phoneticRulesManager;
	}
	
	m_bInitialized = true;
}

// reset the emission probability computation
void HMMManager::resetHMMEmissionProbabilityComputation() {
	
	// estimation
	if (m_iPurpose == HMM_PURPOSE_ESTIMATION) {
		for(int i=0 ; i<m_iHMMStates ; ++i) {
			m_hmmStates[i]->getMixture().resetTimeStamp();
		}
	}
	// decoding 
	else {
		for(int i=0 ; i<m_iHMMStates ; ++i) {
			m_hmmStatesDecoding[i].resetTimeStamp();
		}	
	}
}

// reset the emission probability computation of the given HMM-states
void HMMManager::resetHMMEmissionProbabilityComputation(VHMMState &vHMMState) {

	for(VHMMState::iterator it = vHMMState.begin() ; it != vHMMState.end() ; ++it) {
		(*it)->getMixture().resetTimeStamp();
	}
}

// reset the emission probability computation of the given HMM-states
void HMMManager::resetHMMEmissionProbabilityComputation(VHMMStateDecoding &vHMMStateDecoding) {

	for(VHMMStateDecoding::iterator it = vHMMStateDecoding.begin() ; it != vHMMStateDecoding.end() ; ++it) {
		(*it)->resetTimeStamp();
	}
}

// print models information
void HMMManager::printModelsInfo(bool bVerbose) {

	float fGaussiansPerState = 0.0;
	int iGaussians = getNumberGaussianComponents();
	int iFreeParameters = getNumberFreeParameters();
	
	cout << "- acoustic models information ----------------------------------------" << endl;
	cout << " # phones:                         " << setw(3) << m_iBasePhones << endl;
	cout << " feature dimensionality:           " << setw(3) << m_iDim << endl;
	cout << " HMM-topology:                     " << NUMBER_HMM_STATES << " states left-to-right" << endl;	
	cout << " # physical HMM-states:            " << m_iHMMStates << endl;	
	if (m_iContextModelingOrderHMM != HMM_CONTEXT_MODELING_MONOPHONES) {	
		cout << "# logical HMM-states: " << m_mAccumulatorLogical.size() << endl;
	}	
	cout << " covariance modelling:   " << setw(16) << Gaussian::getCovarianceModellingStringFormat(m_iCovarianceModeling) << endl;
	cout << " context modeling order: " << setw(16) << Accumulator::getContextModelingOrder(m_iContextModelingOrderHMM) << endl;	
	fGaussiansPerState = ((float)iGaussians)/((float)m_iHMMStates);
	cout << "# Gaussian components:         " << iGaussians << " (" << fGaussiansPerState << " Gaussian components per state on average)" << endl;
	cout << "# free parameters:             " << iFreeParameters << endl;
	cout << "----------------------------------------------------------------------" << endl;
	
	if (bVerbose) {
		
		
	}
}

// initialize HMM-states for decoding
void HMMManager::initializeDecoding() {
	
	m_iAccumulatorType = ACCUMULATOR_TYPE_PHYSICAL;
	for(int i=0 ; i<m_iHMMStates ; ++i) {
		m_hmmStatesDecoding[i].initialize();
	}
}

// precompute constants to speed-up emission probability computation
void HMMManager::precomputeConstants() {

	for(int i=0 ; i<m_iHMMStates ; ++i) {
		m_hmmStates[i]->getMixture().precomputeConstant();
	}
}

// dump accumulators (either logical or physical) to disk
void HMMManager::dumpAccumulators(const char *strFile) {

	try {
		// logical accumulators
		if (m_iAccumulatorType == ACCUMULATOR_TYPE_LOGICAL) {	
			Accumulator::storeAccumulators(strFile,m_iDim,m_iCovarianceModeling,
				m_iContextModelingOrderAccumulators,m_iContextModelingOrderAccumulatorsCW,m_mAccumulatorLogical);
		}
		// physical accumulators
		else {		
			Accumulator::storeAccumulators(strFile,m_iDim,m_iCovarianceModeling,
				m_iHMMStates,getNumberGaussianComponents(),m_mAccumulatorPhysical);
		}
	} catch (ExceptionBase) {
		BVC_ERROR << "unable to load the accumulators";
	}	
}

};	// end-of-namespace

