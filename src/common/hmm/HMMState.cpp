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


#include "HMMState.h"
#include "IOBase.h"
#include "PhoneSet.h"
#include "FileInput.h"
#include "FileOutput.h"

namespace Bavieca {

// constructor
HMMState::HMMState(int iDim, int iCovarianceModeling, PhoneSet *phoneSet, int iPhone, 
	int iState, int iPosition, int iGaussianComponents, int iId)
{
	m_iDim = iDim;
	m_iCovarianceModeling = iCovarianceModeling;	
	m_phoneSet = phoneSet;
	
	// compute the number of elements in the covariance of each Gaussian component
	if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_iCovarianceElements = m_iDim;
	} 
	else {
		assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
		m_iCovarianceElements = m_iDim*(m_iDim+1)/2;
	}
	
	assert(iState < NUMBER_HMM_STATES);
	
	// state identity
	m_iPhone = iPhone;
	m_iState = iState;
	m_iPosition = iPosition;
	m_iId = iId;
	
 	// mixture
	m_gaussianMixture = new GaussianMixture(m_iDim,m_iCovarianceModeling,iGaussianComponents);
}

// constructor (to be used when loading the HMM-state from a file)
HMMState::HMMState(int iDim, int iCovarianceModeling, PhoneSet *phoneSet, int iId)
{
	m_iDim = iDim;
	m_iCovarianceModeling = iCovarianceModeling;
	m_phoneSet = phoneSet;
	m_iId = iId;
	
	// compute the number of elements in the covariance of each Gaussian component
	if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_iCovarianceElements = m_iDim;
	} 
	else {
		assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
		m_iCovarianceElements = (m_iDim*(m_iDim+1))/2;
	}
	
 	// mixture
	m_gaussianMixture = new GaussianMixture(m_iDim,m_iCovarianceModeling,0);
}

// destructor
HMMState::~HMMState()
{
	delete m_gaussianMixture;
}

// store the HMM into a file
void HMMState::store(FileOutput &file) {
	
	// phonetic symbol
	char strPhone[MAX_PHONETIC_SYMBOL_LENGTH+1];
	memset(strPhone,0,MAX_PHONETIC_SYMBOL_LENGTH+1);
	strcpy(strPhone,m_phoneSet->getStrPhone(m_iPhone));	
	IOBase::writeBytes(file.getStream(),reinterpret_cast<char*>(strPhone),MAX_PHONETIC_SYMBOL_LENGTH+1);
	
	// state
	IOBase::write(file.getStream(),m_iState);
	
	// within word position (DEPRECATED)
	IOBase::write(file.getStream(),m_iPosition);

	// Gaussian components
	int iGaussianComponents = m_gaussianMixture->getNumberComponents();
	IOBase::write(file.getStream(),iGaussianComponents);
	for(int iGaussian = 0 ; iGaussian < iGaussianComponents ; ++iGaussian) {
		Gaussian *gaussian = (*m_gaussianMixture)(iGaussian); 
		IOBase::write(file.getStream(),gaussian->weight());
		gaussian->mean().writeData(file.getStream());
		if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
			gaussian->covarianceDiag().writeData(file.getStream());
		} else {
			gaussian->covarianceFull().writeData(file.getStream());
		}
	}
}

// load the HMM from a file
void HMMState::load(FileInput &file, unsigned char iEstimationMethod) {

	// phonetic symbol
	char strPhone[MAX_PHONETIC_SYMBOL_LENGTH+1];	
	IOBase::readBytes(file.getStream(),reinterpret_cast<char*>(strPhone),MAX_PHONETIC_SYMBOL_LENGTH+1);
	m_iPhone = m_phoneSet->getPhoneIndex(strPhone);
	assert(m_iPhone != UCHAR_MAX);
	
	// state
	IOBase::read(file.getStream(),&m_iState);
	assert(m_iState < NUMBER_HMM_STATES);
	
	// within word position (DEPRECATED)
	IOBase::read(file.getStream(),&m_iPosition);
	
	// Gaussian components
	int iGaussianComponents = -1;
	IOBase::read(file.getStream(),&iGaussianComponents);
	for(int iGaussian = 0 ; iGaussian < iGaussianComponents ; ++iGaussian) {
	
		Gaussian *gaussian = new Gaussian(m_iDim,m_iCovarianceModeling);	
		
		IOBase::read(file.getStream(),&gaussian->weight());
		gaussian->mean().readData(file.getStream());		
		if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
			gaussian->covarianceDiag().readData(file.getStream());
		} else {
			gaussian->covarianceFull().readData(file.getStream());
		}	
		m_gaussianMixture->addGaussianComponent(gaussian);
	}
}

};	// end-of-namespace

