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


#include "GaussianMixture.h"
#include "Global.h"

namespace Bavieca {

GaussianMixture::GaussianMixture(int iDim, int iCovarianceType, int iComponents)
{
	m_iDim = iDim;
	m_iCovarianceType = iCovarianceType;
	for(int i=0 ; i < iComponents ; ++i) {
		m_vGaussian.push_back(new Gaussian(iDim,iCovarianceType));
	}
}

// destructor
GaussianMixture::~GaussianMixture()
{
	for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
		delete *it;
	}
}

// precompute evaluation constants
void GaussianMixture::precomputeConstant() {

	for(VGaussian::iterator it = m_vGaussian.begin() ; it != m_vGaussian.end() ; ++it) {
		(*it)->precomputeConstant();
	}
}

// remove a Gaussian component
void GaussianMixture::removeGaussianComponent(int g) {

	assert((g >= 0) && (g < (int)m_vGaussian.size()));

	delete m_vGaussian[g];
	// swap
	m_vGaussian[g] = m_vGaussian.back();
	m_vGaussian.pop_back();
}

// add a Gaussian component
void GaussianMixture::addGaussianComponent(Gaussian *gaussian) {

	assert(gaussian != NULL);	
	m_vGaussian.push_back(gaussian);
}

};	// end-of-namespace

