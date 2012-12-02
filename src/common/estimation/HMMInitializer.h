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


#ifndef HMMINITIALIZER_H
#define HMMINITIALIZER_H

using namespace std;

#include <string>

#include "HMMManager.h"
#include "MLFFile.h"

#include "SMatrix.h"
#include "Vector.h"
#include "VectorStatic.h"

namespace Bavieca {

class HMMManager;
class PhoneSet;

#define STATE_MODELS_UNINITIALIZED			0
#define STATE_MODELS_INITIALIZED				1

// initialization methods
#define AM_INITIALIZATION_METHOD_FLATSTART	"flatStart"
#define AM_INITIALIZATION_METHOD_BOOTSTRAP	"bootstrap"
#define AM_INITIALIZATION_METHOD_ALIGNMENT	"alignment"

/**
	@author root <root@localhost.localdomain>
*/
class HMMInitializer {

	private:	
		
		int m_iDim;
		int m_iCovarianceModelling;
		PhoneSet *m_phoneSet;
		int m_iState;
		
		// global distribution
		Vector<float> *m_vMeanGlobal;
		SMatrix<float> *m_mCovarianceGlobal;

		// compute the global mean and covariance
		void computeGlobalDistribution(VMLFUtterance *vMLFUtterance, const char *strFolderFeatures);
	
	public:
	
		// constructor
		HMMInitializer(int iDim, int iCovarianceModelling, PhoneSet *phoneSet);

		// destructor
		~HMMInitializer();	
		
		// initialize the HMMs to have the mean and covariance of the global set of features		
		HMMManager *initializeModelsFlatStart(VMLFUtterance *vMLFUtterance, const char *strFolderFeatures);
		
		// initialize the HMMs from baseline HMMs
		HMMManager *initializeModelsBootstrap(VMLFUtterance *vMLFUtterance, const char *strFolderFeatures, const char *strBootstrapModelsFile, bool bFloorCovariance, float fCovarianceFloor, float *fMeanGlobal, float *fCovarianceGlobal);
	
};

};	// end-of-namespace

#endif
