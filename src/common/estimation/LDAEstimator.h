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


#ifndef LDAESTIMATOR_H
#define LDAESTIMATOR_H

using namespace std;

#include <string>

#include "Vector.h"
#include "SMatrix.h"
#include "TMatrix.h"
#include "Matrix.h"

namespace Bavieca {

class HMMManager;
class LexiconManager;
class PhoneSet;
class Transform;

typedef struct {
	double dOcc;
	//Vector<double> vObs;
} ClassStats;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class LDAEstimator {

	private:
	
		int m_iDimSource;							// original dimensionality	
		int m_iDimReduced;						// reduced dimensionality
		int m_iClasses;							// classes
		
		ClassStats *m_stats;
		SMatrix<double> *m_smObsSquare;
		
	public:

		// constructor
		LDAEstimator(int iDimSource, int iDimReduced, int iClasses);

		// destructor
		~LDAEstimator();
		
		// initialization
		void initialize();
		
		// accumulate data
		void accumulate(Vector<double> &vObs, SMatrix<double> &smObsSquare, 
			double dOcc, int iClass);	
		
		// estimate the transform
		void estimate();

};

};	// end-of-namespace

#endif
