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


#ifndef FMLLRESTIMATOR_H
#define FMLLRESTIMATOR_H

#include "Matrix.h"
#include "Vector.h"
#include "MatrixStatic.h"
#include "VectorStatic.h"

namespace Bavieca {

class Alignment;
class PhoneSet;
class HMMManager;
class Transform;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FMLLREstimator {

	private:
	
		PhoneSet *m_phoneSet;
		HMMManager *m_hmmManager;
		int m_iIterations;	
		bool m_bBestComponentOnly;
		int m_iDim;								// feature dimensionality
		float m_fOccupancyTotal;
		
		// estimation data
		Matrix<double> **m_matrixG;
		Matrix<double> *m_matrixK;
		Matrix<double> *m_matrixAux;

	public:

		// constructor
		FMLLREstimator(PhoneSet *phoneSet, HMMManager *hmmManager, int iIterations, bool bBestComponentOnly);

		// destructor
		~FMLLREstimator();
		
		// initialize the estimation
		void initializeEstimation();
		
		// uninitialize the estimation
		void uninitializeEstimation();
		
		// feed adaptation data from an alignment into the adaptation process
		void feedAdaptationData(MatrixBase<float> &mFeatures, Alignment *alignment, double *dLikelihood);
			
		// feed adaptation data from a batch file containing entries (rawFile alignmentFile)
		void feedAdaptationData(const char *strBatchFile, const char *strAlignmentFormat, double *dLikelihood);
		
		// estimate the feature transform for the given data (typically speaker adaptation data)
		Transform *estimateTransform(Transform *transformInitial);
};

};	// end-of-namespace

#endif
