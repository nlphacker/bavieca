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


#ifndef HLDAESTIMATOR_H
#define HLDAESTIMATOR_H

using namespace std;

#include <string>

#include "Matrix.h"
#include "MatrixStatic.h"
#include "Vector.h"
#include "VectorStatic.h"

namespace Bavieca {

class Gaussian;
class HMMManager;
class LexiconManager;
class PhoneSet;
class Transform;

// auxiliar structure for HLDA computation
typedef struct {
	Vector<float> *vCovDiag;			// diagonal covariance
	Gaussian *gaussian;					// Gaussian component (full covariance)
	double dOccupation;					// Gaussian occupation
} GaussianData;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class HLDAEstimator {

	private: 
	
		HMMManager *m_hmmManager;					// acoustic models
		SMatrix<float> *m_mCovarianceGlobal;	// global covariance (all data used for the estimation)
		const char *m_strFileAccList;				// Gaussian occupation
		string m_strFolderOutput;					// output folder
		int m_iN;										// original dimensionality
		int m_iP;										// reduced dimensionality
		
		// iterations
		int m_iIterationsTransformUpdate;
		int m_iIterationsParametersUpdate;
		
		// compute the likelihood for the given set of Gaussian distributions with associated occupations
		float computeLikelihood(Matrix<float> &matrixA, GaussianData *gaussianData, int iGaussians);
		
		// compute the likelihood given the sufficient statistics
		static float computeLikelihood(VectorBase<float> &vMean, VectorBase<float> &vCovariance, 
			float fOccupation, int iDimensionality);
		
		// load the Gaussian occupation from the accumulators
		void loadGaussianOccupation(double *dOccupationTotal, GaussianData *gaussianData);	

	public:	

		// constructor
		HLDAEstimator(HMMManager *hmmManager, const char *strFileAccList, unsigned int iDimensionalityReduction, int iIterationsTransformUpdate, int iIterationsParameterUpdate, const char *strFolderOutput);

		// destructor
		~HLDAEstimator();
		
		// estimate the HLDA transform
		void estimate();
		
		// apply the given HLDA transform to the given acoustic models
		static void applyTransform(Matrix<float> &mTransform, HMMManager *hmmManager);

		// apply HLDA transform to features
		static MatrixBase<float> *applyTransform(MatrixBase<float> &mTransform, MatrixBase<float> &mFeatures);

};

};	// end-of-namespace

#endif
