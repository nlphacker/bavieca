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


#include <iostream>
#include <iomanip>

#include "FileUtils.h"
#include "HLDAEstimator.h"
#include "HMMManager.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "TimeUtils.h"
#include "Transform.h"

#include "Vector.h"
#include "VectorStatic.h"
#include "Matrix.h"
#include "MatrixStatic.h"

namespace Bavieca {

// constructor
HLDAEstimator::HLDAEstimator(HMMManager *hmmManager, const char *strFileAccList, unsigned int iDimensionalityReduction, int iIterationsTransformUpdate, int iIterationsParameterUpdate, const char *strFolderOutput)
{
	m_hmmManager = hmmManager;
	m_mCovarianceGlobal = NULL;
	m_strFileAccList = strFileAccList;
	m_strFolderOutput = strFolderOutput;
	
	// check parameters
	assert(iDimensionalityReduction > 0);
	assert(iDimensionalityReduction <= m_hmmManager->getFeatureDim());
	m_iN = hmmManager->getFeatureDim();
	m_iP = m_iN-iDimensionalityReduction;
	
	// estimation iterations
	m_iIterationsTransformUpdate = iIterationsTransformUpdate;
	m_iIterationsParametersUpdate = iIterationsParameterUpdate;
}

// destructor
HLDAEstimator::~HLDAEstimator()
{
	if (m_mCovarianceGlobal != NULL) {
		delete m_mCovarianceGlobal;
	}
}

// estimate the HLDA transform
// - the HLDA transform decorrelates the features and reduces their dimensionality
// - the input is a set of models on the extended feature space (N dimensions) with full covariance matrix,
//   occupation counts for each Gaussian are in the accumulators
// - the output is a linear transform in the form of a (NxN)
// - each Gaussian component defines a class, features aligned to that Gaussian component belong to that class
void HLDAEstimator::estimate() {
	
	// get the number of Gaussian components in the system
	int iGaussianComponents = m_hmmManager->getNumberGaussianComponents();
		
	// get the HMM-states in the system
	int iHMMStates = -1;
	HMMState **hmmStates = m_hmmManager->getHMMStates(&iHMMStates);	
	// allocate memory to store the Gaussian components and get the total occupation
	double dOccupationTotal = 0.0;	
	GaussianData *gaussianData = new GaussianData[iGaussianComponents];
	int iGaussian = 0;
	for(int i=0 ; i < iHMMStates ; ++i) {
		for(unsigned int j=0 ; j < hmmStates[i]->getMixture().getNumberComponents() ; ++j, ++iGaussian) {
			gaussianData[iGaussian].gaussian = hmmStates[i]->getMixture()(j);
			gaussianData[iGaussian].vCovDiag = new Vector<float>(m_iN);
			dOccupationTotal += gaussianData[iGaussian].dOccupation;
		}
	}
	
	// load the Gaussian occupation and compute the global covariance
	loadGaussianOccupation(&dOccupationTotal,gaussianData);
	
	// cofactor
	Matrix<float> matrixA(m_iN);
	Matrix<float> matrixACofactor(m_iN);
	// auxiliar vector
	Vector<float> vAux(m_iN);
	// G matrix and its invert
	Matrix<float> matrixG(m_iN);
	Matrix<float> matrixGInverted(m_iN);
	
	// initialize the transform to the identity matrix
	matrixA.setIdentity();
		
	// compute the initial global likelihood
	float fLikelihoodBaseline = computeLikelihood(matrixA,gaussianData,iGaussianComponents);
	float fLikelihoodFrame = (float)(fLikelihoodBaseline/dOccupationTotal);
	printf("baseline models:         likelihood= %.4f (%6.2f)\n",fLikelihoodBaseline,fLikelihoodFrame);	
	
	// actual iterative re-estimation	
	for(int iIterationOut = 1 ; iIterationOut <= m_iIterationsParametersUpdate ; ++iIterationOut) {
	
		printf("iteration: %2d ",iIterationOut);
		fflush(stdout);	
		float fPercentageDisplayed = 0.0;
		
	   double dBegin = TimeUtils::getTimeMilliseconds();		
		
		Matrix<float> matrixCovGlobalFull(*m_mCovarianceGlobal);
	
		// compute the diagonal covariance estimates for every Gaussian component
		for(int g=0 ; g < iGaussianComponents ; ++g) {
			Matrix<float> matrixCovFull(gaussianData[g].gaussian->covarianceFull());
			// compute each element of the diagonal covariance matrix
			for(int k=0 ; k<m_iN ; ++k) {
				vAux.zero();
				// first P rows
				if (k < m_iP) {
					vAux.mul(matrixA.getRow(k),matrixCovFull);
				}
				// last N-P rows
				else {
					vAux.mul(matrixA.getRow(k),matrixCovGlobalFull);
				}
				(*gaussianData[g].vCovDiag)(k) = vAux.mul(matrixA.getRow(k));
				assert((*gaussianData[g].vCovDiag)(k) > 0);
			}	
		}	
		
		// iterative transform update	
		for(int iIterationIn = 1 ; iIterationIn <= m_iIterationsTransformUpdate ; ++iIterationIn) {
		
			// compute the cofactor of the current estimate of A
			matrixACofactor.copy(matrixA);
			matrixACofactor.matrixCofactor();
		
			// for each row:
			for(int r = 0 ; r < m_iN ; ++r) {
				
				// (1) compute the G matrix (n x n) associated to the r-th row
				matrixG.zero();
				for(int g=0 ; g < iGaussianComponents ; ++g) {
					
					Matrix<float> matrixCovFull(gaussianData[g].gaussian->covarianceFull());
					double dConstant = gaussianData[g].dOccupation/gaussianData[g].vCovDiag->getData()[r];
					// first P rows
					if (r < m_iP) {
						matrixG.add((float)dConstant,matrixCovFull);
					}
					// last (N-P) rows
					else {
						matrixG.add((float)dConstant,matrixCovGlobalFull);
					}
				}	
		
				// (2) invert G
				matrixGInverted.copy(matrixG);
				matrixGInverted.invert();
			
				// (3) compute the row in A
				matrixA.getRow(r).mul(matrixACofactor.getRow(r),matrixGInverted);
				double dAux = matrixA.getRow(r).mul(matrixACofactor.getRow(r));	
				matrixA.getRow(r).mul((float)sqrt(dOccupationTotal/dAux));
			}

			// update the progress bar if necessary
			float fPercentage = (((float)iIterationIn)*100.0f)/((float)m_iIterationsTransformUpdate);
			if (fPercentage >= fPercentageDisplayed + 10.0) {
				fPercentageDisplayed += 10.0;
				cout << "*";
				fflush(stdout);
			}
		}
		// update the progress bar if necessary
		if (fPercentageDisplayed < 100.0) {
			cout << "*";
		}
		
		// get the iteration end time
		double dEnd = TimeUtils::getTimeMilliseconds();
		double dMillisecondsInterval = dEnd - dBegin;	
		
		float fLikelihoodGlobal = computeLikelihood(matrixA,gaussianData,iGaussianComponents);
		float fLikelihoodFrame = (float)(fLikelihoodGlobal/dOccupationTotal);
		
		// compute the Real Time Factor of the reestimation process
		float fRTF = (float)((dMillisecondsInterval/10.0f)/((float)dOccupationTotal));	
		// compute audio used for the estimation
		int iFeatureVectorsPerSecond = 100;
		int iFeatureVectorsPerMinute = 60*100;
		int iFeatureVectorsPerHour = 60*60*100;
		int iHours = ((int)dOccupationTotal)/iFeatureVectorsPerHour;
		int iMinutes = ((int)dOccupationTotal)%iFeatureVectorsPerHour;
		iMinutes /= iFeatureVectorsPerMinute;
		int iSeconds = ((int)dOccupationTotal)%iFeatureVectorsPerMinute;
		iSeconds /= iFeatureVectorsPerSecond;	
		// show the iteration information
		printf(" likelihood= %.4f (%6.2f) [%8d Gauss][RTF=%.6f][%d:%02d'%02d''][%6.2f%%]\n",
			fLikelihoodGlobal,fLikelihoodFrame,iGaussianComponents,fRTF,iHours,iMinutes,
			iSeconds,100.0*(1.0-(fLikelihoodGlobal/fLikelihoodBaseline)));	
		
		// normalize the transform so it has determinant equal to one
		Matrix<float> matrixAInverted(matrixA);
		float fDet = matrixAInverted.invert();
		Matrix<float> matrixANormalized(matrixA);
		matrixANormalized.mul(1.0f/pow(fDet,(float)(1.0f/((float)m_iN))));
		/*matrixAInverted.copy(matrixANormalized);
		fDet = matrixAInverted.invert();
		printf("determinant |A|= %f\n",fDet);*/
		
		// dump the transform to disk
		ostringstream strFileTransform;
		strFileTransform << m_strFolderOutput << PATH_SEPARATOR << "transform" << setw(2) << iIterationOut << ".bin";
		Transform *transform = new Transform(TRANSFORM_TYPE_LINEAR,matrixANormalized);
		transform->store(strFileTransform.str().c_str());
		delete transform;
	}
}

// compute the likelihood for the given set of Gaussian distributions with associated occupations
float HLDAEstimator::computeLikelihood(Matrix<float> &matrixA, GaussianData *gaussianData, int iGaussians) {

	float fLikelihood = 0.0;
	
	Vector<float> vMeanX(m_iP);					// only the first P dimensions are kept
	Vector<float> vCovarianceX(m_iP);			// only the first P dimensions are kept
	Matrix<float> matrixCovDiag(m_iP);
	MatrixStatic<float> matrixAReduced(matrixA,0,m_iP,0,m_iN);
		
	for(int g = 0 ; g < iGaussians ; ++g) {

		// compute the transformed mean: A[p] * mean
		vMeanX.mul(matrixAReduced,gaussianData[g].gaussian->mean());

		// compute transformed covariance
		Matrix<float> matrixCovFull(gaussianData[g].gaussian->covarianceFull());
		matrixCovDiag.zero();
		matrixCovDiag.addMul(matrixAReduced,no,matrixCovFull,no,matrixAReduced,yes);		// A[p] * COV * A[p]^T 
		matrixCovDiag.getDiagonal(vCovarianceX);
		assert(vCovarianceX.isPositive());
		
		// compute likelihood
		fLikelihood += (float)computeLikelihood(vMeanX,vCovarianceX,(float)gaussianData[g].dOccupation,m_iP);
	}
	
	return fLikelihood;
}

// compute the likelihood given the sufficient statistics
float HLDAEstimator::computeLikelihood(VectorBase<float> &vMean, VectorBase<float> &vCovariance, 
	float fOccupation, int iDimensionality) {
	
	assert(vCovariance.isPositive());
	double dDeterminant = 0.0;
	for(int i=0 ; i < iDimensionality ; ++i) {
		dDeterminant += log(vCovariance(i));
	}
	assert(finite(dDeterminant));
	float fLikelihood = (float)((-0.5*fOccupation)*(dDeterminant+(iDimensionality*(log(2*PI_NUMBER)+1))));	
	assert(finite(fLikelihood));
	
	return fLikelihood;
}

// apply the given HLDA transform to the given acoustic models
void HLDAEstimator::applyTransform(Matrix<float> &mTransform, HMMManager *hmmManager) {

	unsigned int iP = mTransform.getRows();
	unsigned int iN = mTransform.getCols();
	int iHMMStates = -1;
	HMMState **hmmStates = hmmManager->getHMMStates(&iHMMStates);	

	// (1) checks
	assert(iP <= iN);
	assert(iN == hmmManager->getFeatureDim());
	assert(hmmManager->getCovarianceModelling() != COVARIANCE_MODELLING_TYPE_FULL);
	
	// (2) transform model parameters
	Matrix<float> matrixCovDiag(iP);
	
	// update each state
	for(unsigned int s = 0 ; s < (unsigned int)iHMMStates ; ++s) {

		// update each Gaussian component
		for(unsigned int g = 0 ; g < hmmStates[s]->getMixture().getNumberComponents() ; ++g) {
			
			Gaussian *gaussian = hmmStates[s]->getMixture()(g);
			assert(gaussian != NULL);
			
			// compute the transformed mean: A[p] * mean
			Vector<float> vMeanX(iP);
			vMeanX.mul(mTransform,gaussian->mean());
	
			// compute transformed covariance
			Matrix<float> matrixCovFull(gaussian->covarianceFull());	
			matrixCovDiag.zero();
			matrixCovDiag.addMul(mTransform,no,matrixCovFull,no,mTransform,yes);		// A[p] * COV * A[p]^T 
			Vector<float> vCovarianceX(iP);
			matrixCovDiag.getDiagonal(vCovarianceX);
			assert(vCovarianceX.isPositive());
			
			// replace the original mean and covariance
			gaussian->setMeanCov(vMeanX,vCovarianceX);
		}
		
		// update the model configuration
		hmmStates[s]->setFeatureDimensionality(iP);
		hmmStates[s]->setCovarianceModelling(COVARIANCE_MODELLING_TYPE_DIAGONAL);
	}
	
	// update the models configuration
	hmmManager->setFeatureDimensionality(iP);
	hmmManager->setCovarianceModelling(COVARIANCE_MODELLING_TYPE_DIAGONAL);
}

// apply HLDA transform to features
MatrixBase<float> *HLDAEstimator::applyTransform(MatrixBase<float> &mTransform, MatrixBase<float> &mFeatures) {

	// check parameters
	assert((mTransform.getRows() <= mFeatures.getCols()) && (mTransform.getCols() == mFeatures.getCols()));
	
	// transformed features
	Matrix<float> *mFeaturesX = new Matrix<float>(mFeatures.getRows(),mTransform.getRows());
	// transform each feature vector
	for(unsigned int i=0 ; i<mFeatures.getRows() ; ++i) {
		mFeaturesX->getRow(i).mul(mTransform,mFeatures.getRow(i));
	}
	
	// transform the features
	/*for(int f=0 ; f<iFeatures ; ++f) {
		for(int i=0 ; i<iP ; ++i) {
			(*fFeaturesX)[f*iP+i] = 0.0f;
			for(int j=0 ; j<iN ; ++j) {
				(*fFeaturesX)[f*iP+i] += fA[i*iN+j]*fFeatures[f*iN+j];
			}
		}
	}*/

	return mFeaturesX;
}

// load the Gaussian occupation from the accumulators
void HLDAEstimator::loadGaussianOccupation(double *dOccupationTotal, GaussianData *gaussianData) {

	// 1. load the accumulators
	MAccumulatorPhysical mAccumulatorPhysical;
	AccMetadata metadata;
	Accumulator::loadAccumulatorList(m_strFileAccList,mAccumulatorPhysical,metadata);
		
	// get the HMM-states in the system
	int iHMMStatesAux = -1;
	HMMState **hmmStates = m_hmmManager->getHMMStates(&iHMMStatesAux);
	
	// check consistency
	if ((iHMMStatesAux != metadata.iHMMStates) || 
		(m_hmmManager->getNumberGaussianComponents() != metadata.iGaussianComponents) ||
		(m_hmmManager->getCovarianceModelling() != metadata.iCovarianceModeling) ||
		(m_hmmManager->getFeatureDim() != (unsigned int)metadata.iDim)) {
		BVC_ERROR << "HMMs are not consistent with accumulators";
	}
		
	int iDim = m_hmmManager->getFeatureDim();
	*dOccupationTotal = 0.0;
	
	int iGaussian = 0;
	for(int i=0 ; i < m_hmmManager->getNumberHMMStatesPhysical() ; ++i) {
		for(unsigned int j=0 ; j < hmmStates[i]->getMixture().getNumberComponents() ; ++j, ++iGaussian) {
			unsigned int iKey = Accumulator::getPhysicalAccumulatorKey(i,j);
			MAccumulatorPhysical::iterator it = mAccumulatorPhysical.find(iKey);
			assert(it != mAccumulatorPhysical.end());
			Accumulator *accumulator = it->second;
			gaussianData[iGaussian].dOccupation = accumulator->getOccupation();
			*dOccupationTotal += accumulator->getOccupation();
		}
	}	
	
	cout << "total occupation: " << *dOccupationTotal;
	cout << "# Gaussian components: " << metadata.iGaussianComponents;
	cout << "# Gaussian components with occupation: " << mAccumulatorPhysical.size();
	
	// compute the global distribution
	Vector<double> vGlobalMean(iDim);
	SMatrix<double> mGlobalCovariance(iDim);
	SMatrix<double> mObservationSquare(iDim);
	
	// compute the global mean
	for(MAccumulatorPhysical::iterator it = mAccumulatorPhysical.begin() ; it != mAccumulatorPhysical.end() ; ++it) {
		Accumulator *accumulator = it->second;
		vGlobalMean.add(accumulator->getObservation());
	}
	vGlobalMean.mul(1.0/(*dOccupationTotal));
	
	// compute global covariance
	for(MAccumulatorPhysical::iterator it = mAccumulatorPhysical.begin() ; it != mAccumulatorPhysical.end() ; ++it) {
		Accumulator *accumulator = it->second;
		mObservationSquare.add(accumulator->getObservationSquareFull());
	}	
	
	mGlobalCovariance.add(1.0/(*dOccupationTotal),mObservationSquare);
	mGlobalCovariance.addSquare(1.0,vGlobalMean);
	m_mCovarianceGlobal = new SMatrix<float>(mGlobalCovariance);
	
	cout << "global covariance computed" << endl;
	
	Accumulator::destroy(mAccumulatorPhysical);
}

};	// end-of-namespace
 
 

