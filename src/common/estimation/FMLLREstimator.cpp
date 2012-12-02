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


#include "AlignmentFile.h"
#include "BatchFile.h"
#include "FMLLREstimator.h"
#include "FeatureFile.h"
#include "HMMManager.h"
#include "LogMessage.h"
#include "PhoneSet.h"
#include "Transform.h"

namespace Bavieca {

// constructor
FMLLREstimator::FMLLREstimator(PhoneSet *phoneSet, HMMManager *hmmManager, int iIterations, bool bBestComponentOnly)
{
	m_phoneSet = phoneSet;
	m_hmmManager = hmmManager;
	m_iIterations = iIterations;
	m_bBestComponentOnly = bBestComponentOnly;
	m_iDim = hmmManager->getFeatureDimensionality();
	m_fOccupancyTotal = 0;
}

// destructor
FMLLREstimator::~FMLLREstimator()
{
}

// initialize the estimation
void FMLLREstimator::initializeEstimation() {
	
	m_matrixG = new Matrix<double>*[m_iDim];
	for(int i=0 ; i<m_iDim ; ++i) {
		m_matrixG[i] = new Matrix<double>(m_iDim+1);
		m_matrixG[i]->zero();
	}
	m_matrixK = new Matrix<double>(m_iDim,m_iDim+1);
	m_matrixK->zero();
	m_matrixAux = new Matrix<double>(m_iDim+1);
}

// uninitialize the estimation
void FMLLREstimator::uninitializeEstimation() {
	
	for(int i=0 ; i<m_iDim ; ++i) {
		delete m_matrixG[i];
	}
	delete [] m_matrixG;
	delete m_matrixK;
	delete m_matrixAux;
}

// feed adaptation data from an alignment into the adaptation process
// note: it is possible to accumulate statistics at the Gaussian component or at the transform level
// - Gaussian component level: O(n^2) per Gaussian component
// - Transform level: O(n^3) per transform + O(n) per Gaussian component
void FMLLREstimator::feedAdaptationData(float *fFeatures, int iFeatures, Alignment *alignment, 
	double *dLikelihood) {

	// sanity check
	assert(iFeatures == alignment->getFrames());

	m_fOccupancyTotal += iFeatures;

	*dLikelihood = 0.0;
	for(int t=0 ; t<iFeatures ; ++t) {
		
		float *fFeatureVector = fFeatures+(t*m_iDim);
		
		// extended observation vector
		Vector<double> vObsEx(VectorStatic<float>(fFeatureVector,m_iDim));
		vObsEx.appendFront(1.0);	
		
		// for each HMM-state the observation is assigned to
		FrameAlignment *frameAlignment = alignment->getFrameAlignment(t);
		for(FrameAlignment::iterator it = frameAlignment->begin() ; it != frameAlignment->end() ; ++it) {
			HMMStateDecoding *hmmStateDecoding = m_hmmManager->getHMMStateDecoding((*it)->iHMMState);
			// compute the contribution of each Gaussian 
			// (case 1) all the frame-level adaptation data goes to the best scoring Gaussian component (faster)
			if (m_bBestComponentOnly) {
				float fLikelihood = -FLT_MAX;
				GaussianDecoding *gaussian = hmmStateDecoding->getBestScoringGaussian(fFeatureVector,&fLikelihood);
				*dLikelihood += std::max<float>(fLikelihood,LOG_LIKELIHOOD_FLOOR);
					
				Vector<float> vCovariance(VectorStatic<float>(gaussian->fCovariance,m_iDim));		
			#ifdef OPTIMIZED_COMPUTATION
				vCovariance.mul(2.0);
				vCovariance.invertElements();
			#endif

				// accumulate data for each G(i)
				m_matrixAux->zero();
				m_matrixAux->addVecMul(1.0,vObsEx,vObsEx);
				for(int i=0 ; i < m_iDim; ++i) {
					double dCovInv = 1.0/vCovariance(i);
					m_matrixG[i]->add(dCovInv*1.0,*m_matrixAux);
				}
				
				// accumulate data for each k(i)
				for(int i=0 ; i < m_iDim ; ++i) {	
					double dConstant = (gaussian->fMean[i]/vCovariance(i))*1.0;
					m_matrixK->getRow(i).add(dConstant,vObsEx);
				}	
			}
			// (case 2) adaptation data is shared across all components (slightly more accurate)
			else {			
				double *dProbGaussian = new double[hmmStateDecoding->getGaussianComponents()];
				double dProbTotal = 0.0;
				for(int g=0 ; g < hmmStateDecoding->getGaussianComponents() ; ++g) {	
					dProbGaussian[g] = exp(hmmStateDecoding->computeGaussianProbability(g,fFeatureVector));
					dProbTotal += dProbGaussian[g];
					assert((finite(dProbGaussian[g])) && (finite(dProbTotal)));	
				}
				*dLikelihood += max(log(dProbTotal),LOG_LIKELIHOOD_FLOOR);
				for(unsigned int iGaussian = 0 ; iGaussian < hmmStateDecoding->getMixtureSize() ; ++iGaussian) {
				
					GaussianDecoding *gaussian = hmmStateDecoding->getGaussian(iGaussian);	
					double dGaussianOccupation = dProbGaussian[iGaussian]/dProbTotal;
					
					Vector<float> vCovariance(VectorStatic<float>(gaussian->fCovariance,m_iDim));		
				#ifdef OPTIMIZED_COMPUTATION
					vCovariance.mul(2.0);
					vCovariance.invertElements();
				#endif
	
					// accumulate data for each G(i)
					m_matrixAux->zero();
					m_matrixAux->addVecMul(1.0,vObsEx,vObsEx);
					for(int i=0 ; i < m_iDim; ++i) {
						double dCovInv = 1.0/vCovariance(i);
						m_matrixG[i]->add(dCovInv*dGaussianOccupation,*m_matrixAux);
					}
					
					// accumulate data for each k(i)
					for(int i=0 ; i < m_iDim ; ++i) {	
						double dConstant = (gaussian->fMean[i]/vCovariance(i))*dGaussianOccupation;
						m_matrixK->getRow(i).add(dConstant,vObsEx);
					}
				}
				delete [] dProbGaussian;
			}
		}
	}
}

// feed adaptation data from a batch file containing entries (rawFile alignmentFile)
void FMLLREstimator::feedAdaptationData(const char *strBatchFile, const char *strAlignmentFormat, 
	double *dLikelihood, bool bVerbose) {

	BatchFile batchFile(strBatchFile,"features|alignment");
	batchFile.load();
	
	for(unsigned int i=0 ; i < batchFile.size() ; ++i) {
	//for(int i=0 ; i < 5 ; ++i) {
		
		// load the alignment
		Alignment *alignment = NULL;
		if (strcmp(strAlignmentFormat,"text") == 0) {
			AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet);	
			VPhoneAlignment *vPhoneAlignment = alignmentFile->load(batchFile.getField(i,"alignment"));
			assert(vPhoneAlignment);
			alignment = AlignmentFile::toAlignment(m_phoneSet,m_hmmManager,vPhoneAlignment);
			AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
			delete alignmentFile;
		} else {
			alignment = Alignment::load(batchFile.getField(i,"alignment"),NULL);
			assert(alignment);	
		}
		
		// load the feature vectors
		FeatureFile *featureFile = new FeatureFile(batchFile.getField(i,"features"),MODE_READ);
		featureFile->load();
		int iFeatureVectors = -1;
		float *fFeatures = featureFile->getFeatureVectors(&iFeatureVectors);
		
		// load and apply the transform
		/*Transform *transform = new Transform();
		transform->load("/home/speech/wsj/scripts/test/fmllr/transform.bin");
		float *fFeaturesX = new float[iFeatureVectors*m_iDim];
		for(int j=0;j<iFeatureVectors;++j) {
			transform->apply(fFeatures+j*m_iDim,fFeaturesX+j*m_iDim);
		}
		fFeatures = fFeaturesX;
		delete transform;*/
		
		// check consistency
		if (iFeatureVectors != alignment->getFrames()) {
			BVC_ERROR << "inconsistent number of feature vectors / alignment file";
		}
		
		// accumulate adaptation data
		double dLikelihoodAlignment = 0.0;
		feedAdaptationData(fFeatures,iFeatureVectors,alignment,&dLikelihoodAlignment);
		if (bVerbose) {
			printf("loaded file: %s likelihood: %10.2f (%6d frames)\n",batchFile.getField(i,"alignment"),dLikelihoodAlignment,iFeatureVectors);
		}
		*dLikelihood += dLikelihoodAlignment;
		
		// clean-up
		delete alignment;
		delete featureFile;
		delete [] fFeatures;		
	}
	if (bVerbose) {
		double dLikelihoodFrame = (*dLikelihood)/m_fOccupancyTotal;
		printf("total likelihood: %20.6f (likelihood per frame: %8.4f)\n",*dLikelihood,dLikelihoodFrame);
	}
}

// estimate the feature transform for the given data (typically speaker adaptation data)
Transform *FMLLREstimator::estimateTransform(Transform *transformInitial) {

	Matrix<double> matrixW(m_iDim,m_iDim+1);
	Matrix<float> matrixWF(m_iDim,m_iDim+1);
	Matrix<double> matrixA(m_iDim,m_iDim);
	Vector<double> vAux(m_iDim+1);
	
	// starting from scratch: bias = 0.0 and A = identity
	if (transformInitial == NULL) {
		matrixA.setIdentity();
	} 
	// starting from another transform
	else {
		assert(0);
	}	
	
	// compute the matrix of cofactors of A
	Matrix<double> matrixP(matrixA);
	matrixP.matrixCofactor();
	
	// initialize the transform
	
	// iterative update
	for(int iIteration = 0 ; iIteration < m_iIterations ; ++iIteration) {
	
		// compute the transform row by row
		for(int i=0 ; i < m_iDim ; ++i) {
		
			// invert G(i) (this is suboptimal, the inversion needs to be
			// done just once then it can be reused across iterations)
			Matrix<double> matrixGInverted(*m_matrixG[i]);
			matrixGInverted.invert();
			
			// compute alpha (it implies solving a quadratic equation and selecting the solution
			// that produces a higher increase of the objective function)
			
			// get the zero extended vector of cofactors for the ith row
			Vector<double> vCofactor(matrixP.getRow(i));
			vCofactor.appendFront(0.0);
			
			// compute constants
			double dBeta = m_fOccupancyTotal;
			Vector<double> vAux(m_iDim+1);
			vAux.mul(vCofactor,matrixGInverted);
			double dA = vAux.mul(vCofactor);	
			double dB = vAux.mul(m_matrixK->getRow(i));
			double dC = -1.0*dBeta;
			
			// there are two possible solutions
			double dAlpha1 = (-1.0*dB + sqrt(dB*dB - 4.0*dA*dC))/(2.0*dA);
			double dAlpha2 = (-1.0*dB - sqrt(dB*dB - 4.0*dA*dC))/(2.0*dA);
			
			// we let the objective function decide which solution is better (higher likelihood increase)
			double dAuxiliarFunction1 = dBeta*log(fabs(dAlpha1*dA+dB)) - 0.5*dAlpha1*dAlpha1*dA;	
			double dAuxiliarFunction2 = dBeta*log(fabs(dAlpha2*dA+dB)) - 0.5*dAlpha2*dAlpha2*dA;
			double dAlpha = (dAuxiliarFunction1 > dAuxiliarFunction2) ? dAlpha1 : dAlpha2;
			
			// compute the w(i)
			Vector<double> vAux2(m_matrixK->getRow(i));
			Vector<double> vCofactorRowEx(matrixP.getRow(i));
			vCofactorRowEx.appendFront(0.0);
			vAux2.add(dAlpha,vCofactorRowEx);
			matrixW.getRow(i).mul(vAux2,matrixGInverted);
			
			// update A
			matrixA.getRow(i).copy(VectorStatic<double>(matrixW.getRowData(i)+1,m_iDim));
			
			// compute the matrix of cofactors after updating the row of the transform
			matrixP.copy(matrixA);
			matrixP.matrixCofactor();	
		}
		
		// compute the log(|A|)
		Matrix<double> matrixAInverted(matrixA);	
		float fDet = (float)matrixAInverted.invert();
		if (fDet == 0) {
			BVC_ERROR << "the A matrix is singular!, can't compute inverse";
		}
		printf("log(|A|) = %8.4f\n",log(fDet));	
	}
	
	// create the transform
	Matrix<float> matrixWFloat(matrixW);
	return new Transform(TRANSFORM_TYPE_AFFINE,matrixWFloat);
}

};	// end-of-namespace
