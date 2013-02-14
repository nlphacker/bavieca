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


#include "GMMEditor.h"
#include "MLEstimator.h"

namespace Bavieca {

// constructor
GMMEditor::GMMEditor(HMMManager *hmmManager, MAccumulatorPhysical *mAccumulatorPhysical)
{
	m_hmmManager = hmmManager;
	m_mAccumulatorPhysical = mAccumulatorPhysical;
	m_iDim = hmmManager->getFeatureDimensionality();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iDim,m_iCovarianceModeling);
	m_iHMMStates = -1;
	m_hmmStates = hmmManager->getHMMStates(&m_iHMMStates);
}

// destructor
GMMEditor::~GMMEditor()
{
	if (m_mGaussianInfo != NULL) {
		delete [] m_mGaussianInfo;
	}
}

// initialize the editor
void GMMEditor::initialize() {

	m_mGaussianInfo = new MGaussianInfo[m_iHMMStates];
	int iHMMState = -1;
	int iGaussian = -1;
	for(MAccumulatorPhysical::iterator it = m_mAccumulatorPhysical->begin() ; it != m_mAccumulatorPhysical->end() ; ++it) {
		Accumulator::getPhysicalAccumulatorValues(it->first,iHMMState,iGaussian);
		GaussianInfo gaussianInfo;
		gaussianInfo.dOccupation = it->second->getOccupation();
		gaussianInfo.bMerged = false;
		m_mGaussianInfo[iHMMState].insert(MGaussianInfo::value_type(m_hmmStates[iHMMState]->getMixture()(iGaussian),gaussianInfo));
	}
}

// merge Gaussian components in the mixture
void GMMEditor::mixtureMerge(float fMinimumGaussianOccupation, float fMinimumGaussianWeight, float fCovarianceFlooringRatio) {

	// compute the covariance floor
	Vector<double> vCovarianceFloor(m_iDim);
	MLEstimator::computeCovarianceFloor(m_hmmManager,*m_mAccumulatorPhysical,
		fCovarianceFlooringRatio,vCovarianceFloor);
	
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		mixtureMerge(m_hmmStates[i],fMinimumGaussianOccupation,fMinimumGaussianWeight,vCovarianceFloor);
	}
}

// merge Gaussian components with not enough training data associated
// when a Gaussian is split into two gaussians, and the parameters of the newly created 
// Gaussian components are reestimated (iteration), it is possible that one of the gaussians 
// has very low occupation, in that case we need to merge that Gaussian to the closest one
void GMMEditor::mixtureMerge(HMMState *hmmState, float fMinimumGaussianOccupation, 
	float fMinimumGaussianWeight, Vector<double> &vCovarianceFloor) {

	assert(hmmState->getMixture().getNumberComponents() >= 1);
	
	MGaussianInfo &mGaussianInfo = m_mGaussianInfo[hmmState->getId()];
	
	// check that there are multiple Gaussian components
	if (hmmState->getMixture().getNumberComponents() == 1) {
		mGaussianInfo[hmmState->getMixture()(0)].bMerged = false;	
		return;
	}
	
	// allocate memory for the accumulators
	Vector<float> vMean(m_iDim);
	Vector<float> vCovariance(m_iDim);
	
	// (1) find Gaussian components with occupation/weight below the threshold and merge them to the nearest one
	bool bMerged = false;
	do {	
		bMerged = false;	
		for(int g = 0 ; g < hmmState->getMixture().getNumberComponents() ; ++g) {
			Gaussian *gaussian1 = hmmState->getMixture()(g);
			GaussianInfo *gaussianInfo1 = &mGaussianInfo[gaussian1];
			if ((gaussianInfo1->dOccupation < fMinimumGaussianOccupation) || 
				(gaussian1->weight() < fMinimumGaussianWeight)) {
				int iGaussianClosest = -1;
				// case 1: a single Gaussian (it can happen after all the Gaussian components are merged to one)
				if (hmmState->getMixture().getNumberComponents() == 1) {
					// warning: the HMM-state has just one Gaussian with little data, it will be poorly estimated	
					break;
				}
				// case 2: two Gaussian components, there is no need to compute the distance, just merge them 
				else if (hmmState->getMixture().getNumberComponents() == 2) {
					iGaussianClosest = (g+1)%2;	
				}	
				// case 3: more than two Gaussian -> find the closest gaussian to the poorly trained one
				else {	
					// compute the determinant of the covariance matrix of the Gaussian to merge
					double dCovarianceDet1 = gaussian1->covarianceDiag().sumLog();
					// find the closest gaussian	
					double dDistanceSmallest = DBL_MAX;
					for(int h = 0 ; h < hmmState->getMixture().getNumberComponents() ; ++h) {
						if (h == g) {
							continue;
						}
						Gaussian *gaussian2 = hmmState->getMixture()(h);
						GaussianInfo *gaussianInfo2 = &mGaussianInfo[gaussian2];
						// compute the distance between the gaussians
						// determinant of the second gaussian
						double dCovarianceDet2 = gaussian2->covarianceDiag().sumLog();
						// compute the counts
						double dDenominator = gaussianInfo1->dOccupation+gaussianInfo2->dOccupation;
						double d1 = gaussianInfo1->dOccupation/dDenominator;
						double d2 = gaussianInfo2->dOccupation/dDenominator;
						// compute the determinant of the covariance matrix resulting from merging both gaussians
						double dCovarianceDetMerged = 0.0;
						for(int i = 0 ; i < m_iDim ; ++i) {	
							dCovarianceDetMerged += log((d1*gaussian1->covarianceDiag()(i))+(d2*gaussian2->covarianceDiag()(i))+(d1*d2*(   (gaussian1->mean()(i)-gaussian2->mean()(i))*(gaussian1->mean()(i)-gaussian2->mean()(i)))));
						}
						// compute the distance
						double dDistance = (((gaussianInfo1->dOccupation+gaussianInfo2->dOccupation)/2.0)*dCovarianceDetMerged) - ((gaussianInfo1->dOccupation/2.0)*dCovarianceDet1) - ((gaussianInfo2->dOccupation/2.0)*dCovarianceDet2);
						assert(finite(dDistance));
						if (dDistance <= dDistanceSmallest) {
							dDistanceSmallest = dDistance;
							iGaussianClosest = h;
						}
					}
				}
				// (3) merge the "less-robust" gaussian with the closest one
				assert(iGaussianClosest != -1);
				// merge the gaussians
				int h = iGaussianClosest;
				assert(g != h);
				Gaussian *gaussian2 = hmmState->getMixture()(h);
				GaussianInfo *gaussianInfo2 = &mGaussianInfo[gaussian2];
				// compute the counts
				double dDenominator = gaussianInfo1->dOccupation+gaussianInfo2->dOccupation;
				double d1 = gaussianInfo1->dOccupation/dDenominator;
				double d2 = gaussianInfo2->dOccupation/dDenominator;	
				// compute the new mean and covariance
				vMean.zero();		
				vMean.add((float)d1,gaussian1->mean());
				vMean.add((float)d2,gaussian2->mean());
				vCovariance.zero();
				vCovariance.add((float)d1,gaussian1->covarianceDiag());
				vCovariance.add((float)d2,gaussian2->covarianceDiag());
				Vector<float> vAux(gaussian1->mean());
				vAux.add(-1.0,gaussian2->mean());
				vAux.squareElements();
				vCovariance.add((float)(d1*d2),vAux);
				vCovariance.floor(vCovarianceFloor);
				// copy the mean and covariance
				gaussian1->mean().copy(vMean);
				gaussian1->covarianceDiag().copy(vCovariance);
				// compute the new weight
				gaussian1->weight() += gaussian2->weight();
				// accumulate occupancy (this is an iterative bottom-up clustering, occupancy might still be below threshold)
				gaussianInfo1->dOccupation += gaussianInfo2->dOccupation;
				// mark the new gaussian as "merged" so it wont be split next
				gaussianInfo1->bMerged = true;
				// remove the gaussian component
				hmmState->getMixture().removeGaussianComponent(h);
				mGaussianInfo.erase(mGaussianInfo.find(gaussian2));
				bMerged = true;
				break;
			}
		}
	} while(bMerged);
	
	assert((int)mGaussianInfo.size() == hmmState->getMixture().getNumberComponents());
	
	// sanity check: make sure all the gaussians have an occupancy and weight above the threshold
	if (hmmState->getMixture().getNumberComponents() > 1) {
		for(int g = 0 ; g < hmmState->getMixture().getNumberComponents() ; ++g) {
			Gaussian *gaussian = hmmState->getMixture()(g);
			GaussianInfo *gaussianInfo = &mGaussianInfo[gaussian];	
			assert(gaussian->weight() >= fMinimumGaussianWeight);
			assert(gaussianInfo->dOccupation >= fMinimumGaussianOccupation);
		}
	}
}

// add Gaussian components to each mixture
int GMMEditor::mixtureIncrement(int iIncrement, int iSplittingCriterion, 
	bool bAllowMultipleSplitsPerComponent, float fMinimumGaussianOccupation, float fSplittingEpsilon) {
	
	int iOriginal = 0;
	int iAdded = 0;
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		iOriginal += m_hmmStates[i]->getMixture().getNumberComponents();
		iAdded += mixtureIncrement(m_hmmStates[i],iIncrement,iSplittingCriterion,bAllowMultipleSplitsPerComponent,
			fMinimumGaussianOccupation,fSplittingEpsilon);
	}
	
	if (iAdded > 0) { 
		m_hmmManager->setSingleGaussian(false);
	}
	
	return iAdded;
}

// double the number of Gaussian components in the mixture (if possible)
int GMMEditor::mixtureDouble(float fMinimumGaussianOccupation, float fSplittingEpsilon) {
	
	int iOriginal = 0;
	int iAdded = 0;
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		iOriginal += m_hmmStates[i]->getMixture().getNumberComponents();
		iAdded += mixtureDouble(m_hmmStates[i],fMinimumGaussianOccupation,fSplittingEpsilon);
	}
	
	if (iAdded > 0) { 
		m_hmmManager->setSingleGaussian(false);
	}
	
	return iAdded;
}

// add a Gaussian to the mixture that models the emission probability of the state
// note: there is no need to floor the covariance of the new components, since they stay the same
int GMMEditor::mixtureIncrement(HMMState *hmmState, int iIncrement, int iSplittingCriterion, 
	bool bAllowMultipleSplitsPerComponent, float fMinimumGaussianOccupation, float fSplittingEpsilon) {

	assert(hmmState->getMixture().getNumberComponents() > 0);
	assert(iIncrement > 0);
		
	MGaussianInfo &mGaussianInfo = m_mGaussianInfo[hmmState->getId()];
	int iAdded = 0;			// number of components already added to the mixture
	
	// initialize splitting weights if needed
	if (iSplittingCriterion == SPLITTING_CRITERION_HEAVIEST_MIXTURE_COMPONENT) {
		for(int m = 0 ; m < hmmState->getMixture().getNumberComponents() ; ++m) {
			Gaussian *gaussian = hmmState->getMixture()(m);
			// the weight for mixture increment determines which Gaussian will be split next
			// a Gaussian that has already been split gets its weight substracted by one, what makes it
			// very unlikely to be split again
			gaussian->weightMixtureIncrement() = gaussian->weight();
		}	
	}
	
	while(iAdded < iIncrement) {
	
		int iGaussianToSplit = 0;
		// case 1: there are multiple gaussians: select the mixture to split
		if (hmmState->getMixture().getNumberComponents() > 1) {
			iGaussianToSplit = -1;
			switch(iSplittingCriterion) {
				case SPLITTING_CRITERION_LARGEST_AVERAGE_VARIANCE: {
					double dGeometricMean = 0.0;
					double dGeometricMeanHighest = -DBL_MAX;
					// for each gaussian compute the geometric mean of the variance
					for(int m = 0 ; m < hmmState->getMixture().getNumberComponents() ; ++m) {
						Gaussian *gaussian = hmmState->getMixture()(m);
						GaussianInfo *gaussianInfo = &mGaussianInfo[gaussian];
						// do not split gaussians that were just merged
						if (gaussianInfo->bMerged) {
							continue;
						}
						// the occupation must be at least two times the minimum gaussian occupation
						if (gaussianInfo->dOccupation < 2.0*fMinimumGaussianOccupation) {
							continue;
						}
						// compute the geometric mean of the covariance (log-average)
						dGeometricMean = gaussian->covarianceDiag().sumLog();
						dGeometricMean /= ((float)m_iDim);
						dGeometricMean = exp(dGeometricMean);
						// check if the geometic mean of the covariance is the largest one so far
						if (dGeometricMean > dGeometricMeanHighest) {
							dGeometricMeanHighest = dGeometricMean;
							iGaussianToSplit = m;
						}
					}
					break;
				}
				case SPLITTING_CRITERION_HEAVIEST_MIXTURE_COMPONENT: {
					// find the heaviest mixture
					float fWeightMax = -FLT_MAX;
					int iComponentHeaviest = -1;
					for(int m = 0 ; m < hmmState->getMixture().getNumberComponents() ; ++m) {
						Gaussian *gaussian = hmmState->getMixture()(m);
						GaussianInfo *gaussianInfo = &mGaussianInfo[gaussian];
						// do not split gaussians that were just merged
						if (gaussianInfo->bMerged) {
							continue;
						}
						// the occupation must be at least two times the minimum gaussian occupation
						if (gaussianInfo->dOccupation < 2.0*fMinimumGaussianOccupation) {
							continue;
						}	
						// find the heaviest mixture	
						if (gaussian->weightMixtureIncrement() > fWeightMax) {
							iComponentHeaviest = m;
							fWeightMax = gaussian->weightMixtureIncrement();
						}
					}
					if (iComponentHeaviest != -1) {	
						hmmState->getMixture()(iComponentHeaviest)->weightMixtureIncrement()--;
						iGaussianToSplit = iComponentHeaviest;
					}
					break;
				}
				default: {
					assert(0);
				}
			}	
			// finish in case there are no mixtures to split
			if (iGaussianToSplit == -1) {
				return iAdded;
			}
		}
		// case 2: just one Gaussian: it will be splitted, just check that the occupation count is big enough
		else {
			assert(hmmState->getMixture().getNumberComponents() == 1);
			GaussianInfo *gaussianInfo = &mGaussianInfo[hmmState->getMixture()(0)];
			// do not split Gaussian components that were just merged
			if (gaussianInfo->bMerged) {
				return iAdded;
			}
			// the occupation must be at least two times the minimum Gaussian occupation
			if (gaussianInfo->dOccupation < 2*fMinimumGaussianOccupation) {
				return iAdded;
			}
		}
		
		Gaussian *gaussianToSplit = hmmState->getMixture()(iGaussianToSplit);
		GaussianInfo *gaussianToSplitInfo = &mGaussianInfo[gaussianToSplit];
		GaussianInfo gaussianInfo;
		gaussianInfo.bMerged = false;
		gaussianInfo.dOccupation = 0.0;
		assert(gaussianToSplitInfo->dOccupation >= 2*fMinimumGaussianOccupation);
		
		// (2) split the mixture by creating two new mixtures, the new mixtures have a perturbed mean (+/-0.X standard deviations), identical variance and half the weight
		Gaussian *gaussian = new Gaussian(m_iDim,m_hmmManager->getCovarianceModelling());
		// mark both Gaussian components as merged so in case that multiple splits are 
		// done we do not split a Gaussian that results from a split
		if (bAllowMultipleSplitsPerComponent == false) {
			gaussianInfo.bMerged = true;
			gaussianToSplitInfo->bMerged = true;	
		}
		// compute the new mean and keep the covariance identical
		Vector<double> vPerturbation(gaussianToSplit->covarianceDiag());
		vPerturbation.sqrt();
		gaussian->mean().copy(gaussianToSplit->mean());
		gaussian->mean().add(1.0f*fSplittingEpsilon,vPerturbation);
		gaussianToSplit->mean().add(-1.0f*fSplittingEpsilon,vPerturbation);
		gaussian->covarianceDiag().copy(gaussianToSplit->covarianceDiag());
		// divide the weight and occupation by two	
		gaussianToSplit->weight() /= 2.0;	
		gaussianToSplitInfo->dOccupation /= 2.0;
		gaussian->weight() = gaussianToSplit->weight();
		gaussianInfo.dOccupation = gaussianToSplitInfo->dOccupation;
		hmmState->getMixture().addGaussianComponent(gaussian);	
		mGaussianInfo.insert(MGaussianInfo::value_type(gaussian,gaussianInfo));
		++iAdded; 
	}
	
	assert(hmmState->getMixture().getNumberComponents() == (int)mGaussianInfo.size());
	
	return iAdded;
}

// double the number of Gaussian components in the mixture
// there is no need to select which Gaussian components will be split, since all of them will (except those 
// that were just merged or dont have enough occupation)
int GMMEditor::mixtureDouble(HMMState *hmmState, float fMinimumGaussianOccupation, float fSplittingEpsilon) {

	assert(hmmState->getMixture().getNumberComponents() > 0);
	int iAdded = 0;
	
	MGaussianInfo &mGaussianInfo = m_mGaussianInfo[hmmState->getId()];
	
	// split all the Gaussian components except those that were just merged 
	int iGaussiansOriginal = hmmState->getMixture().getNumberComponents();
	for(int g = 0 ; g < iGaussiansOriginal ; ++g) {
	
		Gaussian *gaussian1 = hmmState->getMixture()(g);
		GaussianInfo *gaussianInfo1 = &mGaussianInfo[gaussian1];	
		
		// skip components with not enough occupation or those that were just merged
		if ((gaussianInfo1->dOccupation < 2.0f*fMinimumGaussianOccupation) ||
			(gaussianInfo1->bMerged)) {
			continue;	
		}
		
		// split the mixture by creating two new mixtures, the new mixtures have a perturbed mean (+/-0.X standard deviations), identical variance and half the weight
		Gaussian *gaussian2 = new Gaussian(m_iDim,m_hmmManager->getCovarianceModelling());
		GaussianInfo gaussianInfo2;
		gaussianInfo2.bMerged = false;
		gaussianInfo2.dOccupation = 0.0;		
		// compute the new mean and keep the covariance identical
		Vector<double> vPerturbation(gaussian1->covarianceDiag());
		vPerturbation.sqrt();
		gaussian2->mean().copy(gaussian1->mean());
		gaussian2->mean().add(1.0f*fSplittingEpsilon,vPerturbation);
		gaussian1->mean().add(-1.0f*fSplittingEpsilon,vPerturbation);
		gaussian2->covarianceDiag().copy(gaussian1->covarianceDiag());
		gaussian1->weight() /= 2.0f;
		gaussianInfo1->dOccupation /= 2.0f;
		gaussian2->weight() = gaussian1->weight();
		gaussianInfo2.dOccupation = gaussianInfo1->dOccupation;	
		hmmState->getMixture().addGaussianComponent(gaussian2);
		mGaussianInfo.insert(MGaussianInfo::value_type(gaussian2,gaussianInfo2));
		++iAdded;
	}	
	
	assert(hmmState->getMixture().getNumberComponents() == (int)mGaussianInfo.size());
	
	return iAdded;
}

};	// end-of-namespace

