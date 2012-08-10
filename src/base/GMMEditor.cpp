/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
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

// constructor
GMMEditor::GMMEditor(HMMManager *hmmManager, MAccumulatorPhysical *mAccumulatorPhysical)
{
	m_hmmManager = hmmManager;
	m_mAccumulatorPhysical = mAccumulatorPhysical;
	m_iFeatureDimensionality = hmmManager->getFeatureDimensionality();
	m_iCovarianceModeling = hmmManager->getCovarianceModelling();
	m_iCovarianceElements = HMMManager::getCovarianceElements(m_iFeatureDimensionality,m_iCovarianceModeling);
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
		m_mGaussianInfo[iHMMState].insert(MGaussianInfo::value_type(m_hmmStates[iHMMState]->getGaussian(iGaussian),gaussianInfo));
	}
		
	return;
}

// merge Gaussian components in the mixture
bool GMMEditor::mixtureMerge(float fMinimumGaussianOccupation, float fMinimumGaussianWeight, float fCovarianceFlooringRatio) {

	// compute the covariance floor
	double *dCovarianceFloor = MLEstimator::computeCovarianceFloor(m_hmmManager,*m_mAccumulatorPhysical,
		fCovarianceFlooringRatio);
	
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		if (mixtureMerge(m_hmmStates[i],fMinimumGaussianOccupation,fMinimumGaussianWeight,dCovarianceFloor) == false) {
			return false;	
		}
	}
	
	delete [] dCovarianceFloor;
	
	return true;
}

// merge Gaussian components with not enough training data associated
// when a Gaussian is split into two gaussians, and the parameters of the newly created 
// Gaussian components are reestimated (iteration), it is possible that one of the gaussians 
// has very low occupation, in that case we need to merge that Gaussian to the closest one
bool GMMEditor::mixtureMerge(HMMState *hmmState, float fMinimumGaussianOccupation, float fMinimumGaussianWeight, double *dCovarianceFloor) {

	assert(hmmState->getGaussianComponents() >= 1);
	
	MGaussianInfo &mGaussianInfo = m_mGaussianInfo[hmmState->getId()];
	
	// check that there are multiple Gaussian components
	if (hmmState->getGaussianComponents() == 1) {
		mGaussianInfo[hmmState->getGaussian(0)].bMerged = false;	
		return true;
	}
	
	// allocate memory for the accumulators
	float *fMean = new float[m_iFeatureDimensionality];
	float *fCovariance = new float[m_iFeatureDimensionality];
	
	// (1) find Gaussian components with occupation/weight below the threshold and merge them to the nearest one
	bool bMerged = false;
	do {	
		bMerged = false;	
		for(int g = 0 ; g < hmmState->getGaussianComponents() ; ++g) {
			Gaussian *gaussian1 = hmmState->getGaussian(g);
			GaussianInfo *gaussianInfo1 = &mGaussianInfo[gaussian1];
			if ((gaussianInfo1->dOccupation < fMinimumGaussianOccupation) || 
				(gaussian1->fWeight < fMinimumGaussianWeight)) {
				int iGaussianClosest = -1;
				// case 1: a single Gaussian (it can happen after all the Gaussian components are merged to one)
				if (hmmState->getGaussianComponents() == 1) {
					// warning: the HMM-state has just one Gaussian with little data, it will be poorly estimated	
					break;
				}
				// case 2: two Gaussian components, there is no need to compute the distance, just merge them 
				else if (hmmState->getGaussianComponents() == 2) {
					iGaussianClosest = (g+1)%2;	
				}	
				// case 3: more than two Gaussian -> find the closest gaussian to the poorly trained one
				else {	
					// compute the determinant of the covariance matrix of the gaussian to merge
					double dCovarianceDet1 = 0.0;
					for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {	
						dCovarianceDet1 += log(gaussian1->fCovariance[i]);
					}
					// find the closest gaussian	
					double dDistanceSmallest = DBL_MAX;
					for(int h = 0 ; h < hmmState->getGaussianComponents() ; ++h) {
						if (h == g) {
							continue;
						}
						Gaussian *gaussian2 = hmmState->getGaussian(h);
						GaussianInfo *gaussianInfo2 = &mGaussianInfo[gaussian2];
						// compute the distance between the gaussians
						// determinant of the second gaussian
						double dCovarianceDet2 = 0.0;
						for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {	
							dCovarianceDet2 += log(gaussian2->fCovariance[i]);
						}
						// compute the covariance matrix determinant of the hypothetical merged gaussian
						double dCovarianceDetMerged = 0.0;
						// compute the counts
						double dDenominator = gaussianInfo1->dOccupation+gaussianInfo2->dOccupation;
						double d1 = gaussianInfo1->dOccupation/dDenominator;
						double d2 = gaussianInfo2->dOccupation/dDenominator;
						// compute the determinant of the covariance matrix resulting from merging both gaussians
						for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {	
							dCovarianceDetMerged += log((d1*gaussian1->fCovariance[i])+(d2*gaussian2->fCovariance[i])+(d1*d2*(   (gaussian1->fMean[i]-gaussian2->fMean[i])*(gaussian1->fMean[i]-gaussian2->fMean[i]))));
						}
						// compute the distance
						double dDistance = (((gaussianInfo1->dOccupation+gaussianInfo2->dOccupation)/2.0)*dCovarianceDetMerged) - ((gaussianInfo1->dOccupation/2.0)*dCovarianceDet1) - ((gaussianInfo2->dOccupation/2.0)*dCovarianceDet2);
						assert(finite(dDistance) != 0);
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
				Gaussian *gaussian2 = hmmState->getGaussian(h);
				GaussianInfo *gaussianInfo2 = &mGaussianInfo[gaussian2];
				// compute the counts
				double dDenominator = gaussianInfo1->dOccupation+gaussianInfo2->dOccupation;
				double d1 = gaussianInfo1->dOccupation/dDenominator;
				double d2 = gaussianInfo2->dOccupation/dDenominator;	
				// compute the new mean and covariance
				for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {	
					fMean[i] = (float)((d1*gaussian1->fMean[i])+(d2*gaussian2->fMean[i]));
					fCovariance[i] = (float)((d1*gaussian1->fCovariance[i])+(d2*gaussian2->fCovariance[i])+(d1*d2*(   (gaussian1->fMean[i]-gaussian2->fMean[i])*(gaussian1->fMean[i]-gaussian2->fMean[i]))));
					// flooring
					if (fCovariance[i] < dCovarianceFloor[i]) {
						fCovariance[i] = dCovarianceFloor[i];
					}
				}
				// copy the mean and covariance
				for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {	
					gaussian1->fMean[i] = fMean[i];
					gaussian1->fCovariance[i] = fCovariance[i];
				}
				// compute the new weight
				gaussian1->fWeight += gaussian2->fWeight;
				// accumulate occupancy (this is an iterative bottom-up clustering, the occupancy might still be below the threshold)
				gaussianInfo1->dOccupation += gaussianInfo2->dOccupation;
				// mark the new gaussian as "merged" so it wont be split next
				gaussianInfo1->bMerged = true;
				// remove the gaussian component
				hmmState->removeGaussianComponent(h);
				mGaussianInfo.erase(mGaussianInfo.find(gaussian2));
				bMerged = true;
				break;
			}
		}
	} while(bMerged == true);
	
	assert(mGaussianInfo.size() == hmmState->getGaussianComponents());
	
	// deallocate memory
	delete [] fMean;
	delete [] fCovariance;	
	
	// sanity check: make sure all the gaussians have an occupancy and weight above the threshold
	if (hmmState->getGaussianComponents() > 1) {
		for(int g = 0 ; g < hmmState->getGaussianComponents() ; ++g) {
			Gaussian *gaussian = hmmState->getGaussian(g);
			GaussianInfo *gaussianInfo = &mGaussianInfo[gaussian];	
			assert(gaussian->fWeight >= fMinimumGaussianWeight);
			assert(gaussianInfo->dOccupation >= fMinimumGaussianOccupation);
		}
	}

	return true;
}

// add Gaussian components to each mixture
int GMMEditor::mixtureIncrement(int iIncrement, int iSplittingCriterion, 
	bool bAllowMultipleSplitsPerComponent, float fMinimumGaussianOccupation, float fSplittingEpsilon) {
	
	int iOriginal = 0;
	int iAdded = 0;
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		iOriginal += m_hmmStates[i]->getGaussianComponents();
		iAdded += mixtureIncrement(m_hmmStates[i],iIncrement,iSplittingCriterion,bAllowMultipleSplitsPerComponent,
			fMinimumGaussianOccupation,fSplittingEpsilon);
	}
	
	if (iAdded > 0) { 
		m_hmmManager->setSingleGaussian(false);
	}
	
	//printf("original: %d\n",iOriginal);
	//printf("added: %d\n",iAdded);
	
	return iAdded;
}

// double the number of Gaussian components in the mixture (if possible)
int GMMEditor::mixtureDouble(float fMinimumGaussianOccupation, float fSplittingEpsilon) {
	
	int iOriginal = 0;
	int iAdded = 0;
	for(int i=0 ; i < m_iHMMStates ; ++i) {
		iOriginal += m_hmmStates[i]->getGaussianComponents();
		iAdded += mixtureDouble(m_hmmStates[i],fMinimumGaussianOccupation,fSplittingEpsilon);
	}
	
	if (iAdded > 0) { 
		m_hmmManager->setSingleGaussian(false);
	}
	
	//printf("original: %d\n",iOriginal);
	//printf("added: %d\n",iAdded);
	
	return iAdded;
}

// add a Gaussian to the mixture that models the emission probability of the state
// note: there is no need to floor the covariance of the new components, since they stay the same
int GMMEditor::mixtureIncrement(HMMState *hmmState, int iIncrement, int iSplittingCriterion, 
	bool bAllowMultipleSplitsPerComponent, float fMinimumGaussianOccupation, float fSplittingEpsilon) {

	assert(hmmState->getGaussianComponents() > 0);
	assert(iIncrement > 0);
		
	MGaussianInfo &mGaussianInfo = m_mGaussianInfo[hmmState->getId()];
	int iAdded = 0;			// number of components already added to the mixture
	
	// initialize splitting weights if needed
	if (iSplittingCriterion == SPLITTING_CRITERION_HEAVIEST_MIXTURE_COMPONENT) {
		for(int m = 0 ; m < hmmState->getGaussianComponents() ; ++m) {
			Gaussian *gaussian = hmmState->getGaussian(m);
			// the weight for mixture increment determines which Gaussian will be split next
			// a Gaussian that has already been split gets its weight substracted by one, what makes it
			// very unlikely to be split again
			gaussian->fWeightMixtureIncrement = gaussian->fWeight;
		}	
	}
	
	while(iAdded < iIncrement) {
	
		int iGaussianToSplit = 0;
		// case 1: there are multiple gaussians: select the mixture to split
		if (hmmState->getGaussianComponents() > 1) {
			iGaussianToSplit = -1;
			switch(iSplittingCriterion) {
				case SPLITTING_CRITERION_LARGEST_AVERAGE_VARIANCE: {
					double dGeometricMean = 0.0;
					double dGeometricMeanHighest = -DBL_MAX;
					// for each gaussian compute the geometric mean of the variance
					for(int m = 0 ; m < hmmState->getGaussianComponents() ; ++m) {
						Gaussian *gaussian = hmmState->getGaussian(m);
						GaussianInfo *gaussianInfo = &mGaussianInfo[gaussian];
						// do not split gaussians that were just merged
						if (gaussianInfo->bMerged == true) {
							continue;
						}
						// the occupation must be at least two times the minimum gaussian occupation
						if (gaussianInfo->dOccupation < 2.0*fMinimumGaussianOccupation) {
							continue;
						}
						// compute the geometric mean of the covariance (log-average)
						dGeometricMean = 0.0;
						for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
							dGeometricMean += log(gaussian->fCovariance[i]);
						}
						dGeometricMean /= ((float)m_iFeatureDimensionality);
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
					for(int m = 0 ; m < hmmState->getGaussianComponents() ; ++m) {
						Gaussian *gaussian = hmmState->getGaussian(m);
						GaussianInfo *gaussianInfo = &mGaussianInfo[gaussian];
						// do not split gaussians that were just merged
						if (gaussianInfo->bMerged == true) {
							continue;
						}
						// the occupation must be at least two times the minimum gaussian occupation
						if (gaussianInfo->dOccupation < 2.0*fMinimumGaussianOccupation) {
							continue;
						}	
						// find the heaviest mixture	
						if (gaussian->fWeightMixtureIncrement > fWeightMax) {
							iComponentHeaviest = m;
							fWeightMax = gaussian->fWeightMixtureIncrement;
						}
					}
					if (iComponentHeaviest != -1) {	
						hmmState->getGaussian(iComponentHeaviest)->fWeightMixtureIncrement--;
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
			assert(hmmState->getGaussianComponents() == 1);
			GaussianInfo *gaussianInfo = &mGaussianInfo[hmmState->getGaussian(0)];
			// do not split Gaussian components that were just merged
			if (gaussianInfo->bMerged == true) {
				return iAdded;
			}
			// the occupation must be at least two times the minimum Gaussian occupation
			if (gaussianInfo->dOccupation < 2*fMinimumGaussianOccupation) {
				return iAdded;
			}
		}
		
		Gaussian *gaussianToSplit = hmmState->getGaussian(iGaussianToSplit);
		GaussianInfo *gaussianToSplitInfo = &mGaussianInfo[gaussianToSplit];
		GaussianInfo gaussianInfo;
		gaussianInfo.bMerged = false;
		gaussianInfo.dOccupation = 0.0;
		assert(gaussianToSplitInfo->dOccupation >= 2*fMinimumGaussianOccupation);
		
		// (2) split the mixture by creating two new mixtures, the new mixtures have a perturbed mean (+/-0.X standard deviations), identical variance and half the weight
		Gaussian *gaussian = HMMState::newGaussian(m_iFeatureDimensionality,m_iCovarianceElements);
		// mark both gaussians as merged so in case that more than one split is done we do not split a Gaussian that results from a split
		if (bAllowMultipleSplitsPerComponent == false) {
			gaussianInfo.bMerged = true;
			gaussianToSplitInfo->bMerged = true;	
		}
		// compute the new mean and keep the covariance identical
		for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {
			float fPerturbation = fSplittingEpsilon*sqrt(gaussianToSplit->fCovariance[i]);
			gaussian->fMean[i] = gaussianToSplit->fMean[i] + fPerturbation;
			gaussianToSplit->fMean[i] -= fPerturbation;
			gaussian->fCovariance[i] = gaussianToSplit->fCovariance[i];
		}
		// divide the weight and occupation by two	
		gaussianToSplit->fWeight /= 2.0;	
		gaussianToSplitInfo->dOccupation /= 2.0;
		gaussian->fWeight = gaussianToSplit->fWeight;
		gaussianInfo.dOccupation = gaussianToSplitInfo->dOccupation;
		hmmState->addGaussianComponent(gaussian);	
		mGaussianInfo.insert(MGaussianInfo::value_type(gaussian,gaussianInfo));
		++iAdded; 
	}
	
	assert(hmmState->getGaussianComponents() == mGaussianInfo.size());
	
	return iAdded;
}

// double the number of Gaussian components in the mixture
// there is no need to select which Gaussian components will be split, since all of them will (except those 
// that were just merged or dont have enough occupation)
int GMMEditor::mixtureDouble(HMMState *hmmState, float fMinimumGaussianOccupation, float fSplittingEpsilon) {

	assert(hmmState->getGaussianComponents() > 0);
	int iAdded = 0;
	
	MGaussianInfo &mGaussianInfo = m_mGaussianInfo[hmmState->getId()];
	
	// split all the Gaussian components except those that were just merged 
	int iGaussiansOriginal = hmmState->getGaussianComponents();
	for(int g = 0 ; g < iGaussiansOriginal ; ++g) {
	
		Gaussian *gaussian1 = hmmState->getGaussian(g);
		GaussianInfo *gaussianInfo1 = &mGaussianInfo[gaussian1];	
		
		// skip components with not enough occupation or those that were just merged
		if ((gaussianInfo1->dOccupation < 2.0f*fMinimumGaussianOccupation) ||
			(gaussianInfo1->bMerged == true)) {
			continue;	
		}
		
		// split the mixture by creating two new mixtures, the new mixtures have a perturbed mean (+/-0.X standard deviations), identical variance and half the weight
		Gaussian *gaussian2 = HMMState::newGaussian(m_iFeatureDimensionality,m_iCovarianceElements);
		GaussianInfo gaussianInfo2;
		gaussianInfo2.bMerged = false;
		gaussianInfo2.dOccupation = 0.0;		
		// compute the new mean and keep the covariance identical
		for(int i = 0 ; i < m_iFeatureDimensionality ; ++i) {	
			double dPerturbation = fSplittingEpsilon*sqrt(gaussian1->fCovariance[i]);
			gaussian2->fMean[i] = gaussian1->fMean[i] + dPerturbation;
			gaussian1->fMean[i] -= dPerturbation;
			gaussian2->fCovariance[i] = gaussian1->fCovariance[i];
		}	
		gaussian1->fWeight /= 2.0f;
		gaussianInfo1->dOccupation /= 2.0f;
		gaussian2->fWeight = gaussian1->fWeight;
		gaussianInfo2.dOccupation = gaussianInfo1->dOccupation;	
		hmmState->addGaussianComponent(gaussian2);
		mGaussianInfo.insert(MGaussianInfo::value_type(gaussian2,gaussianInfo2));
		++iAdded;
	}	
	
	assert(hmmState->getGaussianComponents() == mGaussianInfo.size());
	
	return iAdded;
}
