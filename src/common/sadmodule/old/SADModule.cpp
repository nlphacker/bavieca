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


#include "SADModule.h"

#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "PhoneSet.h"
#include "TimeUtils.h"

// constructor
SADModule::SADModule(PhoneSet *phoneSet, HMMManager *hmmManager, int iMaxGaussianComponentsSilence,
	int iMaxGaussianComponentsSpeech, float fPenaltySilenceToSpeech, int iFramesPadding) {
	
	m_phoneSet = phoneSet;
	m_hmmManager = hmmManager;
	m_iMaxGaussianComponentsSilence = iMaxGaussianComponentsSilence;
	m_iMaxGaussianComponentsSpeech = iMaxGaussianComponentsSpeech;
	m_hmmStateSilence = NULL;
	m_hmmStateSpeech = NULL;
	m_fPenaltySilenceToSpeech = fPenaltySilenceToSpeech;
	m_iFramesPadding = iFramesPadding;
	m_iDim = m_hmmManager->getFeatureDimensionality();
	
	assert(m_iMaxGaussianComponentsSpeech > 0);
}

// destructor
SADModule::~SADModule() {

	if (m_hmmStateSilence != NULL) {
		assert(m_hmmStateSpeech != NULL);
		delete m_hmmStateSilence;
		delete m_hmmStateSpeech;
	}
}

// initialize the SAD system
void SADModule::initialize() {
	
	// get the HMM-states
	int iHMMStates = -1;
	HMMStateDecoding *hmmStates = m_hmmManager->getHMMStatesDecoding( &iHMMStates);
	assert(hmmStates);
	
	// estimate the states needed 
	unsigned int iComponentsSilence = 0;
	unsigned int iComponentsSpeech = 0;
	unsigned char iPhoneSil = m_phoneSet->getPhoneIndexSilence();
	for(int i=0 ; i < iHMMStates ; ++i) {
		if (hmmStates[i].getPhone() == iPhoneSil) {
			iComponentsSilence += hmmStates[i].getMixtureSize();
		} else {
			iComponentsSpeech += min(hmmStates[i].getMixtureSize(),(unsigned int)m_iMaxGaussianComponentsSpeech);
		}
	}
	 
	// create the Gaussian mixtures
	GaussianDecoding *gaussiansSilence = new GaussianDecoding[iComponentsSilence];
	unsigned int iComponentsSilenceFound = 0;
	float fWeightAccumulatedSilence = 0.0;
	
	GaussianDecoding *gaussiansSpeech = new GaussianDecoding[iComponentsSpeech];
	unsigned int iComponentsSpeechFound = 0;
	float fWeightAccumulatedSpeech = 0.0;
	for(int i=0 ; i < iHMMStates ; ++i) {
	
		// get the Gaussian components
		int iGaussianComponents = -1;
		GaussianDecoding *gaussians = hmmStates[i].getGaussians(iGaussianComponents);
		LGaussianDecoding lGaussian;
		for(int j=0 ; j < iGaussianComponents ; ++j) {
			lGaussian.push_back(gaussians+j);
		}
		
		// sort them by weight
		lGaussian.sort(HMMStateDecoding::compareGaussianByWeight);	
	
		if (hmmStates[i].getPhone() == iPhoneSil) {	
		
			// keep the best and accumulate the weight
			unsigned int iCount = 0;
			for(LGaussianDecoding::iterator it = lGaussian.begin() ; (it != lGaussian.end()) && (iCount < m_iMaxGaussianComponentsSilence) ; ++it, ++iCount) {
				HMMStateDecoding::copyGaussian(&(gaussiansSilence[iComponentsSilenceFound]),*it,m_iDim);
				fWeightAccumulatedSilence += (*it)->fWeight;
				++iComponentsSilenceFound;
			}	
		
		} else {
		
			// keep the best and accumulate the weight
			unsigned int iCount = 0;
			for(LGaussianDecoding::iterator it = lGaussian.begin() ; (it != lGaussian.end()) && (iCount < m_iMaxGaussianComponentsSpeech) ; ++it, ++iCount) {
				HMMStateDecoding::copyGaussian(&(gaussiansSpeech[iComponentsSpeechFound]),*it,m_iDim);
				fWeightAccumulatedSpeech += (*it)->fWeight;
				++iComponentsSpeechFound;
			}	
		}
		lGaussian.clear();
	}
	
	// readjust weights	
	for(int i=0 ; i < iComponentsSilence ; ++i) {
		gaussiansSilence[i].fWeight /= fWeightAccumulatedSilence;
	}
	for(int i=0 ; i < iComponentsSpeech ; ++i) {
		gaussiansSpeech[i].fWeight /= fWeightAccumulatedSpeech;
	}
	
	// set Gaussian components
	m_hmmStateSilence = new HMMStateDecoding(m_iDim,m_phoneSet,
		UCHAR_MAX,UCHAR_MAX,UCHAR_MAX,iHMMStates,iComponentsSilenceFound,gaussiansSilence);
	m_hmmStateSpeech = new HMMStateDecoding(m_iDim,m_phoneSet,
		UCHAR_MAX,UCHAR_MAX,UCHAR_MAX,iHMMStates+1,iComponentsSpeech,gaussiansSpeech);	
	
	// initialize the HMM-states (precompute constants for the evaluation)
	m_hmmStateSilence->initialize();
	m_hmmStateSpeech->initialize();
}

// print the segment information to the standard output
void SADModule::printSegments(VSpeechSegment &vSpeechSegment) {

	printf("# segments: %d\n",vSpeechSegment.size());
	for(VSpeechSegment::iterator it = vSpeechSegment.begin() ; it != vSpeechSegment.end() ; ++it) {
		printf("%8d %8d\n",(*it)->iFrameStart,(*it)->iFrameEnd);
	}	
}

// store audio segments to disk
void SADModule::store(const char *strFile, VSpeechSegment &vSpeechSegment) {

	FileOutput file(strFile,false);
	file.open();
	
	for(VSpeechSegment::iterator it = vSpeechSegment.begin() ; it != vSpeechSegment.end() ; ++it) {
		ostringstream ossText;
		ossText << (*it)->iFrameStart << " " << (*it)->iFrameEnd << endl;
		cout << ossText.str();
		IOBase::writeString(file.getStream(),ossText);
	}
	
	file.close();
}

// load audio segments from disk
VSpeechSegment *SADModule::load(const char *strFile) {

	FileInput file(strFile,false);
	file.open();	
	
	string strLine;
	unsigned int iFrameStart = UINT_MAX;
	unsigned int iFrameEnd = UINT_MAX;
	VSpeechSegment *vSpeechSegment = new VSpeechSegment;
	while(std::getline(file.getStream(),strLine).good()) {
		if (strLine.empty()) {
			break;
		}
		std::stringstream s(strLine);
		IOBase::read(s,&iFrameStart);	
		IOBase::read(s,&iFrameEnd);
		if (iFrameStart <= iFrameEnd) {
			BVC_ERROR << "incosistent start and end of segment";	
		}
		vSpeechSegment->push_back(newSpeechSegment(iFrameStart,iFrameEnd));
	}	
	
	file.close();
	
	//printSegments(*vSpeechSegment);

	return vSpeechSegment;
}

// initializes a SAD session
void SADModule::beginSession() {

	m_iTimeFrame = 0; 
}

// process the given features
void SADModule::processFeatures(float *fFeatures, int iFeatures) {

	for(int i=0 ; i < iFeatures ; ++i) {
		m_grid.push_back(new GridElementSAD[2*HMM_STATES_CLASS]);	
	}
	
	int iSize = m_grid.size();

	// reset time-stamps to prevent getting scores from the cache
	m_hmmStateSilence->resetTimeStamp();
	m_hmmStateSpeech->resetTimeStamp();
	int iTimeFrameLast = m_iTimeFrame+iFeatures;
	
	// first time frame?
	if (m_iTimeFrame == 0) {
	
		// allocate memory for the grid
		for(unsigned int i=0 ; i < 2*HMM_STATES_CLASS ; ++i) {
			m_grid[0][i].fScore = -FLT_MAX ;
			m_grid[0][i].iPrev = -1;
		} 
				
		// initialize the first states
		m_grid[0][0].fScore = m_hmmStateSilence->computeEmissionProbabilityNearestNeighborPDE(fFeatures,0);
		m_grid[0][HMM_STATES_CLASS].fScore = m_hmmStateSpeech->computeEmissionProbabilityNearestNeighborPDE(fFeatures,0) + m_fPenaltySilenceToSpeech;
		m_iTimeFrame++;
	}
	
	// fill the grid
	int iPrevSelf = -1;
	int iPrevLeft = -1;
	float *fFeatureVector = fFeatures;
	for(unsigned int i=m_iTimeFrame ; i < iTimeFrameLast ; ++i,++m_iTimeFrame) {
		
		// silence states (all of them share the same mixture)		
		float fSilenceScore = -FLT_MAX;
		#ifdef SIMD
			fSilenceScore = m_hmmStateSilence->computeEmissionProbabilityNearestNeighborSIMD(fFeatureVector,i);	
		#else
			fSilenceScore = m_hmmStateSilence->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,i);
		#endif		
		for(unsigned int j=0 ; j < HMM_STATES_CLASS ; ++j) {
					
			float fLeft = -FLT_MAX;
			float fSelf = -FLT_MAX;
			// self transition
			if (j < i) {
				fSelf = m_grid[i-1][j].fScore;
				iPrevSelf = j;
			}
			// left to right silence to silence transition
			if ((j > 0) && (i > 0)) {
				fLeft = m_grid[i-1][j-1].fScore;
				iPrevLeft = j-1;
			} 
			// left to right speech to silence transition
			else if ((j == 0) && (i >= HMM_STATES_CLASS)) {
				fLeft = m_grid[i-1][2*HMM_STATES_CLASS-1].fScore;
				iPrevLeft = (2*HMM_STATES_CLASS)-1;	
			}
			// get the best backtrace at time t	
			if (fLeft > fSelf) {
				m_grid[i][j].fScore = fLeft+fSilenceScore;
				m_grid[i][j].iPrev = iPrevLeft;
			} else {
				m_grid[i][j].fScore = fSelf+fSilenceScore;
				m_grid[i][j].iPrev = iPrevSelf;
			}	
		}
		
		// speech states (all of them share the same mixture)
		float fSpeechScore = -FLT_MAX;
		#ifdef SIMD
			fSpeechScore = m_hmmStateSpeech->computeEmissionProbabilityNearestNeighborSIMD(fFeatureVector,i);	
		#else
			fSpeechScore = m_hmmStateSpeech->computeEmissionProbabilityNearestNeighborPDE(fFeatureVector,i);
		#endif		
		for(unsigned int j=0 ; j < HMM_STATES_CLASS ; ++j) {
					
			float fLeft = -FLT_MAX;
			float fSelf = -FLT_MAX;
			// self transition
			if (j < i) {
				fSelf = m_grid[i-1][j+HMM_STATES_CLASS].fScore;
				iPrevSelf = j+HMM_STATES_CLASS;
			}
			// left to right silence to silence transition
			if ((j > 0) && (i > 0)) {
				fLeft = m_grid[i-1][j+HMM_STATES_CLASS-1].fScore;
				iPrevLeft = j-1+HMM_STATES_CLASS;
			} 
			// left to right silence to speech transition
			else if ((j == 0) && (i >= HMM_STATES_CLASS)) {
				fLeft = m_grid[i-1][HMM_STATES_CLASS-1].fScore + m_fPenaltySilenceToSpeech;	
				iPrevLeft = HMM_STATES_CLASS-1;
			}
			// get the best backtrace at time t	
			if (fLeft > fSelf) {
				m_grid[i][j+HMM_STATES_CLASS].fScore = fLeft+fSpeechScore;
				m_grid[i][j+HMM_STATES_CLASS].iPrev = iPrevLeft;
			} else {
				m_grid[i][j+HMM_STATES_CLASS].fScore = fSelf+fSpeechScore;
				m_grid[i][j+HMM_STATES_CLASS].iPrev = iPrevSelf;
			}
		}
		fFeatureVector += m_iDim;
	}
}

// recover speech segments by doing back-tracking on the grid
void SADModule::recoverSpeechSegments(VSpeechSegment &vSpeechSegment) {

	// not enough features: return no speech segment
	if (m_iTimeFrame < HMM_STATES_CLASS) {		
		BVC_ERROR << "insufficient number of features, minimum number required: " << HMM_STATES_CLASS;
	}

	// recover the best sequence of states doing back-tracking
	GridElementSAD *elementAux = NULL;
	SpeechSegment *segment = new SpeechSegment;
	unsigned char iTypeCurrent = UCHAR_MAX;
	
	GridElementSAD *elementLastSilence = &m_grid.back()[HMM_STATES_CLASS-1];
	GridElementSAD *elementLastSpeech = &m_grid.back()[2*HMM_STATES_CLASS-1];
	// last state: silence
	if (elementLastSilence->fScore > elementLastSpeech->fScore) {
		iTypeCurrent = AUDIO_SEGMENT_SILENCE;
		elementAux = elementLastSilence;
	} 
	// last state: speech
	else {
		iTypeCurrent = AUDIO_SEGMENT_SPEECH;
		elementAux = elementLastSpeech;
	}	
	segment->iFrameEnd = m_iTimeFrame-1;
	
	unsigned int iFrame = m_iTimeFrame-2;
	while(elementAux->iPrev != -1) {
		unsigned char iTypeBack;
		(elementAux->iPrev >= HMM_STATES_CLASS) ? iTypeBack = AUDIO_SEGMENT_SPEECH : iTypeBack = AUDIO_SEGMENT_SILENCE;
		if (iTypeBack != iTypeCurrent) {
			segment->iFrameStart = iFrame+1;
			assert((segment->iFrameEnd-segment->iFrameStart) >= (HMM_STATES_CLASS-1));
			if (iTypeCurrent == AUDIO_SEGMENT_SPEECH) {
				vSpeechSegment.push_back(segment);
			} else {
				delete segment;
			}
			segment = newSpeechSegment(-1,iFrame);
			iTypeCurrent = iTypeBack;
		}
		elementAux = &m_grid[iFrame-1][elementAux->iPrev];
		--iFrame;
	} 
	segment->iFrameStart = 0;
	assert((segment->iFrameEnd-segment->iFrameStart) >= (HMM_STATES_CLASS-1));
	if (iTypeCurrent == AUDIO_SEGMENT_SPEECH) {
		vSpeechSegment.push_back(segment);
	} else {
		delete segment;
	}
	
	// reverse the elements in the vector
	for(int i=0 ; i < vSpeechSegment.size()/2 ; ++i) {
		SpeechSegment *segmentAux = vSpeechSegment[i];
		vSpeechSegment[i] = vSpeechSegment[vSpeechSegment.size()-(i+1)];
		vSpeechSegment[vSpeechSegment.size()-(i+1)] = segmentAux;
	}
	
	//printSegments(*vSpeechSegment);
	
	// (4) apply padding to the left and right
	unsigned int iSegment = 0;
	SpeechSegment *segmentPrev = NULL;
	for(VSpeechSegment::iterator it = vSpeechSegment.begin() ; it != vSpeechSegment.end() ; ++it, ++iSegment) {
		// left padding (IMP: note that the previous segment was already right padded)
		if (iSegment == 0) {
			if ((*it)->iFrameStart-m_iFramesPadding > 0) {
				(*it)->iFrameStart -= m_iFramesPadding;
			} else {
				(*it)->iFrameStart = 0;
			}
		} else {
			if ((*it)->iFrameStart-m_iFramesPadding > segmentPrev->iFrameEnd) {
				(*it)->iFrameStart -= m_iFramesPadding;
			} else {
				(*it)->iFrameStart = segmentPrev->iFrameEnd+1;
			}
			assert((*it)->iFrameStart > segmentPrev->iFrameEnd);	
		}
		// right padding
		if (iSegment == vSpeechSegment.size()-1) {
			if ((*it)->iFrameEnd+m_iFramesPadding <= m_iTimeFrame-1) {
				(*it)->iFrameEnd += m_iFramesPadding;
			} else {
				(*it)->iFrameEnd = m_iTimeFrame-1;
			}
		} else {
			// get the next segment
			VSpeechSegment::iterator jt = it;
			jt++;
			// apply padding
			unsigned int iFramesSilence = (*jt)->iFrameStart-(*it)->iFrameEnd-1;
			assert(iFramesSilence>=HMM_STATES_CLASS);
			if ((*it)->iFrameEnd+m_iFramesPadding < ((*jt)->iFrameStart-(iFramesSilence/2)+1)) {
				(*it)->iFrameEnd += m_iFramesPadding;
			} else {
				(*it)->iFrameEnd = (*it)->iFrameEnd+(iFramesSilence/2)+1;
			}
			assert((*jt)->iFrameStart > (*it)->iFrameEnd);	
		}
		segmentPrev = *it;
	}
}

// terminates a SAD session
void SADModule::endSession() {

	for(VGridElementSAD::iterator it = m_grid.begin() ; it != m_grid.end() ; ++it) {
		delete [] *it;
	}
}

// print the grid using for dynamic programming
void SADModule::printGrid() {

	for(unsigned int j=0 ; j < m_iTimeFrame ; ++j) {	
		printf("%12d (xx)",j);
	}
	printf("\n");
	for(unsigned int i=0 ; i < 2*HMM_STATES_CLASS ; ++i) {	
		for(unsigned int j=0 ; j < m_iTimeFrame ; ++j) {	
			if (m_grid[j][i].fScore == -FLT_MAX) {
				printf("%12.2f (xx)",0.0);
			} else {
				printf("%12.2f (%2d)",m_grid[j][i].fScore,m_grid[j][i].iPrev);
			}
		}
		printf("\n");
	}
}



