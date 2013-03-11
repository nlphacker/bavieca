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


#include "AudioFile.h"
#include "AlignmentFile.h"
#include "BatchFile.h"
#include "ConfigurationFeatures.h"
#include "FeatureExtractor.h"
#include "FeatureFile.h"
#include "HMMManager.h"
#include "LexiconManager.h"
#include "Viterbi.h"
#include "ViterbiX.h"
#include "VTLEstimator.h"

#include <iostream>
#include <iomanip>

namespace Bavieca {

// constructor
VTLEstimator::VTLEstimator(ConfigurationFeatures *configurationFeatures, HMMManager *hmmManager, 
	PhoneSet *phoneSet, LexiconManager *lexiconManager, const char *strFileFillerPhones, 
	float fWarpFactorFloor, float fWarpFactorCeiling, float fWarpFactorStep, bool bRealignData,
	int iCepstralNormalizationMode, int iCepstralNormalizationMethod) {

	m_configurationFeatures = configurationFeatures;
	m_hmmManager = hmmManager;
	m_phoneSet = phoneSet;
	m_lexiconManager = lexiconManager;
	if (strFileFillerPhones == NULL) {
		m_strFileFillerPhones = "";
	} else {
		m_strFileFillerPhones = strFileFillerPhones;
	}
	m_fWarpFactorFloor = fWarpFactorFloor;
	m_fWarpFactorCeiling = fWarpFactorCeiling;
	m_fWarpFactorStep = fWarpFactorStep;
	m_bRealignData = bRealignData;
	m_iCepstralNormalizationMode = iCepstralNormalizationMode;
	m_iCepstralNormalizationMethod = iCepstralNormalizationMethod;	
	m_iDim = m_hmmManager->getFeatureDimensionality();
}

// destructor
VTLEstimator::~VTLEstimator() {

}


// estimate the warp factor from a batch file containing alignment files
void VTLEstimator::estimateWarpFactor(const char *strBatchFile, unsigned char iAlignmentFormat,
	const char *strOutputFile, float &fLikelihoodGain, float &fWarpFactor) {

	// filler phones	
	bool *bFillerPhones = new bool[m_phoneSet->size()];
	for(unsigned int i=0; i < m_phoneSet->size() ; ++i) {
		bFillerPhones[i] = false;
	}
	
	// (1) load the list of filler phones if any	
	if (m_strFileFillerPhones.compare("") != 0) {
		BatchFile fileFillerPhones(m_strFileFillerPhones.c_str(),"phone");
		fileFillerPhones.load();
		for(unsigned int i=0 ; i < fileFillerPhones.size() ; ++i) {
			const char *strPhone = fileFillerPhones.getField(i,0u);
			int iIndex = m_phoneSet->getPhoneIndex(strPhone);
			if (iIndex == -1) {
				BVC_ERROR << "filler phonetic symbol \"" << strPhone << "\" was not found in the phonetic symbol set";
			}
			bFillerPhones[iIndex] = true;		
		}
	} 
	// default: just silence
	else {
		bFillerPhones[m_phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE)] = true;	
	}

	// (2) load the batch file containing data to estimate the warp factors
	BatchFile batchFile(strBatchFile,"raw|alignment");
	batchFile.load();
	
	// (3) count the number of different warp factors to test
	int iTests = 0;
	for(float f = m_fWarpFactorFloor ; f <= m_fWarpFactorCeiling ; f += m_fWarpFactorStep) {
		++iTests;
	}
	if (iTests < 1) {
		BVC_ERROR << "no estimation needed given warp factors range and step size";
	}
	double *dLikelihoodTest = new double[iTests];	
	for(int i=0 ; i<iTests ; ++i) {
		dLikelihoodTest[i] = 0.0;
	}
	
	// (4) extract session-data
	VUtteranceData vUtteranceData;
	for(unsigned int iUtterance = 0 ; iUtterance < batchFile.size() ; ++iUtterance) {
	
		// load the raw audio
		int iSamples = -1;
		short int *sSamples = AudioFile::load(batchFile.getField(iUtterance,"raw"),&iSamples);
	
		UtteranceData utteranceData;
		utteranceData.samples.sSamples = sSamples;
		utteranceData.samples.iSamples = iSamples;
		utteranceData.features.fFeatures = NULL;
		utteranceData.features.iFeatures = -1;
		vUtteranceData.push_back(utteranceData);
	}	
	
	float m_fBeamWidth = 1000;
	ViterbiX viterbiX(m_phoneSet,m_lexiconManager,m_hmmManager,m_fBeamWidth,1000,true,1000);
	
	VLexUnit vLexUnitOptional;
	vLexUnitOptional.push_back(m_lexiconManager->getLexUnitSilence());
	
	float *fLikelihoodUtteranceBest = new float[vUtteranceData.size()]; 
	float *fWarpFactorUtteranceBest = new float[vUtteranceData.size()];
	for(unsigned int i=0 ; i < vUtteranceData.size() ; ++i) {
		fLikelihoodUtteranceBest[i] = -FLT_MAX;
		fWarpFactorUtteranceBest[i] = -1.0;
	}
	
	m_hmmManager->resetHMMEmissionProbabilityComputation();
	
	// test different warp factors (from lower to higher)
	int iTest = 0;
	for(float f = m_fWarpFactorFloor ; f <= m_fWarpFactorCeiling ; f += m_fWarpFactorStep, ++iTest) {	
			
		// extract warped features
		FeatureExtractor featureExtractor(m_configurationFeatures,f,-1,
			m_iCepstralNormalizationMode,m_iCepstralNormalizationMethod);
		featureExtractor.initialize();
		featureExtractor.extractFeaturesSession(vUtteranceData);
		
		// for each utterance
		for(unsigned int iUtterance = 0 ; iUtterance < batchFile.size() ; ++iUtterance) {	
			
			float fLikelihood;
			float fLikelihoodUtterance = 0.0;
			
		
			// load the alignment information
			VPhoneAlignment *vPhoneAlignmentBase = NULL;
			VLexUnit vLexUnit;
			// text format
			if (iAlignmentFormat == FILE_FORMAT_TEXT) {
				
				AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,m_lexiconManager);
				vPhoneAlignmentBase = alignmentFile->load(batchFile.getField(iUtterance,"alignment"));
				if (vPhoneAlignmentBase == NULL) {
					BVC_ERROR << "unable to load the alignment file";
				}
				delete alignmentFile;
				
				if (m_bRealignData) {
					VLexUnitAlignment *vLexUnitAlignment = AlignmentFile::getVLexUnitAlignment(*vPhoneAlignmentBase);
					assert(vLexUnitAlignment);
					for(VLexUnitAlignment::iterator it = vLexUnitAlignment->begin() ; it != vLexUnitAlignment->end() ; ++it) {
						if ((*it)->lexUnit != m_lexiconManager->getLexUnitSilence()) {
							vLexUnit.push_back((*it)->lexUnit);
						}
					}
					AlignmentFile::destroyLexUnitAlignment(vLexUnitAlignment);
				}
			}
			// binary format
			else {
				assert(iAlignmentFormat == FILE_FORMAT_BINARY);
				
				Alignment *alignment = Alignment::load(batchFile.getField(iUtterance,"alignment"),m_lexiconManager);
				if (alignment == NULL) {
					BVC_ERROR << "unable to load the alignment file";
				}
				vPhoneAlignmentBase = alignment->getPhoneAlignment(m_lexiconManager);
				assert(vPhoneAlignmentBase);
				if (m_bRealignData) {
					VWordAlignment *vWordAlignment = alignment->getWordAlignment();
					for(VWordAlignment::iterator it = vWordAlignment->begin() ; it != vWordAlignment->end() ; ++it) {
						if ((*it)->iLexUnitPron != m_lexiconManager->getLexUnitSilence()->iLexUnitPron) {
							vLexUnit.push_back(m_lexiconManager->getLexUnitPron((*it)->iLexUnitPron));
						}
					}
				}
				delete alignment;	
			}
				
			// align warped features to the given word sequence
			double dUttLikelihood;
			int iErrorCode = -1;
			float *fFeatures = vUtteranceData[iUtterance].features.fFeatures;
			int iFeatures = vUtteranceData[iUtterance].features.iFeatures;
			
			VPhoneAlignment *vPhoneAlignment = vPhoneAlignmentBase; 
			if (m_bRealignData) {
				Alignment *alignment = viterbiX.processUtterance(vLexUnit,false,vLexUnitOptional,
					fFeatures,iFeatures,&dUttLikelihood,iErrorCode);
				if (alignment == NULL) {
					BVC_ERROR << "unable to perform the alignment";
				}
				vPhoneAlignment = alignment->getPhoneAlignment(m_lexiconManager);
				delete alignment;
			}
			
			// check consistency
			if (iFeatures != vPhoneAlignment->back()->iStateEnd[NUMBER_HMM_STATES-1]+1) {
				BVC_ERROR << "timing inconsitency detected between alignment file: " << batchFile.getField(iUtterance,"alignment") << " and raw file: " << 
				batchFile.getField(iUtterance,"raw");
			}	
				
			// compute the utterance likelihood for the given warp factor
			unsigned char iPhoneLeft = m_phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE);
			unsigned char iPhoneRight;
			VHMMStateDecoding vHMMStateDecoding;
			for(VPhoneAlignment::iterator it = vPhoneAlignment->begin() ; it != vPhoneAlignment->end() ; ++it) {
				// determine the right phonetic context
				VPhoneAlignment::iterator jt = it;
				advance(jt,1);
				if (jt == vPhoneAlignment->end()) {
					iPhoneRight = m_phoneSet->getPhoneIndex(PHONETIC_SYMBOL_SILENCE);
				} else {
					iPhoneRight = (*jt)->iPhone;
				}
				if (bFillerPhones[(*it)->iPhone] == false) {
					for(unsigned char iState = 0 ; iState < NUMBER_HMM_STATES; ++iState) {
						HMMStateDecoding *hmmStateDecoding = m_hmmManager->getHMMStateDecoding(&iPhoneLeft,(*it)->iPhone,&iPhoneRight,(*it)->iPosition,iState);
						vHMMStateDecoding.push_back(hmmStateDecoding);
						if (hmmStateDecoding == NULL) {
							BVC_ERROR << "HMM-state was not found, unable to compute acoustic likelihood for warp factor";
						}
						for(int iFrame = (*it)->iStateBegin[iState] ; iFrame <= (*it)->iStateEnd[iState] ; ++iFrame) {
							fLikelihood = hmmStateDecoding->computeEmissionProbabilityNearestNeighborPDE(
								fFeatures+(iFrame*m_iDim),iFrame);
							fLikelihoodUtterance += fLikelihood;
							dLikelihoodTest[iTest] += fLikelihood;
						}
					}	
				} 
				iPhoneLeft = (*it)->iPhone;
			}
				
			m_hmmManager->resetHMMEmissionProbabilityComputation(vHMMStateDecoding);
			vHMMStateDecoding.clear();
				
			
			// keep the best partial warp factor for the given utterance
			if (fLikelihoodUtterance > fLikelihoodUtteranceBest[iUtterance]) {
				fWarpFactorUtteranceBest[iUtterance] = f;
				fLikelihoodUtteranceBest[iUtterance] = fLikelihoodUtterance;
			} 
			
			// clean-up
			delete [] vUtteranceData[iUtterance].features.fFeatures;
			vUtteranceData[iUtterance].features.iFeatures = -1;
			if (m_bRealignData) {
				AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
			}
			AlignmentFile::destroyPhoneAlignment(vPhoneAlignmentBase);	
		}
	}	
	
	// write output data
	assert(strOutputFile);
	FileOutput fileOutput(strOutputFile,false);
	fileOutput.open();
	
	// get the best warp factor for each utterance
	ostringstream oss;
	for(unsigned int i=0 ; i < vUtteranceData.size() ; ++i) {
		oss << "utterance: " << std::setw(4) << i << " best warp factor: " << std::setw(6) << 
		std::setiosflags(ios::fixed) << std::setprecision(2) << fWarpFactorUtteranceBest[i] << 
		" likelihood: " << std::setw(14) << std::setiosflags(ios::fixed) << std::setprecision(4) <<
		fLikelihoodUtteranceBest[i] << endl;
	}
	
	// get the best warp factor globally
	float fWarpFactorBest = -1.0;
	double dLikelihoodBest = -DBL_MAX;
	float fLikelihoodBaseline = -FLT_MAX;
	iTest = 0;
	for(float f = m_fWarpFactorFloor ; f <= m_fWarpFactorCeiling ; f += m_fWarpFactorStep, ++iTest) {	
		if (dLikelihoodTest[iTest] > dLikelihoodBest) {
			dLikelihoodBest = dLikelihoodTest[iTest];
			fWarpFactorBest = f;
		}
		if (fabs(f-1.0) < 0.001) {
			fLikelihoodBaseline = dLikelihoodTest[iTest];
		}
	}

	// keep the best warp factor
	fWarpFactor = fWarpFactorBest;	
	
	fLikelihoodGain = 1.0-(dLikelihoodBest/fLikelihoodBaseline);
	oss << "global warp factor: " << std::setw(6) << std::setiosflags(ios::fixed) << std::setprecision(2) << 
		fWarpFactorBest << ", likelihood gain: " << std::setw(6) << std::setiosflags(ios::fixed) << 
		std::setprecision(2) << fLikelihoodGain*100 << "%" << endl;
	IOBase::writeString(fileOutput.getStream(),oss);
	fileOutput.close();
	
	// clean-up
	for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it) {
		delete [] (*it).samples.sSamples;
	}	
	
	delete [] fLikelihoodUtteranceBest; 
	delete [] fWarpFactorUtteranceBest;	
	
		// clean-up
	delete [] dLikelihoodTest;
	delete [] bFillerPhones;
}

};	// end-of-namespace


