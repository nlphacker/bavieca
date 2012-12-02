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


#ifndef VTLESTIMATOR_H
#define VTLESTIMATOR_H

using namespace std;

#include <string>
#include <vector>

namespace Bavieca {

class ConfigurationFeatures;
class HMMManager;
class LexiconManager;
class PhoneSet;

/**
	@author root <root@localhost.localdomain>
*/
class VTLEstimator {

	private:
	
		ConfigurationFeatures *m_configurationFeatures;
		HMMManager *m_hmmManager;
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		string m_strFileFillerPhones;
		float m_fWarpFactorFloor;
		float m_fWarpFactorCeiling;
		float m_fWarpFactorStep;
		bool m_bRealignData;
		int m_iCepstralNormalizationMode;
		int m_iCepstralNormalizationMethod;
		int m_iDim;
		
		// return whether a string is blank
		bool isBlank(string &string) {
		
			for(unsigned int i=0 ; i<string.length() ; ++i) {
				if (isspace(string.c_str()[i]) == 0) {
					return false;
				}
			}
			
			return true;
		}		

	public:

		// constructor
		// in order to perform VTLN estimation the following items are needed:
		// - state alignment information (to determine the HMM to evaluate each feature frame)
		// - feature frames
		// - set of HMMs: the lower resolution the better, single gaussian triphone models are the best option
		// - set of lexical units that will be ignored during the estimation (non-speech lexical units). By default
		//   it will be silence plus the filler symbols defined in the lexicon
		
		VTLEstimator(ConfigurationFeatures *configurationFeatures, HMMManager *hmmManager, PhoneSet *phoneSet, 
			LexiconManager *lexiconManager, const char *strFileFillerPhones, float fWarpFactorFloor, float fWarpFactorCeiling, float fWarpFactorStep, bool bRealignData, int iCepstralNormalizationMode, int iCepstralNormalizationMethod);	

		// destructor
    	~VTLEstimator();
    	
    	// initialization
    	bool initialize();
    	
    	// estimate the warping factor from a batch file containing alignment files
    	void estimateWarpFactor(const char *strBatchFile, unsigned char iAlignmentFormat, 
    		const char *strOutputFile, float &fLikelihoodGain, float &fWarpFactor);

};

};	// end-of-namespace

#endif
