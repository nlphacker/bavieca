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

#ifndef MLACCUMULATOR_H
#define MLACCUMULATOR_H

#include "ConfigurationFeatures.h"
#include "LexiconManager.h"
#include "PhoneSet.h"
#include "ForwardBackwardX.h"
#include "HMMState.h"
#include "HMMManager.h"
#include "MLFFile.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "Time.h"
#include "Accumulator.h"
#include "Log.h"
#include "PhoneticRulesManager.h"
/**
	@author daniel <dani.bolanos@gmail.com>
*/
class MLAccumulator {

	private:
	
		const char *m_strFilePhoneSet;
		const char *m_strFileConfigurationFeaturesAlignment; 
		const char *m_strFolderFeaturesAlignment;
		const char *m_strFileModelsAlignment;
		unsigned char m_iAccumulatorType;
		unsigned char m_iContextModelingOrderAccumulatorsWW;
		unsigned char m_iContextModelingOrderAccumulatorsCW;
		const char *m_strFileOptionalSymbols;
		bool m_bMultiplePronunciations;
		const char *m_strFileConfigurationFeaturesAcc;
		const char *m_strFolderFeaturesAcc;
		int m_iCovarianceModellingAcc;
		const char *m_strFileLexicon;
		const char *m_strFileMLF;
		const char *m_strFileAccumulators;
		const char *m_strFileGlobalDistribution;
		float m_fForwardPruningBeam;
		float m_fBackwardPruningBeam; 
		int m_iTrellisMaxSize;
		bool m_bTrellisCache;
		int m_iTrellisCacheMaxSize;
		const char *m_strFileLog;
		int m_iLogVerbosity;
		
		// feature stream: it can be either single or double (different features for alignment/estimation)
		bool m_bSingleFeatureStream;
		
		// feature dimensionality
		int m_iFeatureDimensionalityAlignment;
		int m_iFeatureDimensionalityAcc;
		
		// log 
		Log *m_log;
		
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		MLFFile *m_mlfFile;
		HMMManager *m_hmmManagerAlignment;
		HMMManager *m_hmmManagerAccumulation;
		ForwardBackwardX *m_forwardBackwardX;
		
		// optional lex units
		VLexUnit m_vLexUnitOptional;
		
		// accumulators
		MAccumulatorPhysical mAccumulatorPhysical;
		MAccumulatorLogical mAccumulatorLogical;
		
		// dump the global distribution to a file
		//bool dumpGlobalDistribution(const char *strFile);	

	public:

		// constructor
		MLAccumulator(const char *strFilePhoneSet, const char *strFileConfigurationFeaturesAlignment, 
			const char *strFolderFeaturesAlignment, const char *strFileModelsAlignment, 
			unsigned char iAccumulatorType, const char *strFileOptionalSymbols, bool bMultiplePronunciations,
			unsigned char iContextModelingOrderAccumulatorsWW, unsigned char iContextModelingOrderAccumulatorsCW,
			const char *strFileConfigurationFeaturesAcc, const char *strFolderFeaturesAcc, int iCovarianceModellingAcc, const char *strFileLexicon, const char *strFileMLF, const char *strFileAccumulators, const char *strFileGlobalDistribution, float fForwardPruningBeam, float fBackwardPruningBeam, int iTrellisMaxSize, bool bTrellisCache, int iTrellisCacheMaxSize, const char *strFileLog, int iLogVerbosity);

		// destructor
		~MLAccumulator();
		
		// initialize the accumulation
		bool initialize();
		
		// accumulate statistics
		bool accumulate();

};

#endif


