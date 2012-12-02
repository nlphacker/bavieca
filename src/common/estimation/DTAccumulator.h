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


#ifndef DTACCUMULATOR_H
#define DTACCUMULATOR_H

#include "LexiconManager.h"
#include "ForwardBackwardX.h"
#include "ForwardBackward.h"
#include "HMMManager.h"

namespace Bavieca {

class Accumulator;
class ConfigurationFeatures;
class FileUtils;
class TimeUtils;
class PhoneSet;
class MLFFile;
class FeatureFile;

// discriminative training objective functions
#define DISCRIMINATIVE_TRAINING_OBJECTIVE_FUNCTION_MMI			"MMI"				// Maximum Mutual Information
#define DISCRIMINATIVE_TRAINING_OBJECTIVE_FUNCTION_BMMI			"bMMI"			// Boosted Maximum Mutual Information

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class DTAccumulator {

	private: 
	
		const char *m_strFilePhoneSet;
		const char *m_strFileConfigurationFeatures; 
		const char *m_strFolderFeatures;
		const char *m_strFileModels;
		const char *m_strFileOptionalSymbols;
		bool m_bMultiplePronunciations;
		const char *m_strFileLexicon;
		const char *m_strFileMLF;
		const char *m_strFolderLattices;
		const char *m_strFileAccumulatorsNum;
		const char *m_strFileAccumulatorsDen;
		float m_fForwardPruningBeam;
		float m_fBackwardPruningBeam; 
		int m_iTrellisMaxSize;
		bool m_bTrellisCache;
		int m_iTrellisCacheMaxSize;
		
		// boosted MMI
		bool m_bMMI;
		float m_fBoostingFactor;
		
		int m_iCovarianceModeling;
		int m_iGaussianComponents;
		int m_iHMMStates;
		
		// discriminative trainin parameters
		const char *m_strObjectiveFunction;
		bool m_bCanceledStatistics;
		
		// feature dimensionality
		int m_iFeatureDimensionality;
		
		// scaling factors
		float m_fScaleAM;
		float m_fScaleLM;
		
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		MLFFile *m_mlfFile;
		HMMManager *m_hmmManager;
		ForwardBackwardX *m_forwardBackwardX;			// used for numerator statistics
		ForwardBackward *m_forwardBackward;				// used for denominator statistics
		
		// optional lex units
		VLexUnit m_vLexUnitOptional;
		
		// accumulators
		MAccumulatorPhysical m_mAccumulatorNum;
		MAccumulatorPhysical m_mAccumulatorDen;
		
		// statistics cancellation (between numerator and denominator)
		void statisticsCancellation(Alignment *alignmentNum, MOccupation *mOccupationDen);
		
		// accumulate statistics
		void accumulate(Alignment *alignment, float *fFeatures, int iFeatures, bool bNumerator);

	public:

		// constructor
		DTAccumulator(const char *strFilePhoneSet, const char *strFileConfigurationFeatures, 
			const char *strFolderFeatures, const char *strFileModels, const char *strFileOptionalSymbols, 
			bool bMultiplePronunciations, const char *strFileLexicon, const char *strFileMLF, 
			const char *strFolderLattices, float fScaleAM, float fScaleLM, const char *strFileAccumulatorsNum, 
			const char *strFileAccumulatorsDen, const char *strObjectiveFunction, float fBoostingFactor,
			bool bCanceledStatistics, float fForwardPruningBeam, float fBackwardPruningBeam, 
			int iTrellisMaxSize, bool bTrellisCache, int iTrellisCacheMaxSize);

		// destructor
		~DTAccumulator();
		
		// initialize the accumulation
		void initialize();
		
		// accumulate statistics
		void accumulate();

};

};	// end-of-namespace

#endif
