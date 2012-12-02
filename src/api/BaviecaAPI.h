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


#ifndef BAVIECAAPI_H
#define BAVIECAAPI_H

using namespace std;

#include <string>

#include <stdlib.h>
#include <assert.h>

namespace Bavieca {

// under windows it is necessary to export the public interface
#ifdef _WIN32
//#define DLLEXPORT extern "C" __declspec(dllexport)
//#define DLLEXPORT __declspec(dllexport)
//#define DLLEXPORT __declspec(dllimport)
#define DLLEXPORT
#elif __linux__
#define DLLEXPORT
#endif

class BestPath;
class ConfigurationBavieca;
class DynamicDecoderX;
class DynamicNetworkX;
class FeatureExtractor;
class FeatureTransformList;
class FileUtils;
class FillerManager;
class HMMManager;
class LexiconManager;
class LMManager;
class NetworkBuilderX;
class PhoneSet;
class SADModule;
class ViterbiX;

// initialization modes
enum {INIT_SAD=1, INIT_ALIGNER=2, INIT_DECODER=4, INIT_ADAPTATION=8};

// text alignment events
enum {TAE_CORRECT=1, TAE_SUBSTITUTION=2, TAE_DELETION=4, TAE_INSERTION=8};

// keep original alignment and align structure members on 8-byte boundaries
#pragma pack(push,8)

// speech segment
typedef struct {
	int iFrameStart;			// starting time frame (zero-based hundredths of a second)
	int iFrameEnd;				// ending time frame	(zero-based hundredths of a second)
} SpeechSegmentI;

// phone alignment
typedef struct {
	char *strPhone;
	int iFrameStart;
	int iFrameEnd;
} PhoneAlignmentI;

// word alignmnet
typedef struct {
	char *strWord;				// word
	int iFrameStart;
	int iFrameEnd;	
	int iPhones;				// number of phonemes in the word
	PhoneAlignmentI *phoneAlignment;			// phone alignments
} WordAlignmentI;

// configuration parameters (for overriding)
typedef struct {
	const char *strParameter;
	const char *strValue;
} ParamValueI;

// word hypothesis
typedef struct {
	char *strWord;				// word
	int iFrameStart;
	int iFrameEnd;
} WordHypothesisI;

// text alignment
typedef struct {
	unsigned char iEvent;
	int iIndexRef;
	char *strWordRef;
	int iIndexHyp;
	char *strWordHyp;
} TAElementI;

// restore original alignment
#pragma pack(pop)

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class DLLEXPORT BaviecaAPI {

	private:
	
		unsigned char m_iFlags;					// initialization mode
		string m_strFileConfiguration;
		bool m_bInitialized;
		
		ConfigurationBavieca *m_configuration;
		PhoneSet *m_phoneSet;
		LexiconManager *m_lexiconManager;
		SADModule *m_sadModule;
		HMMManager *m_hmmManager;
		FeatureExtractor *m_featureExtractor;
		LMManager *m_lmManager;
		ViterbiX *m_viterbiX;
		DynamicNetworkX *m_network;
		NetworkBuilderX *m_networkBuilder;
		DynamicDecoderX *m_dynamicDecoder;
		bool m_bLatticeGeneration;

	public:

		// constructor (receives default configuration parameters)
		DLLEXPORT BaviecaAPI(const char *strFileConfiguration);

		// destructor
		DLLEXPORT ~BaviecaAPI();
		
		// API INITIALIZATION ------------------------------------------------------------------------------------
		
		// initialize API (overriding parameters as needed)
		DLLEXPORT bool initialize(unsigned char iFlags, ParamValueI *paramValue = NULL, int iParameters = 0);
		
		// uninitialize the API	
		DLLEXPORT void uninitialize();
		
		// FEATURE EXTRACTION -------------------------------------------------------------------------------------
		
		// extract features from the audio
		DLLEXPORT float *extractFeatures(short *sSamples, int iSamples, int *iFeatures);
		
		// return feature dimensionality
		DLLEXPORT int getFeatureDim();
		
		// free features extracted using extractFeatures(...)
		DLLEXPORT void free(float *fFeatures);
		
		// SPEECH ACTIVITY DETECTION ------------------------------------------------------------------------------
		
		// start a SAD session
		DLLEXPORT void sadBeginSession();
		
		// terminate a SAD session
		DLLEXPORT void sadEndSession();
		
		// proces the given features
		DLLEXPORT void sadFeed(float *fFeatures, int iFeatures);
		
		// recover speech segments by doing back-tracking on the grid
		DLLEXPORT SpeechSegmentI *sadRecoverSpeechSegments(int *iSegments);	
		
		// free speech segments returned by sarRecoverSpeechSegments(...)
		DLLEXPORT void free(SpeechSegmentI *speechSegments);
		
		// FORCED ALIGNMENT --------------------------------------------------------------------------------------
		
		// forced alignment between features and audio
		DLLEXPORT WordAlignmentI *align(float *fFeatures, int iFeatures, const char *strText, 
			bool bMultiplePronunciations, int *iWords);
			
		// free word alignments returned by align(...)
		DLLEXPORT void free(WordAlignmentI *wordAlignments, int iWords);
			
		// DECODING ----------------------------------------------------------------------------------------------
		
		// signals beginning of utterance
		DLLEXPORT void decBeginUtterance();
		
		// process feature vectors from an utterance
		DLLEXPORT void decProcess(float *fFeatures, int iFeatures);
		
		// get decoding results
		DLLEXPORT WordHypothesisI *decGetHypothesis(int *iWords, const char *strFileHypothesisLattice = NULL);
		
		// signals end of utterance
		DLLEXPORT void decEndUtterance();
		
		// free word hypotheses returned by getHypothesis(...)
		DLLEXPORT void free(WordHypothesisI *wordHypothesis, int iWords);
		
		// return a word-level assessment given a hypothesis and a reference text
		DLLEXPORT TAElementI *getAssessment(WordHypothesisI *wordHypothesis, int iWordsHyp, const char *strReference,
			int *iTAElements);
			
		// free text alignment elements returned by getAssessment(...)
		DLLEXPORT void free(TAElementI *taElements, int iElements);	
				
		// SPEAKER ADAPTATION ------------------------------------------------------------------------------------
		
		// feed data into speaker adaptation
		DLLEXPORT void mllrFeed(const char *strReference, float *fFeatures, int iFeatures);
		
		// adapt current acoustic models using fed adaptation data
		DLLEXPORT void mllrAdapt();

};

};	// end-of-namespace

#endif
