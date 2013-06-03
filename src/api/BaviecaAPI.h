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
#include <vector>

#include <stdlib.h>
#include <assert.h>

namespace Bavieca {

class BestPath;
class ConfigurationBavieca;
class DynamicDecoderX;
class DynamicNetworkX;
class FeatureExtractor;
class FeatureTransformList;
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

// helper classes
class AlignmentI;
class TextAlignmentI;
class ParamValuesI;
class PhoneAlignmentI;
class WordAlignmentI;
class AlignmentI;
class HypothesisI;
class WordHypothesisI;
class SpeechSegmentsI;


/**
	@author daniel <dani.bolanos@gmail.com>
*/

#ifdef _MSC_VER
class extern "C" __declspec(dllexport) BaviecaAPI {
#elif defined __linux__ || defined __APPLE__ || __MINGW32__ || defined SWIG
class BaviecaAPI {
#endif

	private:
	
		unsigned char m_iFlags;					// initialization mode
		char *m_strFileConfiguration;
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
		BaviecaAPI(const char *strFileConfiguration);

		// destructor
		~BaviecaAPI();
		
		// API INITIALIZATION ------------------------------------------------------------------------------------
		
		// initialize API (overriding parameters as needed)
		bool initialize(unsigned char iFlags, ParamValuesI *paramValues = NULL);
		
		// uninitialize the API	
		void uninitialize();
		
		// FEATURE EXTRACTION -------------------------------------------------------------------------------------
		
		// extract features from the audio
		float *extractFeatures(short *sSamples, unsigned int iSamples, unsigned int *iFeatures);
		
		// return feature dimensionality
		int getFeatureDim();
		
		// free features extracted using extractFeatures(...)
		void free(float *fFeatures);
		
		// SPEECH ACTIVITY DETECTION ------------------------------------------------------------------------------
		
		// start a SAD session
		void sadBeginSession();
		
		// terminate a SAD session
		void sadEndSession();
		
		// proces the given features
		void sadFeed(float *fFeatures, unsigned int iFeatures);
		
		// recover speech segments by doing back-tracking on the grid
		SpeechSegmentsI *sadRecoverSpeechSegments();	
		
		// FORCED ALIGNMENT --------------------------------------------------------------------------------------
		
		// forced alignment between features and audio
		AlignmentI *align(float *fFeatures, unsigned int iFeatures, const char *strText, bool bMultiplePronunciations);
			
		// DECODING ----------------------------------------------------------------------------------------------
		
		// signals beginning of utterance
		void decBeginUtterance();
		
		// process feature vectors from an utterance
		void decProcess(float *fFeatures, unsigned int iFeatures);
		
		// get decoding results
		HypothesisI *decGetHypothesis(const char *strFileHypothesisLattice = NULL);
		
		// signals end of utterance
		void decEndUtterance();
		
		// return a word-level assessment given a hypothesis and a reference text
		TextAlignmentI *getAssessment(HypothesisI *hypothesis, const char *strReference);
				
		// SPEAKER ADAPTATION ------------------------------------------------------------------------------------
		
		// feed data into speaker adaptation
		void mllrFeed(const char *strReference, float *fFeatures, unsigned int iFeatures);
		
		// adapt current acoustic models using fed adaptation data
		void mllrAdapt();

};

// helper-classes

class SpeechSegmentI {

	private:

		int m_iFrameStart;     // starting time frame (zero-based hundredths of a second)
		int m_iFrameEnd;       // ending time frame (zero-based hundredths of a second)

	public:

		SpeechSegmentI(int iFrameStart, int iFrameEnd) {

			m_iFrameStart = iFrameStart;
			m_iFrameEnd = iFrameEnd;
		}

		int getFrameStart() {
			return m_iFrameStart;
		}
		int getFrameEnd() {
			return m_iFrameEnd;
		}
};

class SpeechSegmentsI {

	friend class BaviecaAPI;

	private: 

		vector<SpeechSegmentI*> m_vSpeechSegmentsI;

		SpeechSegmentsI(vector<SpeechSegmentI*> &vSpeechSegmentsI) {			
			m_vSpeechSegmentsI.assign(vSpeechSegmentsI.begin(),vSpeechSegmentsI.end());
		}

	public:

		~SpeechSegmentsI() {
			for(vector<SpeechSegmentI*>::iterator it = m_vSpeechSegmentsI.begin() ; it != m_vSpeechSegmentsI.end() ; ++it) {
				delete *it;
			}
			m_vSpeechSegmentsI.clear();
		}

		unsigned int size() {
			return (unsigned int)m_vSpeechSegmentsI.size();
		}
	
		SpeechSegmentI *getSpeechSegment(unsigned int iSegment) {
			if (iSegment < size()) {
				return m_vSpeechSegmentsI[iSegment];
			}
			return NULL;
		}
};

// configuration parameters (for overriding)
class ParamValueI {

	private:
	
		string m_strParameter;
		string m_strValue;

	public:

		ParamValueI(const char *strParameter, const char *strValue) {
				
			assert((strParameter) && (strValue)); 
			m_strParameter = strParameter;
			m_strValue = strValue;
		}

		const char *getParameter() {
			return m_strParameter.c_str();
		}
		const char *getValue() {
			return m_strValue.c_str();
		}
};

class ParamValuesI {

	friend class BaviecaAPI;

	private: 

		vector<ParamValueI*> m_vParamValuesI;

		ParamValuesI(vector<ParamValueI*> &vParamValuesI) {			
			m_vParamValuesI.assign(vParamValuesI.begin(),vParamValuesI.end());
		}

	public:

		~ParamValuesI() {
			for(vector<ParamValueI*>::iterator it = m_vParamValuesI.begin() ; it != m_vParamValuesI.end() ; ++it) {
				delete *it;
			}
			m_vParamValuesI.clear();
		}

		unsigned int size() {
			return (unsigned int)m_vParamValuesI.size();
		}

		ParamValueI *getParamValue(unsigned int i) {
			if (i < size()) {
				return m_vParamValuesI[i];
			}
			return NULL;
		}
};

class WordHypothesisI {

	private:

		string m_strWord;
		int m_iFrameStart;
		int m_iFrameEnd;

	public:

		WordHypothesisI(const char *strWord, int iFrameStart, int iFrameEnd) {

			m_strWord = (strWord) ? strWord : "";
			m_iFrameStart = iFrameStart;
			m_iFrameEnd = iFrameEnd;
		}

		~WordHypothesisI() {
		}

		const char *getWord() {
			return m_strWord.c_str();
		}
		int getFrameStart() {
			return m_iFrameStart;
		}
		int getFrameEnd() {
			return m_iFrameEnd;
		}
};

class HypothesisI {

	friend class BaviecaAPI;

	private: 

		vector<WordHypothesisI*> m_vWordHypothesisI;

		HypothesisI(vector<WordHypothesisI*> &vWordHypothesisI) {			
			m_vWordHypothesisI.assign(vWordHypothesisI.begin(),vWordHypothesisI.end());
		}

	public:

		~HypothesisI() {
			for(vector<WordHypothesisI*>::iterator it = m_vWordHypothesisI.begin() ; it != m_vWordHypothesisI.end() ; ++it) {
				delete *it;
			}
			m_vWordHypothesisI.clear();
		}

		unsigned int size() {
			return (unsigned int)m_vWordHypothesisI.size();
		}
	
		WordHypothesisI *getWordHypothesis(unsigned int iWord) {
			if (iWord < size()) {
				return m_vWordHypothesisI[iWord];
			}
			return NULL;
		}
};

// phone alignment
class PhoneAlignmentI {

	private:

		string m_strPhone;
		int m_iFrameStart;
		int m_iFrameEnd;

	public:

		PhoneAlignmentI(const char *strPhone, int iFrameStart, int iFrameEnd) {
			m_strPhone = strPhone ? strPhone : "";
			m_iFrameStart = iFrameStart;
			m_iFrameEnd = iFrameEnd;
		}

		const char *getPhone() {
			return m_strPhone.c_str();
		}
		int getFrameStart() {
			return m_iFrameStart;
		}
		int getFrameEnd() {
			return m_iFrameEnd;
		}
};

// word alignmnet
class WordAlignmentI {

	friend class AlignmentI;
	friend class BaviecaAPI;

	private:

		string m_strWord;								// word
		int m_iFrameStart;								// start time frame
		int m_iFrameEnd;								// end time frame
		vector<PhoneAlignmentI*> m_vPhoneAlignmentI;	// phone alignments

		WordAlignmentI(const char *strWord, int iFrameStart, int iFrameEnd, vector<PhoneAlignmentI*> &vPhoneAlignmentI) {
			m_strWord = strWord ? strWord : "";
			m_iFrameStart = iFrameStart;
			m_iFrameEnd = iFrameEnd;
			m_vPhoneAlignmentI.assign(vPhoneAlignmentI.begin(),vPhoneAlignmentI.end());
		}

	public:

		~WordAlignmentI() {
			for(vector<PhoneAlignmentI*>::iterator it = m_vPhoneAlignmentI.begin() ; it != m_vPhoneAlignmentI.end() ; ++it) {
				delete *it;
			}
			m_vPhoneAlignmentI.clear();
		}

		const char *getWord() {
			return m_strWord.c_str();
		}
		int getFrameStart() {
			return m_iFrameStart;
		}
		int getFrameEnd() {
			return m_iFrameEnd;
		}
		unsigned int size() {
			return (unsigned int)m_vPhoneAlignmentI.size();
		}

		PhoneAlignmentI *getPhoneAlignment(unsigned int iPhone) {
			if (iPhone < size()) {
				return m_vPhoneAlignmentI[iPhone];
			} 
			return NULL;			
		}
};

class AlignmentI {

	friend class BaviecaAPI;

	private:

		vector<WordAlignmentI*> m_vWordAlignmentI;

		AlignmentI(vector<WordAlignmentI*> &vWordAlignmentI) {	
			m_vWordAlignmentI.assign(vWordAlignmentI.begin(),vWordAlignmentI.end());
		}

	public:

		~AlignmentI() {
			for(vector<WordAlignmentI*>::iterator it = m_vWordAlignmentI.begin() ; it != m_vWordAlignmentI.end() ; ++it) {
				delete *it;
			}
			m_vWordAlignmentI.clear();
		}

		unsigned int size() {
			return (unsigned int)m_vWordAlignmentI.size();
		}
	
		WordAlignmentI *getWordAlignment(unsigned int iWord) {

			if (iWord < size()) {
				return m_vWordAlignmentI[iWord];
			}
			return NULL;
		}
};

// text alignment element
class TextAlignmentElementI {

	private:

		unsigned char m_iEvent;
		int m_iIndexRef;
		string m_strWordRef;
		int m_iIndexHyp;
		string m_strWordHyp;

	public:

		TextAlignmentElementI(unsigned char iEvent, int iIndexRef, const char *strWordRef, int iIndexHyp, const char *strWordHyp) {
			m_iEvent = iEvent;
			m_iIndexRef = iIndexRef;
			m_strWordRef = (strWordRef) ? strWordRef : "";
			m_iIndexHyp = iIndexHyp;
			m_strWordHyp = (strWordHyp) ? strWordHyp : "";
		}

		unsigned char getEvent() {
			return m_iEvent;
		}
		unsigned int getIndexRef() {
			return m_iIndexRef;
		}
		const char *getWordRef() {
			return m_strWordRef.c_str();
		}
		unsigned int getIndexHyp() {
			return m_iIndexHyp;
		}
		const char *getWordHyp() {
			return m_strWordHyp.c_str();
		}
};

// text alignment
class TextAlignmentI {

	friend class BaviecaAPI;

	private:

		std::vector<TextAlignmentElementI*> m_vElements;

		TextAlignmentI(vector<TextAlignmentElementI*> &vElements) {
			m_vElements.assign(vElements.begin(),vElements.end());
		}

	public:

		~TextAlignmentI() {
			for(std::vector<TextAlignmentElementI*>::iterator it = m_vElements.begin() ; it != m_vElements.end() ; ++it) {
				delete *it;
			}
			m_vElements.clear();
		}

		unsigned int size() {
			return (unsigned int)m_vElements.size();
		}

		TextAlignmentElementI *getElement(unsigned int iElement) {
			if (iElement < size()) {
				return m_vElements[iElement];
			}
			return NULL;
		}
};

};	// end-of-namespace

#endif

