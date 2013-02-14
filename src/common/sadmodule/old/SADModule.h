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


#ifndef SADMODULE_H
#define SADMODULE_H

using namespace std;

#include <vector>

#include "HMMManager.h"

class PhoneSet;

#define HMM_STATES_CLASS		5	// number of HMM-states for Silence/Speech

typedef struct {
	float fScore;
	int iPrev;
} GridElementSAD;

typedef vector<GridElementSAD*> VGridElementSAD;

// types of audio segments
#define AUDIO_SEGMENT_SILENCE			0
#define AUDIO_SEGMENT_SPEECH			1

typedef struct {
	int iFrameStart;			// starting time frame
	int iFrameEnd;				// ending time frame
} SpeechSegment;

typedef vector<SpeechSegment*> VSpeechSegment;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class SADModule {

	private:
	
		PhoneSet *m_phoneSet;						// phonetic symbol set
		HMMManager *m_hmmManager;					// HMMs
		int m_iMaxGaussianComponentsSilence;	// maximum number of Gaussian components
		int m_iMaxGaussianComponentsSpeech;		// maximum number of Gaussian components
		float m_fPenaltySilenceToSpeech;			// penalty applied when moving from a silence segment to a speech segment
		HMMStateDecoding *m_hmmStateSilence;	// silence HMM-state
		HMMStateDecoding *m_hmmStateSpeech;		// speech HMM-state
		int m_iFramesPadding;						// # of frames used to pad the beginning and end of the speech segment
		VGridElementSAD m_grid;						// grid to do the Viterbi search
		int m_iTimeFrame;								// current time frame within the session
		int m_iDim;										// feature dimensionality
		
		static SpeechSegment *newSpeechSegment(int iFrameStart, int iFrameEnd) {
			
			SpeechSegment *segment = new SpeechSegment;
			segment->iFrameStart = iFrameStart;
			segment->iFrameEnd = iFrameEnd;
			return segment;
		}
		
		// print the grid using for dynamic programming
		void printGrid();	

	public:

		// constructor
		SADModule(PhoneSet *phoneSet, HMMManager *hmmManager, int iMaxGaussianComponentsSilence,
			int iMaxGaussianComponentsSpeech, float fPenaltySilenceToSpeech, int iFramesPadding);
		
		// destructor
		~SADModule();
		
		// initialize the SAD system
		void initialize();
		
		// initializes a SAD session
		void beginSession();
		
		// terminates a SAD session
		void endSession();
		
		// proces the given features
		void processFeatures(float *fFeatures, int iFeatures);
		
		// recover speech segments by doing back-tracking on the grid
		void recoverSpeechSegments(VSpeechSegment &vSpeechSegment);	
		
		// print the segment information to the standard output
		static void printSegments(VSpeechSegment &vSpeechSegment);

		// store audio segments to disk
		static void store(const char *strFile, VSpeechSegment &vSpeechSegment);

		// load audio segments from disk
		static VSpeechSegment *load(const char *strFile);
		
		// deallocates memory for a vector of speech segments
		static void deleteVSpeechSegment(VSpeechSegment &vSpeechSegment) {
			
			for(VSpeechSegment::iterator it = vSpeechSegment.begin() ; it != vSpeechSegment.end() ; ++it) {
				delete *it;
			}
		}
};

#endif
