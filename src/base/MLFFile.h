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

#ifndef MLFFILE_H
#define MLFFILE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;

#include <string>
#include <vector>

#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "LogMessage.h"
#include "LexiconManager.h"
#include "Global.h"
#include "Log.h"

#define MODE_READ			0
#define MODE_WRITE		1

typedef struct {
	string strFilePattern;		// file pattern representing an utterance
	VLexUnit vLexUnit;			// transcription: sequence of lexical units (pronunciation sensitive)
} MLFUtterance;

typedef vector<MLFUtterance*> VMLFUtterance;

/**
	@author root <root@localhost.localdomain>
*/
class MLFFile {

	private:
	
		VMLFUtterance m_vMLFUtterance1;						// utterances in the MLF
		VMLFUtterance m_vMLFUtterance2;						// utterances in the MLF
		VMLFUtterance m_vMLFUtterance;						// utterances in the MLF
		LexiconManager *m_lexiconManager;					// lexicon manager
		string m_strFile;											// file name
		Log *m_log;													// log object	
		unsigned char m_iMode;									// access mode
		char m_strMessage[MAX_LOG_MESSAGE_LENGTH];
		
		// return whether the given string is a feature filename
		bool isFeatureFilename(string &strFile);	
		
		// extract the lexical unit identity and the pronunciation number
		LexUnit *parseLexUnitLine(string &strLexUnit, int iMLFLine);
		
		// log an error
		void logError(const char *strError);
	
	public:
	
		// constructor
		MLFFile(LexiconManager *lexiconManager, const char *strFile, Log *log, const unsigned char iMode);

		// destructor
		~MLFFile();
		
		// load the MLF file from disk
		void load();	
		
		// store the MLF utterances to disk
		void store();
		
		// return the utterances
		inline VMLFUtterance *getUtterances() {
		
			return &m_vMLFUtterance;
		}
		
		// add an utterance to the MLF (store mode)
		inline bool addUtterance(MLFUtterance *mlfUtterance) {
		
			if (m_iMode != MODE_WRITE) {
				return false;
			}
		
			m_vMLFUtterance.push_back(mlfUtterance);
			
			return true;
		}
		

};

#endif
