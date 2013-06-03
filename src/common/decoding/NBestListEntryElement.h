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


#ifndef NBESTLISTENTRYELEMENT_H
#define NBESTLISTENTRYELEMENT_H

#include "LexiconManager.h"

namespace Bavieca {

class LexiconManager;

/**
	@author root <dani.bolanos@gmail.com>
*/
class NBestListEntryElement {

	private:
	
		int m_iBegin;
		int m_iEnd;
		LexUnit *m_lexUnit;
		double m_dLikelihoodAM;
		double m_dLikelihoodLM;
		double m_dIP;
		double m_dPP;	

	public:

		// constructor
		NBestListEntryElement(int iBegin, int iEnd, LexUnit *lexUnit, 
			double dAM, double dLM, double dIP, double dPP);

		// destructor
		~NBestListEntryElement();
		
		// store to disk
		void store(ostream &os, bool bTextFormat = false, LexiconManager *lexiconManager = NULL);	
		
		// return the lexical unit
		LexUnit *getLexUnit() {
			return m_lexUnit;
		}
		
		// return the acoustic model likelihood
		double getLikelihoodAM() {
			return m_dLikelihoodAM;
		}
		
		// return the language model likelihood
		double getLikelihoodLM() {
			return m_dLikelihoodLM;
		}
		
		// return the insertion penalty applied
		double getIP() {
			return m_dIP;
		}
		
		// return the posterior probability
		double getPP() {
			return m_dPP;
		}
};

};	// end-of-namespace

#endif
