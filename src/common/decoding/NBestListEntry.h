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
 

#ifndef NBESTLISTENTRY_H
#define NBESTLISTENTRY_H

using namespace std;

#include <vector>

#include "NBestListEntryElement.h"

namespace Bavieca {

class LexiconManager;

typedef vector<NBestListEntryElement*> VNBestListEntryElement;

/**
	@author root <dani.bolanos@gmail.com>
*/
class NBestListEntry {

	private:
	
		LexiconManager *m_lexiconManager;		// lexicon manager
		double m_dLikelihood;						// path likelihood
		double m_dPP;									// path posterior probability
		VNBestListEntryElement m_vElements;		// elements in the entry (typically words)

	public:

		// constructor
		NBestListEntry(LexiconManager *lexiconManager);

		// destructor
		~NBestListEntry();
		
		// store to disk
		void store(ostream &os, bool bTextFormat = false);	
		
		// add an element to the entry
		void add(NBestListEntryElement *element) {
		
			m_vElements.push_back(element);
			m_dLikelihood += element->getLikelihoodAM() + element->getLikelihoodLM() + element->getIP();
			m_dPP += element->getPP();
		}
			
		// return the likelihood
		double getLikelihood() {
		
			return m_dLikelihood;
		}
		
		// set the likelihood
		void setLikelihood(double dLikelihood) {
		
			m_dLikelihood = dLikelihood;
		}
		
		// print the entry
		void print();

};

};	// end-of-namespace

#endif
