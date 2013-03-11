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


#ifndef PHONETICRULESMANAGER_H
#define PHONETICRULESMANAGER_H

#define PHONETIC_RULE_NO		0
#define PHONETIC_RULE_YES		1

#define MAX_RULE_NAME_LENGTH	128

#include <string.h>

using namespace std;

#include <string>
#include <vector>

#include "Global.h"

namespace Bavieca {

class FileInput;
class FileOutput;
class IOBase;
class PhoneSet;

typedef struct {
	char *strName;			// rule name	
	bool *bPhone;			// whether the phone is in the set of phones of the rule
	unsigned int iRule;	// rule identifier
} PhoneticRule;

typedef vector<PhoneticRule*> VPhoneticRule;

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class PhoneticRulesManager {

	private:
	
		string m_strFile;							// file containing phonetic rules
		unsigned char m_iFileFormat;			// file format
		VPhoneticRule m_vPhoneticRule;		// vector of phonetic rules
		PhoneSet *m_phoneSet;

	public:
    
    	// constructor
		PhoneticRulesManager(const char *strFile, PhoneSet *phoneSet);
		
		// constructor
		PhoneticRulesManager();

		// destructor
		~PhoneticRulesManager();
		
		// load the phonetic rules (text format)
		void load();
		
		// store the phonetic rules (binary format)
		void store(FileOutput &file);
		
		// load the phonetic rules (binary format)
		static PhoneticRulesManager *load(FileInput &file, PhoneSet *phoneSet);	
		
		// destroy the phonetic rules
		void destroy();
		
		// print phonetic rules statistics
		void print();
		
		// return the phonetic rules
		inline VPhoneticRule *getRules() {
		
			return &m_vPhoneticRule;
		}
		
		// return a phonetic rule by its index
		inline PhoneticRule *getRule(unsigned int iIndex) {
		
			return m_vPhoneticRule[iIndex];
		}
};

};	// end-of-namespace

#endif
