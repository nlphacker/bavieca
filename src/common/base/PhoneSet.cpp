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


#include "PhoneSet.h"
#include "FileInput.h"
#include "IOBase.h"
#include "ctype.h"
#include "LogMessage.h"

namespace Bavieca {

//constructor
PhoneSet::PhoneSet(const char* strFile) {
   
   m_strFile = strFile;
}

// destructor
PhoneSet::~PhoneSet() {

	for(VPhone::iterator it = m_vPhone.begin() ; it != m_vPhone.end() ; ++it) {
		delete *it;
	}
}
   
// load the phone set
void PhoneSet::load() {

	FileInput file(m_strFile.c_str(),false);
	file.open();	
	
	string strLine;
	int iLine = 0;
	while(std::getline(file.getStream(),strLine)) {
		
		++iLine;
	
		// skip comments and blank lines
		if (strLine.empty() || (strLine.c_str()[0] == '#')) {
			continue;
		}
		
		// read the phone
		std::stringstream s(strLine);
		string strPhone;
		IOBase::readString(s,strPhone);	
		
		Phone *phone = new Phone;
		phone->bContext = true;
		if ((strPhone.c_str()[0] == '(') && (strPhone.c_str()[strPhone.length()-1] == ')')) {	
			phone->bContext = false;
			strPhone = strPhone.substr(1,strPhone.length()-2);
		}	
		
		// check length
		if (strPhone.length() > MAX_PHONETIC_SYMBOL_LENGTH) {
			BVC_ERROR << "incorrect phonetic symbol name \"" << strPhone << "\" found in line " << iLine << ", too many characters";
		}
		
		// check that the phone name is correct
		for(unsigned int i=0 ; i < strPhone.length() ; ++i) {
			if (!isalnum(strPhone.c_str()[i]) && (strPhone.c_str()[i] != '_')) {
				BVC_ERROR << "incorrect phonetic symbol name \"" << strPhone << "\" found in line " << iLine << ", wrong symbol";
			}
		} 
      m_phones.push_back(strPhone);
      m_mPhone.insert(map<string,int>::value_type(strPhone,(int)(m_phones.size()-1)));
   
      phone->iIndex = (unsigned char)(m_phones.size()-1);
      phone->strPhone = strPhone;
      m_vPhone.push_back(phone);	
	}
	
	file.close();
   
   // check that there is at least one phone
   if (m_phones.empty()) {
   	BVC_ERROR << "no phonetic symbols found in the phonetic symbol file";
   }
   
   // check that the silence is defined
   if (m_mPhone.find(PHONETIC_SYMBOL_SILENCE) == m_mPhone.end()) {
   	BVC_ERROR << "phonetic symbol " << PHONETIC_SYMBOL_SILENCE << " must be defined in the phonetic symbol file";
   }	
}

// print the phonetic symbol set
void PhoneSet::print() {

	cout << "-- phonetic symbol set ---------------------------------------------------\n";
	cout << " file: " << m_strFile << "\n";
	cout << " #phones: " << m_vPhone.size() << "\n";
	for(VPhone::iterator it = m_vPhone.begin() ; it != m_vPhone.end() ; ++it) {
		cout << static_cast<int>((*it)->iIndex) << " " << (*it)->strPhone;
		if ((*it)->bContext == false) {
			cout << " (context independent)\n";
		} else {
			cout << "\n";
		}
	}
	cout << "--------------------------------------------------------------------------\n";
}	

};	// end-of-namespace

