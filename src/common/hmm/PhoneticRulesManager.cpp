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


#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"
#include "PhoneticRulesManager.h"
#include "PhoneSet.h"

namespace Bavieca {

// constructor
PhoneticRulesManager::PhoneticRulesManager(const char *strFile, PhoneSet *phoneSet)
{
	m_strFile = strFile;
	m_iFileFormat = FILE_FORMAT_TEXT;
	m_phoneSet = phoneSet;
}

// constructor
PhoneticRulesManager::PhoneticRulesManager()
{
}

// destructor
PhoneticRulesManager::~PhoneticRulesManager()
{
	destroy();
}

// load the phonetic rules (text format)
void PhoneticRulesManager::load() {

   FileInput file(m_strFile.c_str(),false);
   file.open();
   
   int iBasePhones = m_phoneSet->size();
   int iLine = 1;

	string strLine;
   unsigned int iIndex = 0;
   while(std::getline(file.getStream(),strLine)) { 

      // check if the line corresponds to a rule
      if (strLine.c_str()[0] != '$') {
      	++iLine;
         continue;
		}
		
   	char *strLineCopy = new char[strLine.length()+1];
   	strcpy(strLineCopy,strLine.c_str());

      char *str;
      PhoneticRule *phoneticRule = new PhoneticRule;
      phoneticRule->iRule = iIndex++;
      phoneticRule->bPhone = new bool[iBasePhones+1];			// include the context-padding symbol
      
      // get the rule's tokens (l-part and r-part)
      
      // initialization (include the context-padding symbol)
      for(int i=0 ; i < iBasePhones+1 ; ++i) {
         phoneticRule->bPhone[i] = false;
      }
		
		// extract the rule's name
      str = strtok(strLineCopy," \t\n");
      phoneticRule->strName = new char[strlen(&(str[1]))+1];
      strcpy(phoneticRule->strName,&(str[1])); 
      
      // iterate through the right elements (they are either rule names or terminal symbols (phones))
      str = strtok(NULL," \t\n");
      while(str != NULL) {
         if (str[0] != '$') {
         	int iPhone = m_phoneSet->getPhoneIndex(str);
         	if (iPhone == -1) {
         		BVC_WARNING << "unknown phonetic symbol \"" << str << "\" found in file: "
         			<< m_strFile.c_str() << " (line: " << iLine << ")"  << endl;
         		str = strtok(NULL," \t\n");
         		continue;
         	}
            phoneticRule->bPhone[iPhone] = true;
         }
         else {
            PhoneticRule *aux = NULL;
            for(VPhoneticRule::const_iterator it = m_vPhoneticRule.begin() ; it != m_vPhoneticRule.end() ; ++it) {
					if (strcmp(&(str[1]),(*it)->strName) == 0) {
                  aux = *it;
                  break;
               }
            }
            if (aux == NULL) {
            	BVC_ERROR << "phonetic rule " << (str+1) << " used before but not defined yet (line: " << iLine << ")" << endl;
            }
            // copy the phones
            for(int i=0 ; i < iBasePhones ; ++i) {
               if (aux->bPhone[i]) {
                  phoneticRule->bPhone[i] = true;
               }
            }
         }
         str = strtok(NULL," \t\n");
      }
      m_vPhoneticRule.push_back(phoneticRule);
      delete [] strLineCopy;
      ++iLine;
   }
   
   // create a phonetic rule for the context padding symbol
	PhoneticRule *phoneticRule = new PhoneticRule;
	phoneticRule->strName = new char[strlen("contextPadding")+1];
	strcpy(phoneticRule->strName,"contextPadding");
	phoneticRule->iRule = iIndex++;
	phoneticRule->bPhone = new bool[iBasePhones+1];
	for(int i=0 ; i < iBasePhones+1 ; ++i) {
		(i == iBasePhones) ? phoneticRule->bPhone[i] = true : phoneticRule->bPhone[i] = false;
	}
	m_vPhoneticRule.push_back(phoneticRule);

	file.close();
   
   // check that all the phonetic symbols are represented by at least one rule
   bool *bPhone = new bool[iBasePhones+1];
   for(int i=0 ; i<iBasePhones+1 ; ++i) {
   	bPhone[i] = false;
   }
   for(VPhoneticRule::iterator it = m_vPhoneticRule.begin() ; it != m_vPhoneticRule.end() ; ++it) {
		for(int i=0 ; i<iBasePhones+1 ; ++i) {
			if ((*it)->bPhone[i]) {
				bPhone[i] = true;
			}
		}
   }
   for(int i=0 ; i<iBasePhones+1 ; ++i) {
   	if (bPhone[i] == false) {
   		BVC_WARNING << "phonetic symbol " << m_phoneSet->getStrPhone(i) << " not represented by any rule";
   	}
   }
   delete [] bPhone;
}

// store the phonetic rules (binary format)
void PhoneticRulesManager::store(FileOutput &file) {

   // (1) phonetic symbols
   int iPhones = m_phoneSet->size();
   IOBase::write(file.getStream(),iPhones); 
   for(int i=0 ; i<iPhones ; ++i) {
   	IOBase::writeString(file.getStream(),m_phoneSet->getStrPhone(i),(int)strlen(m_phoneSet->getStrPhone(i)));
   }

   // (2) phonetic rules
   int iRules = (int)m_vPhoneticRule.size();
	IOBase::write(file.getStream(),iRules);
   unsigned int iRule = 0;
   for(VPhoneticRule::iterator it = m_vPhoneticRule.begin() ; it != m_vPhoneticRule.end() ; ++it, ++iRule) {
   	IOBase::writeString(file.getStream(),(*it)->strName,(int)strlen((*it)->strName));
		// rule data (includes the context-padding symbol)
		for(int i=0 ; i<iPhones+1 ; ++i) {
			unsigned char iAnswer = ((*it)->bPhone[i]) ? PHONETIC_RULE_YES : PHONETIC_RULE_NO;
			IOBase::write(file.getStream(),iAnswer);
		} 
		assert((*it)->iRule == iRule);
   }
}

// load the phonetic rules (binary format)
PhoneticRulesManager *PhoneticRulesManager::load(FileInput &file, PhoneSet *phoneSet) {

   PhoneticRulesManager *phoneticRulesManager = new PhoneticRulesManager();
   phoneticRulesManager->m_strFile = "";
   phoneticRulesManager->m_phoneSet = phoneSet;
   phoneticRulesManager->m_iFileFormat = FILE_FORMAT_BINARY;

	// (1) phonetic symbols
   int iPhones = -1;
   IOBase::read(file.getStream(),&iPhones); 
   assert(iPhones == (int)phoneSet->size());
   char *strPhone = NULL;
   for(int i=0 ; i<iPhones ; ++i) {
   	IOBase::readString(file.getStream(),&strPhone);
   	assert(!strcmp(strPhone,phoneSet->getStrPhone(i)));
	   delete [] strPhone;
   }
 
	// (2) phonetic rules
   int iRules = -1;
	IOBase::read(file.getStream(),&iRules);
	assert(iRules > 0);
   for(int i=0 ; i<iRules ; ++i) {	
		PhoneticRule *rule = new PhoneticRule();
		rule->iRule = i;
   	IOBase::readString(file.getStream(),&rule->strName);
   	rule->bPhone = new bool[iPhones+1];
		// rule data (includes the context-padding symbol)
		for(int i=0 ; i<iPhones+1 ; ++i) {
			unsigned char iAnswer = 0;
			IOBase::read(file.getStream(),&iAnswer);
			(iAnswer == PHONETIC_RULE_YES) ? rule->bPhone[i] = true : rule->bPhone[i] = false;
		}
		phoneticRulesManager->m_vPhoneticRule.push_back(rule); 
   }

	return phoneticRulesManager;
}

// destroy the phonetic rules
void PhoneticRulesManager::destroy() {

	for(VPhoneticRule::iterator it = m_vPhoneticRule.begin() ; it != m_vPhoneticRule.end() ; ++it) {
		delete [] (*it)->strName;
		delete [] (*it)->bPhone;
		delete (*it);
	}
}

// print phonetic rules statistics
void PhoneticRulesManager::print() {

	cout << m_vPhoneticRule.size() << " phonetic rules loaded from file: " << m_strFile.c_str() << endl;
	for(VPhoneticRule::iterator it = m_vPhoneticRule.begin() ; it != m_vPhoneticRule.end() ; ++it) {
		cout << (*it)->strName;
		for(unsigned int i=0 ; i <= m_phoneSet->size() ; ++i) {	
			if ((*it)->bPhone[i]) {
				cout << " " << m_phoneSet->getStrPhone(i);
			}
		}
		cout << endl;
	}
}

};	// end-of-namespace

