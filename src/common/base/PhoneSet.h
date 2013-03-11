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


#ifndef PHONESET_H
#define PHONESET_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Global.h"

using namespace std;

#include <vector>
#include <string>
#include <map>

namespace Bavieca {

#define MAX_BASEPHONES						255		// maximum # of basephones, the 255 value is reserved for error checking
#define MAX_PHONETIC_SYMBOL_LENGTH			10 

#define PHONETIC_SYMBOL_SILENCE				"SIL"
#define CONTEXT_PADDING						"<>"

typedef struct {
	unsigned char iIndex;			// phone index
	string strPhone;				// phone name
	bool bContext;					// whether context modeling can be applied to the phone 
} Phone;

typedef vector<Phone*> VPhone;

class PhoneSet {

   private: 
      
      string m_strFile;
      vector<string> m_phones;
      map<string,int> m_mPhone;
      VPhone m_vPhone;
     
   public:
   
      //constructor
      PhoneSet(const char *strFile);
   
      //destructor
      ~PhoneSet(); 
      
      //load the phone set
      void load();
 
      // return the number of phones
      inline unsigned int size() {
		  
		  return (unsigned int)m_vPhone.size();
      }
		 
		// return whether the phone is a special phone (can not appear in lexical word transcriptions)
		bool isSpecialPhone(int iPhone);	
		
		inline const char *getStrPhone(int iPhone) {
		
			if (iPhone == (int)m_phones.size()) {
				return CONTEXT_PADDING;
			}
			
			return m_phones[iPhone].c_str();
		}
		
		// return the index of the given phonetic symbol (-1 in case it is not found)
		inline int getPhoneIndex(const char *strPhone) {
			
			map<string,int>::const_iterator it = m_mPhone.find(strPhone);
			if (it == m_mPhone.end()) {
				return -1;
			}
				
			return it->second;
		}	
		
		// return the index of the silence phonetic symbol
		inline int getPhoneIndexSilence() {
			
			map<string,int>::const_iterator it = m_mPhone.find(PHONETIC_SYMBOL_SILENCE);
			if (it == m_mPhone.end()) {
				return -1;
			}
				
			return it->second;
		}	
		
		// return whether context modeling affect the phone
		inline bool isPhoneContextModeled(unsigned char iPhone) {
		
			assert(iPhone < m_vPhone.size());
		
			return m_vPhone[iPhone]->bContext;
		}
		
		// print the phonetic symbol set
		void print();
		
};

};	// end-of-namespace

#endif


