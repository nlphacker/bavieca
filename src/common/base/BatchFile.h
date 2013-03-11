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


#ifndef VTLNBATCHFILE_H
#define VTLNBATCHFILE_H

#include "Global.h"

using namespace std;

#include <string>
#include <vector>
#include <map>

namespace Bavieca {

typedef struct {
	vector<string> vStrElement;	
} BatchEntry;

typedef vector<BatchEntry*> VBatchEntry;


/**
	@author daniel <dani.bolanos@gmail.com>
*/
class BatchFile {

	private:
	
		string m_strFile;							// filename
		VBatchEntry m_vBatchEntry;				// entries
		string m_strType;							// batch file type
		unsigned int m_iColumns;							// columns per entry
		map<string,int> m_mColumnName;		// maps column names to column indices

	public:

		// constructor
		BatchFile(const char *strFile, const char *strType);

		// destructor
		~BatchFile();
		
		// load the content of the batch file
		void load();
		
		// return the number of entries in the batch file
		unsigned int size() {
		
			return (unsigned int)m_vBatchEntry.size();
		}
		
		// return the field in the given entry and column
		const char *getField(unsigned int iEntry, unsigned int iColumn);
		
		// return the field in the given entry by its name
		const char *getField(unsigned int iEntry, const char *strName);

};

};	// end-of-namespace

#endif
