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

#include "MLFFile.h"

// constructor
MLFFile::MLFFile(LexiconManager *lexiconManager, const char *strFile, Log *log, const unsigned char iMode)
{
	m_lexiconManager = lexiconManager;
	m_strFile = strFile;
	m_log = log;
	m_iMode = iMode;
}

// destructor
MLFFile::~MLFFile()
{
	for(VMLFUtterance::iterator it =  m_vMLFUtterance.begin() ; it != m_vMLFUtterance.end() ; ++it) {
		delete *it;
	}
	m_vMLFUtterance.clear();
}

// load the MLF file from disk
void MLFFile::load() {
   
   string strLine;
   int iLine = 0;
   
   assert(m_iMode == MODE_READ);
	
	FileInput file(m_strFile.c_str(),false);
	file.open();
	
	// (2) for each utterance, read the feature file name and the associated lexical units
	MLFUtterance *mlfUtterance = NULL;
   while(std::getline(file.getStream(),strLine).good()) {
   	++iLine;
   	// feature file?
   	if (isFeatureFilename(strLine)) {
   		if (mlfUtterance != NULL) {
   			if (mlfUtterance->vLexUnit.empty() == true) {
   				ERROR << "loading MLF at line " << iLine << "unexpected feature file name";
   			}
   			m_vMLFUtterance.push_back(mlfUtterance);
   		}
   		mlfUtterance = new MLFUtterance;
   		mlfUtterance->strFilePattern = strLine.substr(1,strLine.length()-2);	
   	} 
   	// word?
   	else {
   		if (mlfUtterance == NULL) {
   			ERROR << "loading MLF at line " << iLine << ", unexpected lexical unit found or wrong format of feature file name";
   		}
   		LexUnit *lexUnit = parseLexUnitLine(strLine,iLine);
   		assert(lexUnit);
   		mlfUtterance->vLexUnit.push_back(lexUnit);	
   	}
   }	
	if (mlfUtterance != NULL) {
		if (mlfUtterance->vLexUnit.empty() == true) {	
			ERROR << "loading MLF: line " << iLine;
		}
		m_vMLFUtterance.push_back(mlfUtterance);
	}
	
	file.close();
}

// store the MLF utterances into disk
void MLFFile::store() {

	assert(m_iMode == MODE_WRITE);
	
	FileOutput file(m_strFile.c_str(),false);
	file.open();
	
	// (2) for each utterance, write the feature file name and the associated lexical units
	for(VMLFUtterance::iterator it = m_vMLFUtterance.begin() ; it != m_vMLFUtterance.end() ; ++it) {
		// write the feature file name
		IOBase::writeString(file.getStream(),(*it)->strFilePattern);
		// write each of the lexical units
		for(VLexUnit::const_iterator jt = (*it)->vLexUnit.begin() ; jt != (*it)->vLexUnit.end() ; ++jt) {
			if ((*jt)->iPronunciation == 0) {
				IOBase::writeString(file.getStream(),(*it)->strFilePattern);
			} else {
				assert((*jt)->iPronunciation >= 1);
				std::ostringstream oss;
				oss << (*it)->strFilePattern << "(" << ((*jt)->iPronunciation+1) << ")";
				string str = oss.str();
				IOBase::writeString(file.getStream(),str);
			}
		}
	}
	
	file.close();
}

// return whether the given string is a feature filename
bool MLFFile::isFeatureFilename(string &strFile) {

	if ((strFile.find_first_of("\"") != 0) ||
		(strFile.find_last_of(".fea\"") != strFile.length()-5)) {
		return false;	
	}

	return true;
}

// extract the lexical unit identity and the pronunciation number
LexUnit *MLFFile::parseLexUnitLine(string &strLexUnit, int iMLFLine) {

   string strLexUnitBase = strLexUnit;
   int iPronunciation = 0;
   
   size_t iFirst = strLexUnit.find_first_of("(");
   if (iFirst != -1) {
		size_t iLast = strLexUnit.find_last_of(")");
		if (iLast-iFirst < 2) {
			ERROR << "invalid lexical unit \"" << strLexUnit << "\" read from the MLF: line " << iMLFLine;
		}
		iPronunciation = atoi(strLexUnit.substr(iFirst+1,iLast-iFirst+1).c_str());
		strLexUnitBase = strLexUnit.substr(0,iFirst);
	} 
		
	// get the lexical unit from the lexicon
	LexUnit *lexUnit = m_lexiconManager->getLexUnit(strLexUnitBase.c_str(),iPronunciation);
	if (lexUnit == NULL) {
		ERROR << "lexical unit " << strLexUnit << " read from the MLF was not found in the lexicon: line " << iMLFLine;
	}

	return lexUnit;
}
