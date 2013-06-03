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


#include "FileUtils.h"
#include "LogMessage.h"

namespace Bavieca {

// return the folder of a given file path (path = folder + filename)
void FileUtils::getFolder(char *strFolder, const char *strPath) {

	const char *strPointer = strrchr(strPath,PATH_SEPARATOR);
	if (strPointer == NULL) {
		sprintf(strFolder,".%c",PATH_SEPARATOR);
	} else if (strPointer == strPath) {
		sprintf(strFolder,"%c",PATH_SEPARATOR);
	} else {
		int iLength = (int)(strPointer-strPath);
		assert(iLength>=0);
		strncpy(strFolder,strPath,iLength);
		strFolder[iLength] = 0;
	}
}

// return the filename from a given path
const char *FileUtils::getFileName(const char *strPath) {

	const char *strPointer = strrchr(strPath,PATH_SEPARATOR);
	if (strPointer == NULL) {
		return strPath;
	} else {
		return strPointer+1;
	}	
}

// return the filename from a given path (without the file extension)
void FileUtils::getFileNameWithoutExtension(char *strFileName, const char *strPath) {

	const char *strPointerSeparator = strrchr(strPath,PATH_SEPARATOR);
	const char *strPointerExtension = strrchr(strPath,'.');
	if (strPointerSeparator == NULL) {
		if (strPointerExtension == NULL) {
			strcpy(strFileName,strPath);
		} else {
			unsigned int iLength = (unsigned int)(strPointerExtension-strPath);
			strncpy(strFileName,strPath,iLength);
			strFileName[iLength] = 0;
		}	
	} else {
		if (strPointerExtension == NULL) {
			unsigned int iLength = (unsigned int)(strlen(strPath)-(strPointerSeparator-strPath)-1);
			strncpy(strFileName,strPointerSeparator+1,iLength);
			strFileName[iLength] = 0;
		} else {
			unsigned int iLength = (unsigned int)(strPointerExtension-strPointerSeparator-1);
			strncpy(strFileName,strPointerSeparator+1,iLength);
			strFileName[iLength] = 0;
		}	
	}
}

// return the extension from a given path
const char *FileUtils::getFileExtension(const char *strPath) {

	const char *strPointer = strrchr(strPath,'.');
	if (strPointer == NULL) {
		return NULL;
	} else {
		return strPointer+1;
	}	
}

// returns a copy of the given path with the new extension (the extension is everything after the last period)
void FileUtils::replaceExtension(char *strDest, const char *strPath, const char *strExtension) {

	const char *strPointer = strrchr(strPath,'.');
	if (strPointer == NULL) {
		sprintf(strDest,"%s.%s",strPath,strExtension);
	} else {	
		int iChars = (int)(strPointer-strPath+1);
		strncpy(strDest,strPath,iChars);
		strcpy(strDest+iChars,strExtension);
	}	
}

// returns a copy of the given path with the new folder (the folder is everything before the last path separator)
void FileUtils::replaceFolder(char *strDest, const char *strPath, const char *strFolder) {

	const char *strPointer = strrchr(strPath,PATH_SEPARATOR);
	if (strPointer == NULL) {
		sprintf(strDest,"%s%c%s",strFolder,PATH_SEPARATOR,strPath);
	} else {
		sprintf(strDest,"%s%s",strFolder,strPointer);
	}	
}

// truncate a file (set it size to zero)
void FileUtils::truncateFile(const char *strFile) {

	// truncate the hypothesis file in case it exists
	FILE *file = fopen(strFile,"wb");
	if (file) {
		if (fclose(file) == EOF) {
			BVC_ERROR << "Unable to truncate the file: " << strFile;
		}
	}
}

// return whether the file exists
bool FileUtils::isFile(const char *strFile) {

	struct stat stFileInfo;

	// get the file attributes
	int iStat = stat(strFile,&stFileInfo);
	if (iStat == 0) {
		// unable to get the file attributes: the file does not exist
		return true;
	} 

	return false;
}

#if defined __linux__ || defined __APPLE__ || __MINGW32__

// create the given folder
int FileUtils::createFolder(const char *strFolder, mode_t mode) {

	struct stat st;

	// check if the folder exists
	if (stat(strFolder, &st) != 0) {
		// the folder does not exist: create it		
		if (mkdir(strFolder, mode) != 0) {
			return RETURN_CODE_ERROR;
		}
	}
	// folder exists already!
	else if (!S_ISDIR(st.st_mode)) {
		return RETURN_CODE_NOT_FOLDER;
	}
			
	return RETURN_CODE_SUCCESS;
}

// create the given path 
// it goes from left to right trying to create the folder
int FileUtils::createPath(const char *strPath, mode_t mode) {
	
	int iReturnCode = RETURN_CODE_SUCCESS;
	char *strAux = new char[strlen(strPath)+1];
	strcpy(strAux,strPath);
		
	char *strP1 = strAux;
	char *strP2;
	while((iReturnCode == RETURN_CODE_SUCCESS) && ((strP2 = strchr(strP1,'/')) != 0)) {
		if (strP2 != strP1) {
			*strP2 = '\0';
			iReturnCode = createFolder(strAux,mode);
			*strP2 = '/';
		}
		strP1 = strP2 + 1;
	}
	if (iReturnCode == RETURN_CODE_SUCCESS) {
		iReturnCode = createFolder(strPath);
	}	
	delete [] strAux;

	return iReturnCode;
}

#endif

};	// end-of-namespace
