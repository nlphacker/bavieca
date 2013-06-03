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


#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>

namespace Bavieca {

#if defined __linux__ || defined __APPLE__ || __MINGW32__
#define PATH_SEPARATOR		'/'
#elif _MSC_VER
#define PATH_SEPARATOR	'\\'
#else
	#error "unsupported platform"
#endif

// return codes for createFolder
#define RETURN_CODE_ERROR					0
#define RETURN_CODE_SUCCESS				1
#define RETURN_CODE_NOT_FOLDER			2

/**
	@author daniel <dani.bolanos@gmail.com>
*/
class FileUtils {

	public:
   
		// return the folder of a given file path (path = folder + filename)
		static void getFolder(char *strFolder, const char *strPath);
		
		// return the filename from a given path
		static const char *getFileName(const char *strPath);
		
		// return the filename from a given path (without the file extension)
		static void getFileNameWithoutExtension(char *strFileName, const char *strPath);
		
		// return the extension from a given path
		static const char *getFileExtension(const char *strPath);	
		
		// returns a copy of the given path with the new extension (the extension is everything after the last period)
		static void replaceExtension(char *strDest, const char *strPath, const char *strExtension);
		
		// returns a copy of the given path with the new folder (the folder is everything before the last path separator)
		static void replaceFolder(char *strDest, const char *strPath, const char *strFolder);
		
		// truncate a file (set it size to zero)
		static void truncateFile(const char *strFile);

		// return whether the file exists
		static bool isFile(const char *strFile);

	#if defined __linux__ || defined __APPLE__ || __MINGW32__
		
		// create the given folder
		static int createFolder(const char *strFolder, mode_t mode = 0777);
		
		// create the given path 
		// it goes from left to right trying to create the folder
		static int createPath(const char *strPath, mode_t mode = 0777);

	#endif
		
};

};	// end-of-namespace

#endif
