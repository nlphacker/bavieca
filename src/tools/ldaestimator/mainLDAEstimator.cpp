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


#include "CommandLineManager.h"
#include "LDAEstimator.h"

using namespace Bavieca;
 
// main for the tool "ldaestimator"
int main(int argc, char *argv[]) {

	// define the parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("hldaestimator",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-acc","input accumulator filelist",PARAMETER_TYPE_FILE,true);	
	m_commandLineManager->defineParameter("-dim","target dimensionality",PARAMETER_TYPE_INTEGER,false);	
	m_commandLineManager->defineParameter("-out","output transform",PARAMETER_TYPE_FILE,false);	
	
	// parse the parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
		
	// fixing-bugs

	return 0;
}

