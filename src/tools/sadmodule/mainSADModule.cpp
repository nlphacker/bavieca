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
#include "FeatureFile.h"
#include "FeatureExtractor.h"
#include "Global.h"
#include "HMMManager.h"
#include "PhoneSet.h"
#include "SADModule.h"
#include "TimeUtils.h"

using namespace std;

#include <string>
#include <map>

using namespace Bavieca;

// main for the Speech Activity Detection tool: "sadmodule"
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("sadmodule",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-fea","feature file to process",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-sil","max Gaussian components for silence",
			PARAMETER_TYPE_INTEGER,false);
		commandLineManager.defineParameter("-sph","max Gaussian components for speech",
			PARAMETER_TYPE_INTEGER,false);
		commandLineManager.defineParameter("-pad","speech padding (# frames)",PARAMETER_TYPE_INTEGER,true,NULL,"10");	
		commandLineManager.defineParameter("-pen","speech insertion penalty",PARAMETER_TYPE_FLOAT,false);	
		commandLineManager.defineParameter("-out","speech segmentation output",PARAMETER_TYPE_FILE,false);
		
		// (2) parse command line parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// get parameters
		const char *strFileFeatures = commandLineManager.getParameterValue("-fea");
		const char *strFilePhoneSet = commandLineManager.getParameterValue("-pho");
		const char *strFileModels = commandLineManager.getParameterValue("-mod");	
		int iMaxGaussianComponentsSilence = atoi(commandLineManager.getParameterValue("-sil"));	
		int iMaxGaussianComponentsSpeech = atoi(commandLineManager.getParameterValue("-sph"));	
		int iFramesPadding = atoi(commandLineManager.getParameterValue("-pad"));	
		float fPenaltySilenceToSpeech = atof(commandLineManager.getParameterValue("-pen"));	
		const char *strFileOutput = commandLineManager.getParameterValue("-out");	
		
		// load the phone set
		PhoneSet phoneSet(strFilePhoneSet);
		phoneSet.load();
		
		// load the acoustic models
		HMMManager hmmManager(&phoneSet,HMM_PURPOSE_EVALUATION);
		hmmManager.load(strFileModels);
		
		// load the file to process
		FeatureFile featureFile(strFileFeatures,MODE_READ);
		featureFile.load();
		
		unsigned int iFeatures = 0;
		float *fFeatures = featureFile.getFeatureVectors(&iFeatures);
		
		// create the SAD module
		SADModule sadModule(&phoneSet,&hmmManager,iMaxGaussianComponentsSilence,
			iMaxGaussianComponentsSpeech,fPenaltySilenceToSpeech,iFramesPadding);
		sadModule.initialize();
		
		double dTimeBegin = TimeUtils::getTimeMilliseconds();
	
		// perform the actual segmentation
		sadModule.beginSession();	
		sadModule.processFeatures(fFeatures,iFeatures);	
		VSpeechSegment vSpeechSegment;
		sadModule.recoverSpeechSegments(vSpeechSegment);
		//sadModule.printSegments(vSpeechSegment);	
		sadModule.endSession();
		
		double dTimeEnd = TimeUtils::getTimeMilliseconds();	
		double dSecondsProcessing = (dTimeEnd-dTimeBegin)/1000.0;
		double dSecondsAudio = ((float)iFeatures)/100.0;
		printf("feature frames: %u (RTF: %.2f)\n",iFeatures,dSecondsProcessing/dSecondsAudio);
		
		// write the segmentation to file
		SADModule::store(strFileOutput,vSpeechSegment);
		
		// clean-up
		SADModule::deleteVSpeechSegment(vSpeechSegment);
		delete [] fFeatures;
		
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
	return 0;
}
