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

	// (1) define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("sadmodule",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-pho","phonetic symbol set",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-mod","acoustic models",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-fea","feature file to process",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineParameter("-sil","max Gaussian components for silence",
		PARAMETER_TYPE_INTEGER,false);
	m_commandLineManager->defineParameter("-sph","max Gaussian components for speech",
		PARAMETER_TYPE_INTEGER,false);
	m_commandLineManager->defineParameter("-pad","speech padding (# frames)",PARAMETER_TYPE_INTEGER,true,NULL,"10");	
	m_commandLineManager->defineParameter("-pen","speech insertion penalty",PARAMETER_TYPE_FLOAT,false);	
	m_commandLineManager->defineParameter("-out","speech segmentation output",PARAMETER_TYPE_FILE,false);
	
	// (2) parse command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get parameters
	const char *m_strFileFeatures = m_commandLineManager->getParameterValue("-fea");
	const char *m_strFilePhoneSet = m_commandLineManager->getParameterValue("-pho");
	const char *m_strFileModels = m_commandLineManager->getParameterValue("-mod");	
	int m_iMaxGaussianComponentsSilence = atoi(m_commandLineManager->getParameterValue("-sil"));	
	int m_iMaxGaussianComponentsSpeech = atoi(m_commandLineManager->getParameterValue("-sph"));	
	int m_iFramesPadding = atoi(m_commandLineManager->getParameterValue("-pad"));	
	float m_fPenaltySilenceToSpeech = atof(m_commandLineManager->getParameterValue("-pen"));	
	const char *m_strFileOutput = m_commandLineManager->getParameterValue("-out");	
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneSet);
   m_phoneSet->load();
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	
	// load the file to process
	FeatureFile *m_featureFile = new FeatureFile(m_strFileFeatures,MODE_READ);
	m_featureFile->load();
	
	int m_iFeatures = -1;
	float *m_fFeatures = m_featureFile->getFeatureVectors(&m_iFeatures);
	
	// create the SAD module
	SADModule *m_sadModule = new SADModule(m_phoneSet,m_hmmManager,m_iMaxGaussianComponentsSilence,
		m_iMaxGaussianComponentsSpeech,m_fPenaltySilenceToSpeech,m_iFramesPadding);
	m_sadModule->initialize();
	
	double dTimeBegin = TimeUtils::getTimeMilliseconds();

	// perform the actual segmentation
	m_sadModule->beginSession();	
	m_sadModule->processFeatures(m_fFeatures,m_iFeatures);	
	VSpeechSegment m_vSpeechSegment;
	m_sadModule->recoverSpeechSegments(m_vSpeechSegment);
	//m_sadModule->printSegments(m_vSpeechSegment);	
	m_sadModule->endSession();
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();	
	double dSecondsProcessing = (dTimeEnd-dTimeBegin)/1000.0;
	double dSecondsAudio = ((float)m_iFeatures)/100.0;
	printf("feature frames: %u (RTF: %.2f)\n",m_iFeatures,dSecondsProcessing/dSecondsAudio);
	
	// write the segmentation to file
	SADModule::store(m_strFileOutput,m_vSpeechSegment);
	
	// clean-up
	SADModule::deleteVSpeechSegment(m_vSpeechSegment);
	delete m_sadModule;
	delete m_featureFile;
	delete [] m_fFeatures;
	delete m_hmmManager;
	delete m_phoneSet;	
	delete m_commandLineManager;
	
	return 0;
}
