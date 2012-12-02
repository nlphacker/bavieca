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


#include "AudioFile.h"
#include "BatchFile.h"
#include "BestPath.h"
#include "ConfigurationFeatures.h"
#include "ConfigurationWFSADecoder.h"
#include "CommandLineManager.h"
#include "FeatureExtractor.h"
#include "FeatureFile.h"
#include "FileUtils.h"
#include "LogMessage.h"
#include "SADModule.h"
#include "Transform.h"
#include "TimeUtils.h"
#include "Viterbi.h"
#include "WFSADecoder.h"

using namespace Bavieca;

// main for the tool "wfsdecoder"
int main(int argc, char *argv[]) {

	// define command line parameters
	CommandLineManager *m_commandLineManager = new CommandLineManager("wfsadecoder",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
	m_commandLineManager->defineParameter("-cfg","configuration",PARAMETER_TYPE_FILE,false);	
	m_commandLineManager->defineParameter("-bat","batch file containing entries [rawFile/featureFile utteranceId]",PARAMETER_TYPE_FILE,false);
	m_commandLineManager->defineIncompatibility("-list|-audio",true);
	m_commandLineManager->defineParameter("-hyp","hypothesis",PARAMETER_TYPE_FILE,false);

	// parse command line parameters
	if (m_commandLineManager->parseParameters(argc,argv) == false) {
		return -1;
	}
	
	// get the command-line parameters
	const char *m_strFileConfiguration = m_commandLineManager->getParameterValue("-cfg");
	
	// load the configuration file
	ConfigurationWFSADecoder *m_configuration = new ConfigurationWFSADecoder(m_strFileConfiguration);
	m_configuration->load();
	
	// get the configuration parameters
	
	// decoding network
	const char *m_strFileNetwork = m_configuration->getStrParameterValue("decodingNetwork.file");
	
	// input
	const char *m_strInputType = m_configuration->getStrParameterValue("input.type");	

	// phone-set
	const char *m_strFilePhoneticSymbolSet = m_configuration->getStrParameterValue("phoneticSymbolSet.file");	
	
	// feature extraction
	const char *m_strFileConfigurationFeatures = 
		m_configuration->getStrParameterValue("feature.configurationFile");
	const char *m_strCepstralNormalizationMode =
		m_configuration->getStrParameterValue("feature.cepstralNormalization.mode");
	const char *m_strCepstralNormalizationMethod =
		m_configuration->getStrParameterValue("feature.cepstralNormalization.method");
	int m_iCepstralNormalizationBufferSize = 
		m_configuration->getIntParameterValue("feature.cepstralNormalization.bufferSize");
	float m_fWarpFactor = 1.0;
	if (m_configuration->isParameterSet("feature.warpFactor")) {
		m_fWarpFactor = atof(m_configuration->getParameterValue("feature.warpFactor"));
	}
	const char *m_strFileFeatureTransform = NULL;
	if (m_configuration->isParameterSet("feature.transformFile")) {
		m_strFileFeatureTransform = m_configuration->getParameterValue("feature.transformFile");
	}
	
	// acoustic models
	const char *m_strFileModels = m_configuration->getStrParameterValue("acousticModels.file");

	// lexicon	
	const char *m_strFileLexicon = m_configuration->getStrParameterValue("lexicon.file");
	
	// pruning
	int m_iMaxActiveStates = m_configuration->getIntParameterValue("pruning.maxActiveStates");
	float m_fLikelihoodBeam = m_configuration->getFloatParameterValue("pruning.likelihoodBeam");	
	
	// output lattice?
	bool m_bLatticeGeneration = m_configuration->isParameterSet("output.lattice.folder");	
	const char *m_strFolderLattices = NULL;
	int m_iMaxWordSequencesState = -1;
	if (m_bLatticeGeneration) {
		m_strFolderLattices = m_configuration->getStrParameterValue("output.lattice.folder"); 
		m_iMaxWordSequencesState = 
			m_configuration->getIntParameterValue("output.lattice.maxWordSequencesState");
	}
	
	// output features?
	bool m_bOutputFeatures = m_configuration->isParameterSet("output.features.folder");
	const char *m_strFolderFeatures = NULL;
	if (m_bOutputFeatures) {
		m_strFolderFeatures = m_configuration->getStrParameterValue("output.features.folder"); ;
	}
	
	// output alignment?
	bool m_bOutputAlignment = m_configuration->isParameterSet("output.alignment.folder");
	const char *m_strFolderAlignments = NULL;
	if (m_bOutputAlignment) {
		m_strFolderAlignments = m_configuration->getStrParameterValue("output.alignment.folder"); ;
	}

	// output audio?
	bool m_bOutputAudio = m_configuration->isParameterSet("output.audio.folder");
	const char *m_strFolderAudio = NULL;
	if (m_bOutputAudio) {
		m_strFolderAudio = m_configuration->getStrParameterValue("output.audio.folder");
	}
	
	// (3) initialize the system
	
   // load the phone set
   PhoneSet *m_phoneSet = new PhoneSet(m_strFilePhoneticSymbolSet);
   m_phoneSet->load();
   
   // load the feature transforms
   VTransform m_vTransformFeatures;
   if (m_strFileFeatureTransform) {
		BatchFile *m_batchFile = new BatchFile(m_strFileFeatureTransform,"transform");
		m_batchFile->load();
		for(unsigned int i=0 ; i < m_batchFile->size() ; ++i) {
			Transform *transform = new Transform();
			transform->load(m_batchFile->getField(i,0u));
			m_vTransformFeatures.push_back(transform);
		}	
		delete m_batchFile;	
   }   
   
	// load the acoustic models
	HMMManager *m_hmmManager = new HMMManager(m_phoneSet,HMM_PURPOSE_EVALUATION);
	m_hmmManager->load(m_strFileModels);
	m_hmmManager->initializeDecoding();	
   
   // load the lexicon
   LexiconManager *m_lexiconManager = new LexiconManager(m_strFileLexicon,m_phoneSet);
   m_lexiconManager->load();
   m_lexiconManager->print(false);
	
	// create the aligner object (if needed)
	Viterbi *m_viterbi = NULL;
	if (m_bOutputAlignment == true) {
		m_viterbi = new Viterbi(m_phoneSet,m_hmmManager,m_lexiconManager);
	}
	
   // load the feature configuration
   ConfigurationFeatures *configurationFeatures = new ConfigurationFeatures(m_strFileConfigurationFeatures);
   configurationFeatures->load();
   
   // get the feature normalization mode and method
   int m_iCepstralNormalizationMode = FeatureExtractor::getNormalizationMode(m_strCepstralNormalizationMode);
   int m_iCepstralNormalizationMethod = FeatureExtractor::getNormalizationMethod(m_strCepstralNormalizationMethod);
   
   // create the feature extractor
   FeatureExtractor *m_featureExtractor = new FeatureExtractor(configurationFeatures,m_fWarpFactor,
   	m_iCepstralNormalizationBufferSize,m_iCepstralNormalizationMode,m_iCepstralNormalizationMethod);
   m_featureExtractor->initialize();

	// load the decoding network
	WFSAcceptor *m_wfsAcceptor = WFSAcceptor::load(m_lexiconManager,m_strFileNetwork);
	assert(m_wfsAcceptor);

	// create the decoding object
	WFSADecoder *m_wfsaDecoder = new WFSADecoder(m_phoneSet,m_hmmManager,m_lexiconManager,
		m_wfsAcceptor,m_iMaxActiveStates,m_fLikelihoodBeam,m_bLatticeGeneration,m_iMaxWordSequencesState);
	m_wfsaDecoder->initialize();
	
   string m_strFileHypothesis = m_commandLineManager->getParameterValue("-hyp");
	
   double dLikelihoodTotal = 0.0;
	unsigned int iFeatureVectorsTotal = 0;
	unsigned int iUtterances = 0;
   double dTimeBegin = TimeUtils::getTimeMilliseconds();
   
   BatchFile *batchFile = NULL;
	
	// BATCH MODE (batch file containing pairs [ featureFile  utteranceID ])
	if (m_commandLineManager->isParameterSet("-bat")) {
	
		// best single path options
		bool bBestSinglePathOutputFillers = false;
		bool bBestSinglePathOutputConfidenceValues = false;
		bool bBestSinglePathOutputSentenceDelimiters = false;
	
		string strFileList = m_commandLineManager->getParameterValue("-bat");
		VFeaturesUtterance vFeaturesUtterance;
		
		// audio input
		if (strcmp(m_strInputType,"audio") == 0) {
		
			// load the batch file
			batchFile = new BatchFile(strFileList.c_str(),"audio|id");
			batchFile->load();
			
			// (4) extract session-data
			VUtteranceData vUtteranceData;
			for(unsigned int iUtterance = 0 ; iUtterance < batchFile->size() ; ++iUtterance) {
			
				// load the raw audio
				int iSamples = -1;
				short int *sSamples = AudioFile::load(batchFile->getField(iUtterance,"audio"),&iSamples);
			
				UtteranceData utteranceData;
				utteranceData.samples.sSamples = sSamples;
				utteranceData.samples.iSamples = iSamples;
				utteranceData.features.fFeatures = NULL;
				utteranceData.features.iFeatures = -1;
				vUtteranceData.push_back(utteranceData);
			}	
			
			// extract features
			m_featureExtractor->extractFeaturesSession(vUtteranceData,true);
			
			// apply feture transforms
			int iDimFea = m_featureExtractor->getFeatureDimensionality();
			for(VTransform::iterator it = m_vTransformFeatures.begin() ; it != m_vTransformFeatures.end() ; ++it) {
				printf("%d -> %d\n",iDimFea,(*it)->getRows());
				for(VUtteranceData::iterator jt = vUtteranceData.begin() ; jt != vUtteranceData.end() ; ++jt) {
					float *fFeaturesX = new float[jt->features.iFeatures*(*it)->getRows()];
					for(int i=0 ; i < jt->features.iFeatures ; ++i) {
						(*it)->apply(jt->features.fFeatures+(i*iDimFea),fFeaturesX+(i*(*it)->getRows()));
					}
					delete [] jt->features.fFeatures;
					jt->features.fFeatures = fFeaturesX;
				}
				iDimFea = (*it)->getRows();
			}
			
			FeaturesUtterance featuresUtterance;
			for(VUtteranceData::iterator it = vUtteranceData.begin() ; it != vUtteranceData.end() ; ++it) {
				featuresUtterance.fFeatures = it->features.fFeatures;
				featuresUtterance.iFeatures = it->features.iFeatures;
				delete [] it->samples.sSamples;
				vFeaturesUtterance.push_back(featuresUtterance);	
			}
		}
		// feature input
		else {
		
			// load the batch file
			batchFile = new BatchFile(strFileList.c_str(),"features|id");
			batchFile->load();
		
			iUtterances = batchFile->size();
			for(unsigned int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) {
			
				string strUtteranceId = batchFile->getField(iUtterance,"id");
				string strFile = batchFile->getField(iUtterance,"features");
			
				FeatureFile *featureFile = new FeatureFile(strFile.c_str(),MODE_READ);
				featureFile->load();
				FeaturesUtterance featuresUtterance;
				featuresUtterance.fFeatures = featureFile->getFeatureVectors(&featuresUtterance.iFeatures);
				vFeaturesUtterance.push_back(featuresUtterance);	
				delete featureFile;
			}
		}
		
		// actual decoding
		FileOutput fileHypothesis(m_strFileHypothesis.c_str(),false);
		fileHypothesis.open();
		iUtterances = batchFile->size();
		for(unsigned int iUtterance = 0 ; iUtterance < iUtterances ; ++iUtterance) {
		
			string strUtteranceId = batchFile->getField(iUtterance,"id");
			printf("Processing utterance: %s\n",strUtteranceId.c_str());
			
			int iFeatures = vFeaturesUtterance[iUtterance].iFeatures;
			float *fFeatures = vFeaturesUtterance[iUtterance].fFeatures;;
			m_wfsaDecoder->viterbi(fFeatures,iFeatures);
			iFeatureVectorsTotal += iFeatures;
			
			// hypothesis lattice
			if (m_bLatticeGeneration) {
				HypothesisLattice *hypothesisLattice = m_wfsaDecoder->getHypothesisLattice();
				if (hypothesisLattice != NULL) {
					ostringstream ossText,ossBin;
					ossText << m_strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".txt";
					hypothesisLattice->store(ossText.str().c_str(),FILE_FORMAT_TEXT);
					ossBin << m_strFolderLattices << PATH_SEPARATOR << strUtteranceId << ".bin";
					hypothesisLattice->store(ossBin.str().c_str(),FILE_FORMAT_BINARY);
					delete hypothesisLattice;	
				}	
			}
			
			// best path
			BestPath *bestPathUtterance = m_wfsaDecoder->getBestPath();		
			if (bestPathUtterance != NULL) {
				dLikelihoodTotal += bestPathUtterance->getPathScore();
				bestPathUtterance->print();
				bestPathUtterance->write(fileHypothesis.getStream(),strUtteranceId.c_str(),
					bBestSinglePathOutputSentenceDelimiters,bBestSinglePathOutputFillers,bBestSinglePathOutputConfidenceValues);
				// state-alignment
				if (m_bOutputAlignment) {
					VPhoneAlignment *vPhoneAlignment = m_viterbi->align(fFeatures,iFeatures,bestPathUtterance);
					if (vPhoneAlignment != NULL) {
						AlignmentFile *alignmentFile = new AlignmentFile(m_phoneSet,m_lexiconManager);
						char strFileAlignment[1024+1];
						sprintf(strFileAlignment,"%s%c%s.ali",m_strFolderAlignments,PATH_SEPARATOR,strUtteranceId.c_str());
						alignmentFile->store(*vPhoneAlignment,strFileAlignment);
						delete alignmentFile;
						AlignmentFile::destroyPhoneAlignment(vPhoneAlignment);
					}	
				}
				delete bestPathUtterance;
			}
			
			// clean-up
			if (fFeatures != NULL) {
				delete [] fFeatures;	
			}
			
			//if (iUtterance == 10)
			//	break;
		}	
		fileHypothesis.close();
		delete batchFile;	
	} 
	
	double dTimeEnd = TimeUtils::getTimeMilliseconds();
	double dTimeSeconds = (dTimeEnd-dTimeBegin)/1000.0;
	double dRTF = dTimeSeconds/(((float)iFeatureVectorsTotal)/100.0); 
   
	printf("- summary ------------------------------------\n");
	printf("# utterances: %d speech time: %.2f seconds\n",iUtterances,((float)iFeatureVectorsTotal)/100.0);
	printf("decoding time: %.2f seconds (RTF: %5.2f)\n",dTimeSeconds,dRTF);
	printf("likelihood: %.4f (per frame: %8.4f)\n",dLikelihoodTotal,dLikelihoodTotal/((float)iFeatureVectorsTotal));
	printf("----------------------------------------------\n");	

	// clean up
	delete m_wfsaDecoder;
	delete m_lexiconManager;
	delete m_hmmManager;
	delete m_phoneSet;
	if (m_viterbi) {
		delete m_viterbi;
	}
	delete m_commandLineManager;
	delete m_wfsAcceptor;
	delete m_configuration;

	return 0;
}

