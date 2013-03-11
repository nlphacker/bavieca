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


#include "BatchFile.h"
#include "CommandLineManager.h"
#include "ConfigurationFeatures.h"
#include "FeatureFile.h"
#include "Global.h"
#include "TimeUtils.h"
#include "Transform.h"

using namespace std;

#include <string>
#include <map>

using namespace Bavieca;

float *applyTransform(Transform *transform, float *fFeatures, int iFeatures, int iDim, int &iDimX);

// main for the feature transformation tool: "paramX"
int main(int argc, char *argv[]) {

	try {

		// (1) define command line parameters
		CommandLineManager commandLineManager("paramx",SYSTEM_VERSION,SYSTEM_AUTHOR,SYSTEM_DATE);
		commandLineManager.defineParameter("-cfg","feature configuration",PARAMETER_TYPE_FILE,false);	
		commandLineManager.defineParameter("-tra","feature transformation",PARAMETER_TYPE_FILE,false);
		commandLineManager.defineParameter("-bat","batch file containing pairs (feaIn feaOut)",
			PARAMETER_TYPE_FILE,true);	
		commandLineManager.defineParameter("-in","input feature vectors",PARAMETER_TYPE_FILE,true);
		commandLineManager.defineParameter("-out","output feature vectors",PARAMETER_TYPE_FILE,true);
		
		// parse the parameters
		if (commandLineManager.parseParameters(argc,argv) == false) {
			return -1;
		}
		
		// load the parameters
		const char *strFileConfiguration = commandLineManager.getParameterValue("-cfg");
		const char *strFileTransform = commandLineManager.getParameterValue("-tra");
		const char *strFileBatch = commandLineManager.getParameterValue("-bat");	
		const char *strFileInput = commandLineManager.getParameterValue("-in");
		const char *strFileOutput = commandLineManager.getParameterValue("-out");
		
		// load the feature configuration
		ConfigurationFeatures configurationFeatures(strFileConfiguration);
		configurationFeatures.load();
			
		// load the transformation
		Transform transform;
		transform.load(strFileTransform);
		//transform->print(true);	
		
		// single feature file
		if (strFileBatch == NULL) {
		
			// load the features
			int iDimensionality = atoi(configurationFeatures.getStrParameterValue("cepstralCoefficients"))+1;
			FeatureFile featureFile(strFileInput,MODE_READ,FORMAT_FEATURES_FILE_DEFAULT,iDimensionality);
			featureFile.load();
			
			unsigned int iFeatures = 0;
			float *fFeatures = featureFile.getFeatureVectors(&iFeatures);	
		
			// apply the transformation		
			int iDimensionalityX = -1;
			float *fFeaturesX = applyTransform(&transform,fFeatures,iFeatures,iDimensionality,iDimensionalityX);
			if (!fFeaturesX) {
				BVC_ERROR << "unable to apply the transform";
			}
			
			// create the transformed feature file
			FeatureFile featureFileX(strFileOutput,MODE_WRITE,FORMAT_FEATURES_FILE_DEFAULT,iDimensionalityX);
			featureFileX.store(fFeaturesX,iFeatures);
			
			delete [] fFeatures;
			delete [] fFeaturesX;
		} 
		// batch mode
		else {	
		
			BatchFile batchFile(strFileBatch,"featuresIn|featuresOut");
			batchFile.load();
			for(unsigned int i=0 ; i < batchFile.size() ; ++i) {
				
				// load the features
				int iDimensionality = atoi(configurationFeatures.getStrParameterValue("cepstralCoefficients"))+1;
				FeatureFile featureFile(batchFile.getField(i,"featuresIn"),MODE_READ,
					FORMAT_FEATURES_FILE_DEFAULT,iDimensionality);
				featureFile.load();
				
				unsigned int iFeatures = 0;
				float *fFeatures = featureFile.getFeatureVectors(&iFeatures);	
			
				// apply the transformation		
				int iDimensionalityX = -1;
				float *fFeaturesX = applyTransform(&transform,fFeatures,iFeatures,iDimensionality,iDimensionalityX);
				if (!fFeaturesX) {
					BVC_ERROR << "unable to apply the transform";
				}
				
				// create the transformed feature file
				FeatureFile featureFileX(batchFile.getField(i,"featuresOut"),MODE_WRITE,
					FORMAT_FEATURES_FILE_DEFAULT,iDimensionalityX);
				featureFileX.store(fFeaturesX,iFeatures);
				
				delete [] fFeatures;
				delete [] fFeaturesX;
			}	
		}
	
	} catch (ExceptionBase &e) {
	
		std::cerr << e.what() << std::endl;
		return -1;
	}	
			
	return 0;	
}

// apply the transform to the features
// case A: columns in transform = feature dimension (linear transform)
// case B: columns in transform = feature dimension+1 (affine transform, append 1 to the feature vector)
float *applyTransform(Transform *transform, float *fFeatures, int iFeatures, int iDim, int &iDimX) {

	float *fFeaturesX = NULL;
	int iRows = -1;
	int iColumns = -1;
	float *fA = transform->getTransform().getData();
	iDimX = iRows;
	
	// linear transform (for example HLDA)
	if (iColumns == iDim) {
		
		fFeaturesX = new float[iFeatures*iDimX];		
		
		// transform the features
		for(int f=0 ; f<iFeatures ; ++f) {
			for(int i=0 ; i<iDimX ; ++i) {
				fFeaturesX[f*iDimX+i] = 0.0f;
				for(int j=0 ; j<iDim ; ++j) {
					fFeaturesX[f*iDimX+i] += fA[i*iDim+j]*fFeatures[f*iDim+j];
				}
			}
		}
	}
	// affine transform (for example fMLLR)
	else if (iColumns == iDim+1) {
	
		for(int i=0 ; i < iFeatures ; ++i) {
			float *fFea = &fFeatures[i*iDim];
			for(int j=0 ; j < iDim ; ++j) {
				fFeaturesX[j] = fA[j*(iDim+1)];
				for(int k=1 ; k < iDim+1 ; ++k) {
					fFeaturesX[j] += fA[j*(iDim+1)+k]*fFea[k-1];
				}	
			}	
		}
	}
	else {
		return NULL;
	}

	return fFeaturesX;
}


