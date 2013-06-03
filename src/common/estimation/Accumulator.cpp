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


#include "Accumulator.h"
#include "BatchFile.h"
#include "Numeric.h"

namespace Bavieca {

// constructor (logical accumulator)
Accumulator::Accumulator(int iDim, int iCovarianceModeling, unsigned char *iIdentity, unsigned char iContextModelingOrder) {

	// accumulator type
	m_iType = ACCUMULATOR_TYPE_LOGICAL;

	m_iDim = iDim;
	m_iCovarianceModeling = iCovarianceModeling;
	
	// identity
	int iLength = iContextModelingOrder+3;
	m_iIdentity = new unsigned char[iLength];
	for(int i=0 ; iIdentity[i] != UCHAR_MAX ; ++i) {
		m_iIdentity[i] = iIdentity[i];
	}
	m_iIdentity[iLength-1] = UCHAR_MAX;
	m_iContextModelingOrder = iContextModelingOrder;
	m_iContextSize = (iContextModelingOrder-1)/2;
	
	// statistics
	m_vObservation = new Vector<double>(m_iDim);
	if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_vObservationSquare = new Vector<double>(m_iDim);
		m_mObservationSquare = NULL;
	} else {
		m_vObservationSquare = NULL;
		m_mObservationSquare = new SMatrix<double>(m_iDim);
	}
	m_dOccupation = 0.0;
	m_bDataAllocated = true;	
	// next 
	m_accumulatorNext = NULL;
	
	// unused fields
	m_iHMMState = -1;
	m_iGaussianComponent = -1;	
}

// constructor (physical accumulator)
Accumulator::Accumulator(int iDim, int iCovarianceModeling, int iHMMState, int iGaussianComponent) {

	// accumulator type
	m_iType = ACCUMULATOR_TYPE_PHYSICAL;
	
	m_iDim = iDim;
	m_iCovarianceModeling = iCovarianceModeling;
	
	// identity
	m_iHMMState = iHMMState;
	m_iGaussianComponent = iGaussianComponent;
	
	// data
	m_vObservation = new Vector<double>(m_iDim);
	if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_vObservationSquare = new Vector<double>(m_iDim);
		m_mObservationSquare = NULL;
	} else {
		m_vObservationSquare = NULL;
		m_mObservationSquare = new SMatrix<double>(m_iDim);
	}
	m_dOccupation = 0.0;
	m_bDataAllocated = true;		
	// next 
	m_accumulatorNext = NULL;
	
	// unused fields
	m_iIdentity	= NULL;
	m_iContextSize = UCHAR_MAX;
	m_iContextModelingOrder = UCHAR_MAX;
}

// copy constructor
Accumulator::Accumulator(Accumulator *accumulator) {

	m_iType = accumulator->m_iType;
	m_iDim = accumulator->m_iDim;
	m_iCovarianceModeling = accumulator->m_iCovarianceModeling;
		
	// logical accumulator
	if (m_iType == ACCUMULATOR_TYPE_LOGICAL) {
	
		assert(accumulator->m_iIdentity != NULL);
		// create a copy of the identity
		m_iIdentity = getCopyIdentity(accumulator->m_iIdentity);	
		m_iContextModelingOrder = accumulator->m_iContextModelingOrder;
		m_iContextSize = accumulator->m_iContextSize;
		assert(accumulator->m_accumulatorNext == NULL);
		m_accumulatorNext = NULL;
	
		// unused fields
		m_iHMMState = -1;
		m_iGaussianComponent = -1;
	} 
	// physical accumulator
	else {
		assert(m_iType == ACCUMULATOR_TYPE_PHYSICAL);
		m_iHMMState = accumulator->m_iHMMState;
		m_iGaussianComponent = accumulator->m_iGaussianComponent;
		
		// unused fields
		m_iContextModelingOrder = UCHAR_MAX;
		m_iContextSize = UCHAR_MAX;
		m_iIdentity = NULL;
		m_accumulatorNext = NULL;
	}
	
	// statistics
	m_vObservation = new Vector<double>(accumulator->getObservation());
	if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_vObservationSquare = new Vector<double>(accumulator->getObservationSquareDiag());
		m_mObservationSquare = NULL;
	} else {
		m_vObservationSquare = NULL;
		m_mObservationSquare = new SMatrix<double>(accumulator->getObservationSquareFull());
	}
	m_dOccupation = accumulator->m_dOccupation;			
	m_bDataAllocated = true;	
}

// constructor
Accumulator::Accumulator() {

	m_bDataAllocated = false;
}

// destructor
Accumulator::~Accumulator() {

	if (m_bDataAllocated) {
		delete m_vObservation;
		if (m_vObservationSquare) {
			delete m_vObservationSquare;
		} else {
			delete m_mObservationSquare;
		}
		m_vObservation = NULL;
		m_vObservationSquare = NULL;
		m_mObservationSquare = NULL;
	}
	if (m_iIdentity != NULL) {
		delete [] m_iIdentity;
		m_iIdentity = NULL;
		m_accumulatorNext = NULL;
	}
}

// store the accumulator
void Accumulator::store(FileOutput &file) {

	// logical accumulator
	if (m_iType == ACCUMULATOR_TYPE_LOGICAL) {

		// identity	
		for(int i = 0 ; m_iIdentity[i] != UCHAR_MAX ; ++i) {
			IOBase::write(file.getStream(),m_iIdentity[i]);
		}
		unsigned char iEnd = UCHAR_MAX;
		IOBase::write(file.getStream(),iEnd);
	}
	// physical accumulator
	else {
		assert(m_iType == ACCUMULATOR_TYPE_PHYSICAL);
	
		IOBase::write(file.getStream(),m_iHMMState);
		IOBase::write(file.getStream(),m_iGaussianComponent);
	}
	
	// accumulator data
	m_vObservation->writeData(file.getStream());
	if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_vObservationSquare->writeData(file.getStream());
	} else {
		m_mObservationSquare->writeData(file.getStream());	
	}	
	IOBase::write(file.getStream(),m_dOccupation);
}

// load an accumulator from disk
Accumulator *Accumulator::load(FileInput &file, int iDim, int iCovarianceModeling, unsigned char iType,
	unsigned char iContextModelingOrder) {

	Accumulator *accumulator = new Accumulator();
	accumulator->m_iDim = iDim;
	accumulator->m_iCovarianceModeling = iCovarianceModeling;
	accumulator->m_iType = iType;
	accumulator->m_iContextModelingOrder = iContextModelingOrder;
	
	// logical accumulator
	if (iType == ACCUMULATOR_TYPE_LOGICAL) {	

		accumulator->m_iContextSize = (iContextModelingOrder-1)/2;
		accumulator->m_accumulatorNext = NULL;
	
		// identity	
		accumulator->m_iIdentity = new unsigned char[iContextModelingOrder+3];
		for(int i=0 ; ; ++i) {
			IOBase::read(file.getStream(),&accumulator->m_iIdentity[i]);
			if (accumulator->m_iIdentity[i] == UCHAR_MAX) {
				assert(i == iContextModelingOrder+2);
				break;
			}
		}
		
		// unused fields
		accumulator->m_iHMMState = -1;
		accumulator->m_iGaussianComponent = -1;	
	}
	// logical accumulator
	else {
		assert(iType == ACCUMULATOR_TYPE_PHYSICAL);	
		
		IOBase::read(file.getStream(),&accumulator->m_iHMMState);
		IOBase::read(file.getStream(),&accumulator->m_iGaussianComponent);
		
		// unused fields
		accumulator->m_iIdentity = NULL;
		accumulator->m_iContextSize = UCHAR_MAX;
		accumulator->m_accumulatorNext = NULL;
	}
	
	// read data
	accumulator->m_bDataAllocated = true;
	accumulator->m_vObservation = new Vector<double>(iDim);
	accumulator->m_vObservation->readData(file.getStream());
	if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		accumulator->m_vObservationSquare = new Vector<double>(iDim);
		accumulator->m_vObservationSquare->readData(file.getStream());
		accumulator->m_mObservationSquare = NULL;
	} else {
		accumulator->m_vObservationSquare = NULL;
		accumulator->m_mObservationSquare = new SMatrix<double>(iDim);
		accumulator->m_mObservationSquare->readData(file.getStream());
	}
	IOBase::read(file.getStream(),&accumulator->m_dOccupation);		
	
	return accumulator;
}

// return the context modeling order in numeric format
unsigned char Accumulator::getContextModelingOrder(const char *strContextModelingOrder) {

	if (strcmp(strContextModelingOrder,HMM_CONTEXT_MODELING_MONOPHONES_STR) == 0) {
		return HMM_CONTEXT_MODELING_MONOPHONES;
	} 
	else if (strcmp(strContextModelingOrder,HMM_CONTEXT_MODELING_TRIPHONES_STR) == 0) {
		return HMM_CONTEXT_MODELING_TRIPHONES;
	} 
	else if (strcmp(strContextModelingOrder,HMM_CONTEXT_MODELING_PENTAPHONES_STR) == 0) {
		return HMM_CONTEXT_MODELING_PENTAPHONES;
	} 
	else if (strcmp(strContextModelingOrder,HMM_CONTEXT_MODELING_HEPTAPHONES_STR) == 0) {
		return HMM_CONTEXT_MODELING_HEPTAPHONES;
	} 
	else if (strcmp(strContextModelingOrder,HMM_CONTEXT_MODELING_NONAPHONES_STR) == 0) {
		return HMM_CONTEXT_MODELING_NONAPHONES;
	} 
	else if (strcmp(strContextModelingOrder,HMM_CONTEXT_MODELING_ENDECAPHONES_STR) == 0) {
		return HMM_CONTEXT_MODELING_ENDECAPHONES;
	} 
	else {
		return UCHAR_MAX; 
	}
}

// return the context modeling order in numeric format
const char *Accumulator::getContextModelingOrder(unsigned char iContextModelingOrder) {

	switch(iContextModelingOrder) {
		case HMM_CONTEXT_MODELING_MONOPHONES: {
			return HMM_CONTEXT_MODELING_MONOPHONES_STR;
		}
		case HMM_CONTEXT_MODELING_TRIPHONES: {
			return HMM_CONTEXT_MODELING_TRIPHONES_STR;
		}
		case HMM_CONTEXT_MODELING_PENTAPHONES: {
			return HMM_CONTEXT_MODELING_PENTAPHONES_STR;
		}
		case HMM_CONTEXT_MODELING_HEPTAPHONES: {
			return HMM_CONTEXT_MODELING_HEPTAPHONES_STR;
		}
		case HMM_CONTEXT_MODELING_NONAPHONES: {
			return HMM_CONTEXT_MODELING_NONAPHONES_STR;
		}
		case HMM_CONTEXT_MODELING_ENDECAPHONES: {
			return HMM_CONTEXT_MODELING_ENDECAPHONES_STR;
		}
		default: {
			return NULL;
		}
	}
}

// dump logical accumulators to disk
void Accumulator::storeAccumulators(const char *strFile, int iDim, 
	int iCovarianceModeling, unsigned char iContextModelingOrderWW, 
	unsigned char iContextModelingOrderCW, MAccumulatorLogical &mAccumulator) {

	FileOutput file(strFile,true);
	file.open();
	
	// logical accumulators: one accumulator for each context dependent unit
	
	// write the accumulator type
	unsigned char iType = ACCUMULATOR_TYPE_LOGICAL;
	IOBase::write(file.getStream(),iType);
	
	// write the number of accumulators
	int iAccumulators = 0;
	for(MAccumulatorLogical::iterator it = mAccumulator.begin() ; it != mAccumulator.end() ; ++it) {
		assert(it->second->m_iType == ACCUMULATOR_TYPE_LOGICAL);
		// skip if no occupation (this is possible since there can be optional symbols in the hmm-graph
		// used for the forward-backward alignment)
		if (it->second->getOccupation() != 0.0) {
			++iAccumulators;
		}	
	}
	IOBase::write(file.getStream(),iAccumulators);
	IOBase::write(file.getStream(),iDim);
	IOBase::write(file.getStream(),iCovarianceModeling);
	IOBase::write(file.getStream(),iContextModelingOrderWW);
	IOBase::write(file.getStream(),iContextModelingOrderCW);

	// write the accumulators to disk one by one
	for(MAccumulatorLogical::iterator it = mAccumulator.begin() ; it != mAccumulator.end() ; ++it) {
		assert(it->second->m_iType == ACCUMULATOR_TYPE_LOGICAL);
		// skip if no occupation (this is possible since there can be optional symbols in the hmm-graph
		// used for the  alignment)
		if (it->second->getOccupation() == 0.0) {
			continue;
		}	
		it->second->store(file);
	}	
	
	file.close();
}

// dump physical accumulators to disk
void Accumulator::storeAccumulators(const char *strFile, int iDim, 
	int iCovarianceModeling, int iHMMStates, int iGaussianComponents, MAccumulatorPhysical &mPhysicalAccumulator) {

	FileOutput file(strFile,true);
	file.open();	
	
	// write the accumulator type
	unsigned char iType = ACCUMULATOR_TYPE_PHYSICAL;
	IOBase::write(file.getStream(),iType);
	
	int iAccumulators = (int)mPhysicalAccumulator.size();
	IOBase::write(file.getStream(),iAccumulators);
	IOBase::write(file.getStream(),iDim);
	IOBase::write(file.getStream(),iCovarianceModeling);	
	IOBase::write(file.getStream(),iHMMStates);	
	IOBase::write(file.getStream(),iGaussianComponents);	

	// write the accumulators to disk one by one
	for(MAccumulatorPhysical::iterator it = mPhysicalAccumulator.begin() ; it != mPhysicalAccumulator.end() ; ++it) {
		assert(it->second->m_iType == ACCUMULATOR_TYPE_PHYSICAL);
		it->second->store(file);
	}	
	
	file.close();
}

// add accumulators
void Accumulator::add(Accumulator *accumulator) {

	m_vObservation->add(accumulator->getObservation());
	if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
		m_vObservationSquare->add(accumulator->getObservationSquareDiag());
	} else {
		m_mObservationSquare->add(accumulator->getObservationSquareFull());
	}
	m_dOccupation += accumulator->m_dOccupation;	
}

// add logical accumulators
void Accumulator::addAccumulators(MAccumulatorLogical &mAccumulator1, MAccumulatorLogical &mAccumulator2) {
	
	for(MAccumulatorLogical::iterator it = mAccumulator2.begin() ; it != mAccumulator2.end() ; ++it) {
		MAccumulatorLogical::iterator jt = mAccumulator1.find(it->second->m_iIdentity);
		if (jt != mAccumulator1.end()) {
			jt->second->add(it->second);
		} else {
			Accumulator *accumulator = new Accumulator(it->second);	
			mAccumulator1.insert(MAccumulatorLogical::value_type(accumulator->m_iIdentity,accumulator));
		}
	}
}

// add physical accumulators
void Accumulator::addAccumulators(MAccumulatorPhysical &mAccumulator1, MAccumulatorPhysical &mAccumulator2) {
	
	for(MAccumulatorPhysical::iterator it = mAccumulator2.begin() ; it != mAccumulator2.end() ; ++it) {
		unsigned int iKey = getPhysicalAccumulatorKey(it->second->m_iHMMState,it->second->m_iGaussianComponent);
		MAccumulatorPhysical::iterator jt = mAccumulator1.find(iKey);
		if (jt != mAccumulator1.end()) {
			jt->second->add(it->second);
		} else {
			Accumulator *accumulator = new Accumulator(it->second);
			mAccumulator1.insert(MAccumulatorPhysical::value_type(iKey,accumulator));
		}
	}
}

// load accumulators from a file
void Accumulator::loadAccumulators(const char *strFile, MAccumulatorPhysical &mAccumulatorPhysical, AccMetadata &metadata) {

	FileInput file(strFile,true);
	file.open();
	
	int iAccumulators;
	unsigned char iType;
	
	IOBase::read(file.getStream(),&iType);
	if (iType != ACCUMULATOR_TYPE_PHYSICAL) {
		BVC_ERROR << "physical accumulators expected on file: " << strFile;
	} 
	
	// get the number of accumulators and their properties
	IOBase::read(file.getStream(),&iAccumulators);
	if (iAccumulators < 1) {
		BVC_ERROR << "no accumulators found in file: " << strFile;	
	}
	IOBase::read(file.getStream(),&metadata.iDim);
	IOBase::read(file.getStream(),&metadata.iCovarianceModeling);	
	IOBase::read(file.getStream(),&metadata.iHMMStates);	
	IOBase::read(file.getStream(),&metadata.iGaussianComponents);
	metadata.iContextModelingOrderWW = metadata.iContextModelingOrderCW = UCHAR_MAX;
	
	// read the accumulators from disk one by one
	for(int i=0 ; i<iAccumulators ; ++i) {
		Accumulator *accumulator = Accumulator::load(file,metadata.iDim,metadata.iCovarianceModeling,
			ACCUMULATOR_TYPE_PHYSICAL,UCHAR_MAX);
		assert(accumulator);
		assert(accumulator->m_iType == ACCUMULATOR_TYPE_PHYSICAL);
		mAccumulatorPhysical.insert(MAccumulatorPhysical::value_type(
			Accumulator::getPhysicalAccumulatorKey(accumulator->getHMMState(),
			accumulator->getGaussianComponent()),accumulator));
	}	

	file.close();
}

// load accumulators from file
void Accumulator::loadAccumulators(const char *strFile, MAccumulatorLogical &mAccumulatorLogical, AccMetadata &metadata) {

	FileInput file(strFile,true);
	file.open();
	
	int iAccumulators;
	unsigned char iType;
	
	IOBase::read(file.getStream(),&iType);
	if (iType != ACCUMULATOR_TYPE_LOGICAL) {
		BVC_ERROR << "logical accumulators expected, physical accumulators found";
	}
	
	// read the number of accumulators
	IOBase::read(file.getStream(),&iAccumulators);
	if (iAccumulators < 1) {
		BVC_ERROR << "no accumulators found";
	}
	
	IOBase::read(file.getStream(),&metadata.iDim);
	IOBase::read(file.getStream(),&metadata.iCovarianceModeling);	
	IOBase::read(file.getStream(),&metadata.iContextModelingOrderWW);	
	IOBase::read(file.getStream(),&metadata.iContextModelingOrderCW);
	metadata.iHMMStates = metadata.iGaussianComponents = -1;	
	
	// read the accumulators from disk one by one
	for(int i=0 ; i < iAccumulators ; ++i) {
		Accumulator *accumulator = Accumulator::load(file,metadata.iDim,metadata.iCovarianceModeling,
			ACCUMULATOR_TYPE_LOGICAL,metadata.iContextModelingOrderWW);
		assert(accumulator);
		assert(accumulator->m_iType == ACCUMULATOR_TYPE_LOGICAL);
		mAccumulatorLogical.insert(MAccumulatorLogical::value_type(accumulator->getIdentity(),accumulator));
	}	
		
	file.close();
}

// load and combine physical accumulators from multiple files
void Accumulator::loadAccumulatorList(const char *strFileList, MAccumulatorPhysical &mAccumulatorPhysical, AccMetadata &metadata) {

	BatchFile batchFile(strFileList,"acc");
	batchFile.load();

	assert(batchFile.size() > 0);
	for(unsigned int i=0 ; i < batchFile.size() ; ++i) {
		
		MAccumulatorPhysical mAccumulatorPhysicalAux;
		
		// load the accumulators
		loadAccumulators(batchFile.getField(i,0u),mAccumulatorPhysicalAux,metadata);
		// combine the accumulators
		addAccumulators(mAccumulatorPhysical,mAccumulatorPhysicalAux);
		destroy(mAccumulatorPhysicalAux);
	}
}

// load and combine logical accumulators from multiple files
void Accumulator::loadAccumulatorList(const char *strFileList, MAccumulatorLogical &mAccumulatorLogical, AccMetadata &metadata) {

	BatchFile batchFile(strFileList,"acc");
	batchFile.load();

	for(unsigned int i=0 ; i < batchFile.size() ; ++i) {
		
		MAccumulatorLogical mAccumulatorLogicalAux;
		
		// load the accumulators
		loadAccumulators(batchFile.getField(i,0u),mAccumulatorLogicalAux,metadata);
		// combine the accumulators
		addAccumulators(mAccumulatorLogical,mAccumulatorLogicalAux);
		destroy(mAccumulatorLogicalAux);	
	}
}

// destroy the accumulators
void Accumulator::destroy(MAccumulatorLogical &mAccumulatorLogical) {

	VAccumulator vAccumulator;

	// accumulators cannot be directly removed from the map without corrupting it
	for(MAccumulatorLogical::iterator it = mAccumulatorLogical.begin() ; it != mAccumulatorLogical.end() ; ++it) {
		vAccumulator.push_back(it->second);
	}
	mAccumulatorLogical.clear();
	
	for(VAccumulator::iterator it = vAccumulator.begin() ; it != vAccumulator.end() ; ++it) {
		delete *it;	
	}
	vAccumulator.clear();
}

// destroy the accumulators
void Accumulator::destroy(MAccumulatorPhysical &mAccumulatorPhysical) {
 
	// accumulators cannot be directly removed from the map without corrupting it
	for(MAccumulatorPhysical::iterator it = mAccumulatorPhysical.begin() ; it != mAccumulatorPhysical.end() ; ++it) {
		delete it->second;
	}
	mAccumulatorPhysical.clear();
}

// print accumulator info
void Accumulator::print(MAccumulatorPhysical &mAccumulatorPhysical) {

	double dOccupation = 0.0;

	// accumulators cannot be directly removed from the map without corrupting it
	for(MAccumulatorPhysical::iterator it = mAccumulatorPhysical.begin() ; it != mAccumulatorPhysical.end() ; ++it) {
		dOccupation += it->second->m_dOccupation;
	}
	
	cout << "# accumulators:   " << mAccumulatorPhysical.size() << endl;
	cout << "total occupation: " << dOccupation << endl;
}

// print accumulator info
void Accumulator::print(MAccumulatorLogical &mAccumulatorLogical) {

	double dOccupation = 0.0;

	// accumulators cannot be directly removed from the map without corrupting it
	for(MAccumulatorLogical::iterator it = mAccumulatorLogical.begin() ; it != mAccumulatorLogical.end() ; ++it) {
		dOccupation += it->second->m_dOccupation;
	}
	
	cout << "# accumulators:   " << mAccumulatorLogical.size() << endl;
	cout << "total occupation: " << dOccupation << endl;
}

// adapt the accumulators to the given within-word and cross-word context length
void Accumulator::adaptContextWidth(MAccumulatorLogical &mAccumulator, unsigned char iContextSize, unsigned char iContextSizeNew) {

	assert(iContextSize < iContextSizeNew);
	
	// (1) shorten accumulators because of within-word context
	int iAccumulatorsOriginal = (int)mAccumulator.size();
	MAccumulatorLogical mAccumulatorTemp;
	if (iContextSizeNew < iContextSize) {
		for(MAccumulatorLogical::iterator it = mAccumulator.begin() ; it != mAccumulator.end() ; ++it) {
			Accumulator *accumulatorSrc = it->second;
			// shorten the identity
			unsigned char *iIdentityShort = getShorterIdentity(accumulatorSrc->getIdentity(),iContextSize,iContextSizeNew);
			MAccumulatorLogical::iterator jt = mAccumulatorTemp.find(iIdentityShort);
			// create a new entry
			if (jt == mAccumulatorTemp.end()) {
				Accumulator *accumulatorDst = new Accumulator(accumulatorSrc);
				accumulatorDst->shortenContext(iContextSizeNew);
				mAccumulatorTemp.insert(MAccumulatorLogical::value_type(accumulatorDst->getIdentity(),accumulatorDst));	
			} 
			// add the data to the already existing accumulator
			else {
				jt->second->add(accumulatorSrc);	
			}
			delete [] iIdentityShort;
			delete it->second;
		}
		mAccumulator.swap(mAccumulatorTemp);
	}
	int iAccumulatorsCompacted = (int)mAccumulator.size();
	
	cout << "original: " << iAccumulatorsOriginal << endl;
	cout << "original: " << iAccumulatorsCompacted << endl;
}

};	// end-of-namespace



