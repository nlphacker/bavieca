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

#ifndef ACCUMULATORFILE_H
#define ACCUMULATORFILE_H

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "Global.h"
#include "Alignment.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "IOBase.h"

#define MODE_READ		0
#define MODE_WRITE	1

// accumulator type
#define ACCUMULATOR_TYPE_LOGICAL					0
#define ACCUMULATOR_TYPE_PHYSICAL				1

// context modeling order (numeric format)
#define HMM_CONTEXT_MODELING_MONOPHONES		1	
#define HMM_CONTEXT_MODELING_TRIPHONES			3
#define HMM_CONTEXT_MODELING_PENTAPHONES		5
#define HMM_CONTEXT_MODELING_HEPTAPHONES		7
#define HMM_CONTEXT_MODELING_NONAPHONES		9
#define HMM_CONTEXT_MODELING_ENDECAPHONES		11

// context modeling order (string format)
#define HMM_CONTEXT_MODELING_MONOPHONES_STR		"monophones"
#define HMM_CONTEXT_MODELING_TRIPHONES_STR		"triphones"
#define HMM_CONTEXT_MODELING_PENTAPHONES_STR		"pentaphones"
#define HMM_CONTEXT_MODELING_HEPTAPHONES_STR		"heptaphones"
#define HMM_CONTEXT_MODELING_NONAPHONES_STR		"nonaphones"
#define HMM_CONTEXT_MODELING_ENDECAPHONES_STR	"endecaphones"

#define HMM_CONTEXT_SIZE_MAX						5       // endecaphones have 5 on each side

// context type attending to word location
#define HMM_CONTEXT_TYPE_WITHIN_WORD			0
#define HMM_CONTEXT_TYPE_CROSS_WORD				1

#define MAX_IDENTITY_LENGTH		32

// covariance modeling type
#define COVARIANCE_MODELLING_TYPE_DEFAULT			-1
#define COVARIANCE_MODELLING_TYPE_DIAGONAL		0
#define COVARIANCE_MODELLING_TYPE_FULL				1

using namespace std;

#include <string>
#include <vector>

#ifdef __linux__
using namespace __gnu_cxx;
#include <ext/hash_map>
#endif

#ifdef _WIN32
#include <hash_map>
#endif

#include "Allocator.h"

class Accumulator;

typedef vector<Accumulator*> VAccumulator;

typedef struct {
	int iFeatureDimensionality;
	int iCovarianceModeling;
	int iHMMStates;								// physical
	int iGaussianComponents;					//	physical
	int iContextModelingOrderWW;				// logical acc
	int iContextModelingOrderCW;				// logical acc
} AccMetadata;

// ad-hoc functions to use Duple as the key in a hash_map data structure
struct MAccumulatorFunctions
{

#ifdef _WIN32
	static const size_t bucket_size = 4;
	static const size_t min_buckets = 8;
#endif

	// comparison function
	bool operator()(const unsigned char *contextUnit1, const unsigned char *contextUnit2) const {
	
		for(int i=0 ; contextUnit1[i] != UCHAR_MAX ; ++i) {
			assert(contextUnit2[i] != UCHAR_MAX);
			if (contextUnit1[i] != contextUnit2[i]) {
				return false;
			}
		}
		
		return true;
	}
	
	// hash function
	size_t operator()(const unsigned char *contextUnit) const {
	
		unsigned int iAcc = 0;
		unsigned int iAux = 0;
		for(int i=0 ; contextUnit[i] != UCHAR_MAX ; ++i) {
			//printf("%u ",contextUnit[i]);
			if (i <= 3) {	
				iAcc <<= (8*i);
				iAcc += contextUnit[i];
			} else {
				iAux = contextUnit[i];
				iAux <<= (8*(i%4));
				iAcc ^= iAux;
			}
		}
		//printf("-> %u\n",iAcc);
	
		return iAcc;
	}
};

// structure for easy access to logical accumulators (maps context dependent phones to accumulators)
#ifdef __linux__
typedef hash_map<unsigned char*,Accumulator*,MAccumulatorFunctions,MAccumulatorFunctions> MAccumulatorLogical;
#elif _WIN32
typedef hash_map<unsigned char*,Accumulator*,MAccumulatorFunctions> MAccumulatorLogical;
#endif

// structure to keep physical accumulators
typedef hash_map<unsigned int,Accumulator*> MAccumulatorPhysical;

/**
	@author root <root@localhost.localdomain>
*/
class Accumulator {

	private:
	
		unsigned char m_iType;							// accumulator type (physical or logical)
		int m_iFeatureDimensionality;					// feature dimensionality
		int m_iCovarianceModeling;						// covariance modeling type
		bool m_bDataAllocated;							// data allocated	
		
		// only for logical accumulators
		unsigned char *m_iIdentity;						// accumulator identity
		unsigned char m_iContextModelingOrder;			// context modeling order
		unsigned char m_iContextSize;						// context size
		Accumulator *m_accumulatorNext;					// next accumulator in the linked list of accumulators
		
		// only for physical accumulators
		int m_iHMMState;									// HMM-state
		int m_iGaussianComponent;						// Gaussian-component
		
		// data
		double *m_dObservation;				// statistics needed to compute the mean 
		double *m_dObservationSquare;		//	statistics needed to compute the covariance
		double m_dOccupation;				// gaussian occupation
		
		// return the number of relevant elements in the covariance matrix
		inline int getCovarianceElements() {
		
			return getCovarianceElements(m_iFeatureDimensionality,m_iCovarianceModeling);
		}

		// return the number of relevant elements in the covariance matrix
		static inline int getCovarianceElements(int iFeatureDimensionality, int iCovarianceModeling) {
		
			if (iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				return iFeatureDimensionality;
			} else {
				assert(iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
				return (iFeatureDimensionality*(iFeatureDimensionality+1))/2;
			}
		}

	public:

		// constructor (logical accumulator)
		Accumulator(int iFeatureDimensionality, int iCovarianceModeling, unsigned char *iIdentity, unsigned char iContextModelingOrder);
		
		// constructor (physical accumulator) (used when loading accumulators from disk)
		Accumulator(int iFeatureDimensionality, int iCovarianceModeling, int iHMMState, int iGaussianComponent);
		
		// constructor (physical accumulator) (used when building an accumulator from a Gaussian component in HMMState)
		Accumulator(int iFeatureDimensionality, int iCovarianceModeling, int iHMMState, int iGaussianComponent, double *dObservation, double *dObservationSquare, double dOccupation);
		
		// copy constructor
		Accumulator(Accumulator *accumulator);	
		
		// constructor
		Accumulator();
		
		// destructor
		~Accumulator();
		
		// accumulate an observation
		inline void accumulateObservation(float *fFeatureVector, double dOccupation) {
		
			for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
				m_dObservation[i] += dOccupation*fFeatureVector[i];
			}
			// diagonal
			if (m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_DIAGONAL) {
				for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
					m_dObservationSquare[i] += dOccupation*(fFeatureVector[i]*fFeatureVector[i]);
				}
			}
			// full
			else {
				assert(m_iCovarianceModeling == COVARIANCE_MODELLING_TYPE_FULL);
				int iIndex = 0;
				for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
					for(int j=i ; j < m_iFeatureDimensionality ; ++j, ++iIndex) {
						m_dObservationSquare[iIndex] += dOccupation*(fFeatureVector[i]*fFeatureVector[j]);
					}
				}	
			}
			m_dOccupation += dOccupation;
			
			return;
		}
		
		// return the HMM-state and Gaussian component given a physical accumulator key
		static inline void getPhysicalAccumulatorValues(unsigned int iKey, int &iHMMState, int &iGaussianComponent) {
		
			iGaussianComponent = iKey%65536;
			iHMMState = iKey/65536;
			
			return;
		}		
		
		// build a physical accumulator key from the HMM-state and Gaussian-component number
		static inline unsigned int getPhysicalAccumulatorKey(int iHMMState, int iGaussianComponent) {
			
			assert(sizeof(unsigned int) >= 4);
			assert((iHMMState >= 0) && (iHMMState < 65536));
			unsigned int iKey = iHMMState*65536+iGaussianComponent;
	
			return iKey;
		}
	
		// build the identity
		static inline void buildIdentity(unsigned char *iIndentity, unsigned char *iPhoneLeft, unsigned char iPhone, 
			unsigned char *iPhoneRight, unsigned char iPosition, unsigned char iState, unsigned char iContextModelingOrder) {
			
			unsigned char iContextSize = (iContextModelingOrder-1)/2;
			for(int i=0 ; i < iContextSize ; ++i) {
				iIndentity[i] = iPhoneLeft[i];
				iIndentity[iContextSize+1+i] = iPhoneRight[i];
			}
			iIndentity[iContextSize] = iPhone;
			iIndentity[2*iContextSize+1] = iPosition;
			iIndentity[2*iContextSize+2] = iState;
			iIndentity[2*iContextSize+3] = UCHAR_MAX;
			
			return;
		}	
		
		// return a copy of the given identity
		static inline unsigned char *getCopyIdentity(unsigned char *iIdentity) {
		
			// count the elements in the given identity
			int i=0;
			while(iIdentity[i] != UCHAR_MAX) {
				i++;
			}
			assert(i > 0);
			unsigned char *iIdentityCopy = new unsigned char[i+1];
			int j;
			for(j = 0 ; j < i ; ++j) {
				iIdentityCopy[j] = iIdentity[j];
			}
			iIdentityCopy[j] = UCHAR_MAX;
		
			return iIdentityCopy;
		}
		
		// return the identity
		inline unsigned char *getIdentity() {
		
			return m_iIdentity;
		}
		
		// return the left context for the given position
		inline unsigned char getLeftPhone(unsigned char iPosition) {
		
			assert(iPosition < m_iContextSize);	
		
			return m_iIdentity[iPosition];
		}
		
		// return the right context for the given position
		inline unsigned char getRightPhone(unsigned char iPosition) {
		
			assert(iPosition < m_iContextSize);
		
			return m_iIdentity[m_iContextSize+iPosition+1];
		}
		
		// return the central phone
		inline unsigned char getPhone() {
		
			return m_iIdentity[m_iContextSize];
		}
		
		// return the within-word position
		inline unsigned char getPosition() {
		
			return m_iIdentity[m_iContextSize*2+1];
		}
		
		// return the HMM-state
		inline unsigned char getState() {
		
			return m_iIdentity[m_iContextSize*2+2];
		}
		
		inline unsigned char getContextModelingOrder() {
		
			return m_iContextModelingOrder;
		}
		
		inline double *getObservation() {
		
			return m_dObservation;
		}
		
		inline double *getObservationSquare() {
		
			return m_dObservationSquare;
		}
		
		inline double getOccupation() {
		
			return m_dOccupation;
		}
		
		inline void setNext(Accumulator *accumulator) {
		
			m_accumulatorNext = accumulator;
		}
		
		inline Accumulator *getNext() {
		
			return m_accumulatorNext;
		}
		
		inline int getHMMState() {
		
			return m_iHMMState;
		}
		
		inline int getGaussianComponent() {
		
			return m_iGaussianComponent;
		}
		
		// reset the accumulator
		inline void reset() {
		
			for(int i=0 ; i < m_iFeatureDimensionality ; ++i) {
				m_dObservation[i] = 0.0;
			}
			for(int i=0 ; i < getCovarianceElements() ; ++i) {
				m_dObservationSquare[i] = 0.0;
			}			
			m_dOccupation = 0.0;
			
			m_accumulatorNext = NULL;
				
			return;
		}
		
		// print the accumulator
		inline void print(PhoneSet *phoneSet) {
		
			unsigned char iPhoneContextPadding = phoneSet->getSize();
			
			if (m_iType == ACCUMULATOR_TYPE_LOGICAL) {
				for(int i=0 ; i < m_iContextSize ; ++i) {
					printf("%3s-",phoneSet->getStrPhone(getLeftPhone(i)));
				}
				printf("%3s(%d)[%d]",phoneSet->getStrPhone(getPhone()),getState(),getPosition());
				for(int i=0 ; i < m_iContextSize ; ++i) {
					printf("+%3s",phoneSet->getStrPhone(getRightPhone(i)));
				}	
				printf(" %12.4f\n",m_dOccupation);
			}
			
			return;
		}
		
		// return the context modeling order in numeric format
		static unsigned char getContextModelingOrder(const char *strContextModelingOrder);
		
		// return the context modeling order in numeric format
		static const char *getContextModelingOrder(unsigned char iContextModelingOrder);	
		
		// return the dimensionality
		inline int getDimensionality() {
		
			return m_iFeatureDimensionality;
		}
		
		// return the covariance modeling type
		inline int getCovarianceModeling() {
		
			return m_iCovarianceModeling;
		}
		
		// return whether the given context modeling order is valid
		static inline bool isValid(unsigned char iContextModelingOrder) {
		
			switch(iContextModelingOrder) {
				case HMM_CONTEXT_MODELING_MONOPHONES: 
				case HMM_CONTEXT_MODELING_TRIPHONES: 
				case HMM_CONTEXT_MODELING_PENTAPHONES: 
				case HMM_CONTEXT_MODELING_HEPTAPHONES: 
				case HMM_CONTEXT_MODELING_NONAPHONES: 
				case HMM_CONTEXT_MODELING_ENDECAPHONES: {		
					return true;
				}
				default : {
					return false;
				}
			}
		}	
		
		inline void shortenContext(unsigned char iContextSizeNew) {
		
			unsigned char *iIdentityNew = getShorterIdentity(m_iIdentity,m_iContextSize,iContextSizeNew);
			delete [] m_iIdentity;
			m_iIdentity = iIdentityNew;
			
			m_iContextSize = iContextSizeNew;
			m_iContextModelingOrder = iContextSizeNew*2+1;
		}
		
		// arithmetics -----------------------------------------------------------------------------------
		
		// add logical accumulators
		static void addAccumulators(MAccumulatorLogical &mAccumulator1, MAccumulatorLogical &mAccumulator2);
		
		// add physical accumulators
		static void addAccumulators(MAccumulatorPhysical &mPhysicalAccumulator1, MAccumulatorPhysical &mPhysicalAccumulator2);
		
		// add accumulators
		void add(Accumulator *accumulator);		
		
		// input/output from file ------------------------------------------------------------------------
		
		// store the accumulator to the given file
		bool store(FileOutput &file);
		
		// load the accumulator from the given file
		static Accumulator *load(FileInput &file, int iFeatureDimensionality, int iCovarianceModeling, unsigned char iType,
			unsigned char iContextModelingOrder);	
		
		// store logical accumulators to disk
		static void storeAccumulators(const char *strFile, int iFeatureDimensionality, int iCovarianceModeling,
			unsigned char iContextModelingOrderWW, unsigned char iContextModelingOrderCW, 
			MAccumulatorLogical &mAccumulatorLogical);
		
		// store physical accumulators to disk
		static void storeAccumulators(const char *strFile, int iFeatureDimensionality, int iCovarianceModeling,
			int iHMMStates, int iGaussianComponents, MAccumulatorPhysical &mAccumulatorPhysical);
		
		// load logical accumulators from a file
		static void loadAccumulators(const char *strFile, MAccumulatorLogical &mAccumulatorLogical, AccMetadata &metadata);	
		
		// load physical accumulators from a file
		static void loadAccumulators(const char *strFile, MAccumulatorPhysical &mAccumulatorPhysical, AccMetadata &metadata);
		
		// load and combine physical accumulators from multiple files
		static void loadAccumulatorList(const char *strFileList, MAccumulatorPhysical &mAccumulatorPhysical, AccMetadata &metadata);
		
		// load and combine logical accumulators from multiple files
		static void loadAccumulatorList(const char *strFileList, MAccumulatorLogical &mAccumulatorLogical, AccMetadata &metadata);
		
		// destroy the accumulators
		static void destroy(MAccumulatorLogical &mAccumulatorLogical);
		
		// destroy the accumulators
		static void destroy(MAccumulatorPhysical &mAccumulatorPhysical);	
		
		// print accumulator info
		static void print(MAccumulatorPhysical &mAccumulatorPhysical);
		
		// print accumulator info
		static void print(MAccumulatorLogical &mAccumulatorLogical);		
		
		// accumulate data from the alignment
		static void accumulate(MAccumulatorPhysical &mAccumulator, Alignment *alignment, 
			float *fFeatures, int iFeatures, int iFeatureDimensionality);
			
		// adapt the accumulators to the given within-word and cross-word context length
		static bool adaptContextWidth(MAccumulatorLogical &mAccumulator, unsigned char iContextSize, unsigned char iContextSizeNew);
		
		// allocate memory for the identity
		static unsigned char *newIdentity(int iContextModelingOrder) {
				
			// identity
			int iLength = iContextModelingOrder+3;
			unsigned char *iIdentity = new unsigned char[iLength];
			for(int i=0 ; iIdentity[i] != UCHAR_MAX ; ++i) {
				iIdentity[i] = iIdentity[i];
			}
			iIdentity[iLength-1] = UCHAR_MAX;
			
			return iIdentity;
		}
		
		// shorten the identity
		static unsigned char *getShorterIdentity(unsigned char *iIdentity, unsigned char iContextSize, 
			unsigned char iContextSizeNew) {
			
			int iLength = (iContextSize*2+1)+3;
			unsigned char *iIdentityNew = new unsigned char[iLength];
			// copy the left and right context and the phone too (which is in the middle)
			for(int i=0 ; i<(iContextSize-1)*2+1 ; ++i) {
				iIdentityNew[i] = iIdentity[i+1];
			}
			iIdentity[2*iContextSizeNew+1] = iIdentity[2*iContextSize+1];
			iIdentity[2*iContextSizeNew+2] = iIdentity[2*iContextSize+2];
			iIdentity[2*iContextSizeNew+3] = UCHAR_MAX;
			
			return iIdentityNew;
		}
};

#endif
