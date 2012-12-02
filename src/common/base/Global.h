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


#ifndef GLOBAL_H
#define GLOBAL_H

/**
	@author root <root@localhost.localdomain>
*/

// author, version and date
#define SYSTEM_AUTHOR       "Daniel Bolanos"
#define SYSTEM_VERSION      "0014"
#define SYSTEM_DATE         "July 1st 2012" 

#define SYSTEM_VERSION_FIELD_WIDTH	  	8	

// HMM-state topology
#define NUMBER_HMM_STATE_POSITIONS		4			// beg/mid/end/monophone
#define NUMBER_HMM_STATES				3

// file format
#define FILE_FORMAT_TEXT		0	
#define FILE_FORMAT_BINARY		1

#define LOG_LIKELIHOOD_FLOOR		-70.0

// emission probability computation
#define OPTIMIZED_COMPUTATION					// precomputed constants and inverted covariance
//#define SIMD										// special memory alignment to work with simd instructions

// open-mp
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#undef min
#undef max
#include <float.h>
#define finite _finite
#endif

// asserts
//#define NDEBUG 
#include <assert.h> 
#include <math.h>
#include <float.h>

const double PI_NUMBER = 2.0*acos(0.0);

// inlining
#if defined __linux__
     #define FORCE_INLINE __attribute__((always_inline))
     #define NO_INLINE __attribute__((noinline))
#elif defined _WIN32
     #define FORCE_INLINE __forceinline
     #define NO_INLINE __declspec(noinline)
#else
     #error "no inlining support"
     #define FORCE_INLINE inline
     #define NO_INLINE
#endif 


#endif
