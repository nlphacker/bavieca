/* BaviecaAPI.i */
%module BaviecaAPI_SWIG
%{
#include "BaviecaAPI.h"
%}


// typemaps for function parameters
%include <typemaps.i>
%apply float *INOUT { float *fFeatures }
%apply short *INOUT { short *sSamples }
%apply unsigned int *INOUT { unsigned int *iFeatures }
%apply int *INOUT { int *iWords }
%apply int *INOUT { int *iSegments }
%apply int *INOUT { int *iTAElements }

# typemap for return parameters (for example extractFeatures())
%typemap(jstype) float *extractFeatures "float[]"
%typemap(jtype) float *extractFeatures "float[]"

%typemap(jni) float *extractFeatures "jfloatArray"
%typemap(javaout) float *extractFeatures {
  return $jnicall;
}

%typemap(in,numinputs=0,noblock=1) int *iFeatures {
   int temp_iFeatures;
   $1 = &temp_iFeatures;
}

%typemap(out) float *extractFeatures {
  $result = JCALL1(NewFloatArray, jenv, (jsize)(*arg4*(arg1)->getFeatureDim()));
  JCALL4(SetFloatArrayRegion, jenv, $result, 0, (jsize)(*arg4*(arg1)->getFeatureDim()), $1);
  delete [] $1;
}

%include "BaviecaAPI.h"

//  (jsize)(*arg4*(arg1)->getFeatureDim())

