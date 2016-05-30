#ifndef CUDAMATCOPYSUBFROMHOST_H
#define CUDAMATCOPYSUBFROMHOST_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cudaMatCopySubFromHost_types.h"

extern void cudaMatCopySubFromHost(int nrows, int ncols, const emxArray_real_T
  *mat, const struct0_T *cuMat, int *errCode, boolean_T *toplevel);
extern void cudaMatCopySubFromHost_initialize(void);
extern void cudaMatCopySubFromHost_terminate(void);
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif
