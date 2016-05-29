#ifndef CUDAVECCOPYFROMHOST_H
#define CUDAVECCOPYFROMHOST_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cudaVecCopyFromHost_types.h"

extern void cudaVecCopyFromHost(const emxArray_real_T *vec, const struct0_T
  *cuVec, int *errCode, boolean_T *toplevel);
extern void cudaVecCopyFromHost_initialize(void);
extern void cudaVecCopyFromHost_terminate(void);
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif
