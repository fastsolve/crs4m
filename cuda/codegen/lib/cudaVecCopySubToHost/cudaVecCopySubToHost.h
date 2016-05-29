#ifndef CUDAVECCOPYSUBTOHOST_H
#define CUDAVECCOPYSUBTOHOST_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cudaVecCopySubToHost_types.h"

extern void cudaVecCopySubToHost(int n, const struct0_T *cuVec, int incCuVec,
  const emxArray_real_T *vec, int istart, int inc, int *errCode, boolean_T
  *toplevel);
extern void cudaVecCopySubToHost_initialize(void);
extern void cudaVecCopySubToHost_terminate(void);
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);

#endif
