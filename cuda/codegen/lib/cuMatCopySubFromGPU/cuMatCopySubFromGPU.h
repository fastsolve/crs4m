#ifndef CUMATCOPYSUBFROMGPU_H
#define CUMATCOPYSUBFROMGPU_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuMatCopySubFromGPU_types.h"

extern void cuMatCopySubFromGPU(int nrows, int ncols, const struct0_T *cuMat,
  const emxArray_real_T *mat, int *errCode, boolean_T *toplevel);
extern void cuMatCopySubFromGPU_initialize(void);
extern void cuMatCopySubFromGPU_terminate(void);
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);

#endif
