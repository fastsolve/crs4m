#ifndef PRODATX_H
#define PRODATX_H
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "prodAtx_types.h"

extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                    emxArray_real_T *b, int nthreads);
extern void prodAtx_initialize(void);
extern void prodAtx_ser(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
extern void prodAtx_ser1(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
extern void prodAtx_terminate(void);

#endif
