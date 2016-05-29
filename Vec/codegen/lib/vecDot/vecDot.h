#ifndef VECDOT_H
#define VECDOT_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "vecDot_types.h"

extern emxArray_char_T *emxCreateND_char_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_uint8_T *emxCreateND_uint8_T(int numDimensions, int *size);
extern emxArray_char_T *emxCreateWrapperND_char_T(char *data, int numDimensions,
  int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_uint8_T *emxCreateWrapperND_uint8_T(unsigned char *data, int
  numDimensions, int *size);
extern emxArray_char_T *emxCreateWrapper_char_T(char *data, int rows, int cols);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_uint8_T *emxCreateWrapper_uint8_T(unsigned char *data, int rows,
  int cols);
extern emxArray_char_T *emxCreate_char_T(int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern emxArray_uint8_T *emxCreate_uint8_T(int rows, int cols);
extern void emxDestroyArray_char_T(emxArray_char_T *emxArray);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxDestroyArray_uint8_T(emxArray_uint8_T *emxArray);
extern void emxDestroy_struct1_T(struct1_T emxArray);
extern void emxInitArray_char_T(emxArray_char_T **pEmxArray, int numDimensions);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_struct1_T(struct1_T *pStruct);
extern void vecDot(const emxArray_real_T *x, const emxArray_real_T *y, const
                   emxArray_real_T *buf, double *prod, boolean_T *toplevel);
extern double vecDot_cublas(const struct0_T *u, const struct0_T *v, const
  struct0_T *buf, const emxArray_char_T *mode, const struct1_T *cublasHdl);
extern double vecDot_cublas_sync(const struct0_T *u, const struct0_T *v, const
  emxArray_real_T *buf, const emxArray_char_T *mode, const struct1_T *cublasHdl,
  const emxArray_char_T *sync);
extern void vecDot_initialize(void);
extern double vecDot_omp(const emxArray_real_T *u, const emxArray_real_T *v,
  emxArray_real_T *buf, int nthreads);
extern double vecDot_ser(const emxArray_real_T *u, const emxArray_real_T *v);
extern void vecDot_terminate(void);

#endif
