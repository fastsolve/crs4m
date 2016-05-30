#ifndef CUVECCREATE_H
#define CUVECCREATE_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuVecCreate_types.h"

extern void cuVecCreate(int n, int type, struct0_T *vec, int *errCode, boolean_T
  *toplevel);
extern void cuVecCreate_1arg(int n, struct1_T *vec, int *errCode, boolean_T
  *toplevel);
extern void cuVecCreate_initialize(void);
extern void cuVecCreate_terminate(void);
extern emxArray_int32_T *emxCreateND_int32_T(int numDimensions, int *size);
extern emxArray_int32_T *emxCreateWrapperND_int32_T(int *data, int numDimensions,
  int *size);
extern emxArray_int32_T *emxCreateWrapper_int32_T(int *data, int rows, int cols);
extern emxArray_int32_T *emxCreate_int32_T(int rows, int cols);
extern void emxDestroyArray_int32_T(emxArray_int32_T *emxArray);
extern void emxDestroy_struct0_T(struct0_T emxArray);
extern void emxInit_struct0_T(struct0_T *pStruct);

#endif
