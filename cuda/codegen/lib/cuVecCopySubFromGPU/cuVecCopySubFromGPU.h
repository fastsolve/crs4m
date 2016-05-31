#ifndef CUVECCOPYSUBFROMGPU_H
#define CUVECCOPYSUBFROMGPU_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuVecCopySubFromGPU_types.h"

extern void cuVecCopySubFromGPU(int32_T n, const struct0_T *cuVec, int32_T
  incCuVec, const emxArray_real_T *vec, int32_T istart, int32_T inc, int32_T
  *errCode, boolean_T *toplevel);
extern void cuVecCopySubFromGPU_async(int32_T n, const struct0_T *cuVec, int32_T
  incCuVec, const emxArray_real_T *vec, int32_T istart, int32_T inc, const
  struct1_T *strm, int32_T *errCode, boolean_T *toplevel);
extern void cuVecCopySubFromGPU_initialize(void);
extern void cuVecCopySubFromGPU_terminate(void);
extern emxArray_char_T *emxCreateND_char_T(int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size);
extern emxArray_uint8_T *emxCreateND_uint8_T(int32_T numDimensions, int32_T
  *size);
extern emxArray_char_T *emxCreateWrapperND_char_T(char_T *data, int32_T
  numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T
  numDimensions, int32_T *size);
extern emxArray_uint8_T *emxCreateWrapperND_uint8_T(uint8_T *data, int32_T
  numDimensions, int32_T *size);
extern emxArray_char_T *emxCreateWrapper_char_T(char_T *data, int32_T rows,
  int32_T cols);
extern emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows,
  int32_T cols);
extern emxArray_uint8_T *emxCreateWrapper_uint8_T(uint8_T *data, int32_T rows,
  int32_T cols);
extern emxArray_char_T *emxCreate_char_T(int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols);
extern emxArray_uint8_T *emxCreate_uint8_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_char_T(emxArray_char_T *emxArray);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxDestroyArray_uint8_T(emxArray_uint8_T *emxArray);
extern void emxDestroy_struct1_T(struct1_T emxArray);
extern void emxInit_struct1_T(struct1_T *pStruct);

#endif
