#ifndef CUMATCOPYSUBTOGPU_H
#define CUMATCOPYSUBTOGPU_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuMatCopySubToGPU_types.h"

extern void cuMatCopySubToGPU(int32_T nrows, int32_T ncols, const
  emxArray_real_T *mat, const struct0_T *cuMat, int32_T *errCode, boolean_T
  *toplevel);
extern void cuMatCopySubToGPU_async(int32_T nrows, int32_T ncols, const
  emxArray_real_T *mat, const struct0_T *cuMat, const struct1_T *strm, int32_T
  *errCode, boolean_T *toplevel);
extern void cuMatCopySubToGPU_initialize(void);
extern void cuMatCopySubToGPU_terminate(void);
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
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int32_T
  numDimensions);
extern void emxInit_struct1_T(struct1_T *pStruct);

#endif
