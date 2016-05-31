#ifndef CUBLASDESTROY_H
#define CUBLASDESTROY_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuBlasDestroy_types.h"

extern void cuBlasDestroy(const struct0_T *hdl, int32_T *errCode, boolean_T
  *toplevel);
extern void cuBlasDestroy_initialize(void);
extern void cuBlasDestroy_terminate(void);
extern emxArray_char_T *emxCreateND_char_T(int32_T numDimensions, int32_T *size);
extern emxArray_uint8_T *emxCreateND_uint8_T(int32_T numDimensions, int32_T
  *size);
extern emxArray_char_T *emxCreateWrapperND_char_T(char_T *data, int32_T
  numDimensions, int32_T *size);
extern emxArray_uint8_T *emxCreateWrapperND_uint8_T(uint8_T *data, int32_T
  numDimensions, int32_T *size);
extern emxArray_char_T *emxCreateWrapper_char_T(char_T *data, int32_T rows,
  int32_T cols);
extern emxArray_uint8_T *emxCreateWrapper_uint8_T(uint8_T *data, int32_T rows,
  int32_T cols);
extern emxArray_char_T *emxCreate_char_T(int32_T rows, int32_T cols);
extern emxArray_uint8_T *emxCreate_uint8_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_char_T(emxArray_char_T *emxArray);
extern void emxDestroyArray_uint8_T(emxArray_uint8_T *emxArray);
extern void emxDestroy_struct0_T(struct0_T emxArray);
extern void emxInit_struct0_T(struct0_T *pStruct);

#endif
