#ifndef CUBLASGETENUM_H
#define CUBLASGETENUM_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuBlasGetEnum_types.h"

extern int32_T cuBlasGetEnum(const emxArray_char_T *str);
extern void cuBlasGetEnum_initialize(void);
extern void cuBlasGetEnum_terminate(void);
extern emxArray_char_T *emxCreateND_char_T(int32_T numDimensions, int32_T *size);
extern emxArray_char_T *emxCreateWrapperND_char_T(char_T *data, int32_T
  numDimensions, int32_T *size);
extern emxArray_char_T *emxCreateWrapper_char_T(char_T *data, int32_T rows,
  int32_T cols);
extern emxArray_char_T *emxCreate_char_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_char_T(emxArray_char_T *emxArray);
extern void emxInitArray_char_T(emxArray_char_T **pEmxArray, int32_T
  numDimensions);

#endif
