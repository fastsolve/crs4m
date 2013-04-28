#ifndef __CRS_TRIU_H__
#define __CRS_TRIU_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_triu_types.h"
extern void crs_triu(const struct_T *A, struct_T *U);
extern void crs_triu1(const struct_T *A, int32_T k, struct_T *U);
extern void crs_triu_initialize(void);
extern void crs_triu_terminate(void);
extern emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T cols);
extern emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_int32_T(emxArray_int32_T *emxArray);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
#endif
