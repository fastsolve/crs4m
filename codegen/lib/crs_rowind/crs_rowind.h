#ifndef __CRS_ROWIND_H__
#define __CRS_ROWIND_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_rowind_types.h"
extern void crs_rowind(const emxArray_int32_T *row_ptr, const emxArray_int32_T *col_ind, emxArray_int32_T *row_ind);
extern void crs_rowind_initialize(void);
extern void crs_rowind_terminate(void);
extern emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T cols);
extern emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_int32_T(emxArray_int32_T *emxArray);
#endif
