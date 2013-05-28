#ifndef __CRS_DIAG_H__
#define __CRS_DIAG_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "plctypes.h"
#include "crs_diag_types.h"
extern void crs_diag(const struct_T *A, plcArray_real_T *D);
extern void crs_diag1(const struct_T *A, int32_T k, plcArray_real_T *D);
extern void crs_diag_initialize(void);
extern void crs_diag_terminate(void);
extern plcArray_int32_T *plcCreateND_int32_T( plcShort numDimensions, plcSize *size);
extern plcArray_real_T *plcCreateND_real_T( plcShort numDimensions, plcSize *size);
extern plcArray_int32_T *plcCreateWrapperND_int32_T(int32_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_real_T *plcCreateWrapperND_real_T(real_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_int32_T *plcCreateWrapper_int32_T(int32_T *data, plcSize rows, plcSize cols);
extern plcArray_real_T *plcCreateWrapper_real_T(real_T *data, plcSize rows, plcSize cols);
extern plcArray_int32_T *plcCreate_int32_T( plcSize rows, plcSize cols);
extern plcArray_real_T *plcCreate_real_T( plcSize rows, plcSize cols);
extern void plcDestroyArray_int32_T(plcArray_int32_T *plcArray);
extern void plcDestroyArray_real_T(plcArray_real_T *plcArray);
#endif
