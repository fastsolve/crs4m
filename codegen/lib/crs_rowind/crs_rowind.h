#ifndef __CRS_ROWIND_H__
#define __CRS_ROWIND_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "plctypes.h"
#include "crs_rowind_types.h"
extern void crs_rowind(const plcArray_int32_T *row_ptr, const plcArray_int32_T *col_ind, plcArray_int32_T *row_ind);
extern void crs_rowind_initialize(void);
extern void crs_rowind_terminate(void);
extern plcArray_int32_T *plcCreateND_int32_T( plcShort numDimensions, plcSize *size);
extern plcArray_int32_T *plcCreateWrapperND_int32_T(int32_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_int32_T *plcCreateWrapper_int32_T(int32_T *data, plcSize rows, plcSize cols);
extern plcArray_int32_T *plcCreate_int32_T( plcSize rows, plcSize cols);
extern void plcDestroyArray_int32_T(plcArray_int32_T *plcArray);
#endif
