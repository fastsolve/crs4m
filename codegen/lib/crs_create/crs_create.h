#ifndef __CRS_CREATE_H__
#define __CRS_CREATE_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "plctypes.h"
#include "crs_create_types.h"
extern void crs_create(const plcArray_int32_T *rows, const plcArray_int32_T *cols, const plcArray_real_T *vs, struct_T *A);
extern b_struct_T crs_create0(int32_T ni, int32_T nj);
extern void crs_create1(const plcArray_int32_T *is, const plcArray_int32_T *js, const plcArray_real_T *vs, int32_T ni, int32_T nj, struct_T *A);
extern void crs_create_initialize(void);
extern void crs_create_terminate(void);
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
