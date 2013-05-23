#ifndef __CRS_ROWIND_H__
#define __CRS_ROWIND_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_rowind_types.h"
extern void crs_rowind(const m2cArray_int32_T *row_ptr, const m2cArray_int32_T *col_ind, m2cArray_int32_T *row_ind);
extern void crs_rowind_initialize(void);
extern void crs_rowind_terminate(void);
extern m2cArray_int32_T *m2cCreateND_int32_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_int32_T *m2cCreateWrapperND_int32_T(int32_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_int32_T *m2cCreateWrapper_int32_T(int32_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_int32_T *m2cCreate_int32_T( m2cSize rows, m2cSize cols);
extern void m2cDestroyArray_int32_T(m2cArray_int32_T *m2cArray);
#endif
