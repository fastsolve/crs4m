#ifndef __CRS_SORT_H__
#define __CRS_SORT_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_sort_types.h"
extern void crs_sort(const m2cArray_int32_T *row_ptr, m2cArray_int32_T *col_ind, m2cArray_real_T *val);
extern void crs_sort0(const m2cArray_int32_T *row_ptr, m2cArray_int32_T *col_ind);
extern void crs_sort_initialize(void);
extern void crs_sort_terminate(void);
extern m2cArray_int32_T *m2cCreateND_int32_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_real_T *m2cCreateND_real_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_int32_T *m2cCreateWrapperND_int32_T(int32_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_real_T *m2cCreateWrapperND_real_T(real_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_int32_T *m2cCreateWrapper_int32_T(int32_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_real_T *m2cCreateWrapper_real_T(real_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_int32_T *m2cCreate_int32_T( m2cSize rows, m2cSize cols);
extern m2cArray_real_T *m2cCreate_real_T( m2cSize rows, m2cSize cols);
extern void m2cDestroyArray_int32_T(m2cArray_int32_T *m2cArray);
extern void m2cDestroyArray_real_T(m2cArray_real_T *m2cArray);
#endif
