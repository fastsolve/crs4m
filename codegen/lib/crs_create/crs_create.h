#ifndef __CRS_CREATE_H__
#define __CRS_CREATE_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_create_types.h"
extern void crs_create(const m2cArray_int32_T *rows, const m2cArray_int32_T *cols, const m2cArray_real_T *vs, struct_T *A);
extern b_struct_T crs_create0(int32_T ni, int32_T nj);
extern void crs_create1(const m2cArray_int32_T *is, const m2cArray_int32_T *js, const m2cArray_real_T *vs, int32_T ni, int32_T nj, struct_T *A);
extern void crs_create_initialize(void);
extern void crs_create_terminate(void);
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
