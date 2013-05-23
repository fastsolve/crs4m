#ifndef __CRS_PRODAB_H__
#define __CRS_PRODAB_H__
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_prodAB_types.h"
extern void crs_prodAB(const struct_T *A, const struct_T *B, struct_T *C);
extern void crs_prodAB_initialize(void);
extern void crs_prodAB_terminate(void);
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
