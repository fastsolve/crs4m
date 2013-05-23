#ifndef __CRS_PRODAX_H__
#define __CRS_PRODAX_H__
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_prodAx_types.h"
extern void crs_prodAx(const struct_T *A, const m2cArray_real_T *x, m2cArray_real_T *b, const m2cArray_int32_T *nthreads);
extern void crs_prodAx_initialize(void);
extern void crs_prodAx_mpi(const struct_T *A, const m2cArray_real_T *x, m2cArray_real_T *b, const m2cArray_int32_T *nthreads, const b_struct_T *comm);
extern void crs_prodAx_mpip(const struct_T *A, const m2cArray_real_T *x, m2cArray_real_T *b, const m2cArray_int32_T *nthreads, const b_struct_T *comm, const m2cArray_real_T *pbmsg);
extern void crs_prodAx_mpip1(const struct_T *A, const m2cArray_real_T *x, m2cArray_real_T *b, const m2cArray_int32_T *nthreads, const b_struct_T *comm, const m2cArray_real_T *pbmsg, int32_T pbsz);
extern void crs_prodAx_ser(const struct_T *A, const m2cArray_real_T *x, m2cArray_real_T *b);
extern void crs_prodAx_ser1(const struct_T *A, const m2cArray_real_T *x, m2cArray_real_T *b);
extern void crs_prodAx_terminate(void);
extern m2cArray_char_T *m2cCreateND_char_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_int32_T *m2cCreateND_int32_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_real_T *m2cCreateND_real_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_uint8_T *m2cCreateND_uint8_T( m2cShort numDimensions, m2cSize *size);
extern m2cArray_char_T *m2cCreateWrapperND_char_T(char_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_int32_T *m2cCreateWrapperND_int32_T(int32_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_real_T *m2cCreateWrapperND_real_T(real_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_uint8_T *m2cCreateWrapperND_uint8_T(uint8_T *data, m2cShort numDimensions, m2cSize *size);
extern m2cArray_char_T *m2cCreateWrapper_char_T(char_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_int32_T *m2cCreateWrapper_int32_T(int32_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_real_T *m2cCreateWrapper_real_T(real_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_uint8_T *m2cCreateWrapper_uint8_T(uint8_T *data, m2cSize rows, m2cSize cols);
extern m2cArray_char_T *m2cCreate_char_T( m2cSize rows, m2cSize cols);
extern m2cArray_int32_T *m2cCreate_int32_T( m2cSize rows, m2cSize cols);
extern m2cArray_real_T *m2cCreate_real_T( m2cSize rows, m2cSize cols);
extern m2cArray_uint8_T *m2cCreate_uint8_T( m2cSize rows, m2cSize cols);
extern void m2cDestroyArray_char_T(m2cArray_char_T *m2cArray);
extern void m2cDestroyArray_int32_T(m2cArray_int32_T *m2cArray);
extern void m2cDestroyArray_real_T(m2cArray_real_T *m2cArray);
extern void m2cDestroyArray_uint8_T(m2cArray_uint8_T *m2cArray);
#endif
