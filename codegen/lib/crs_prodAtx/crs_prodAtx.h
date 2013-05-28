#ifndef __CRS_PRODATX_H__
#define __CRS_PRODATX_H__
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "plctypes.h"
#include "crs_prodAtx_types.h"
extern void crs_prodAtx(const struct_T *A, const plcArray_real_T *x, plcArray_real_T *b, const plcArray_int32_T *nthreads);
extern void crs_prodAtx_initialize(void);
extern void crs_prodAtx_mpi(const struct_T *A, const plcArray_real_T *x, plcArray_real_T *b, const plcArray_int32_T *nthreads, const b_struct_T *comm);
extern void crs_prodAtx_mpip(const struct_T *A, const plcArray_real_T *x, plcArray_real_T *b, const plcArray_int32_T *nthreads, const b_struct_T *comm, const plcArray_real_T *pbmsg);
extern void crs_prodAtx_mpip1(const struct_T *A, const plcArray_real_T *x, plcArray_real_T *b, const plcArray_int32_T *nthreads, const b_struct_T *comm, const plcArray_real_T *pbmsg, int32_T pbsz);
extern void crs_prodAtx_ser(const struct_T *A, const plcArray_real_T *x, plcArray_real_T *b);
extern void crs_prodAtx_ser1(const struct_T *A, const plcArray_real_T *x, plcArray_real_T *b);
extern void crs_prodAtx_terminate(void);
extern plcArray_char_T *plcCreateND_char_T( plcShort numDimensions, plcSize *size);
extern plcArray_int32_T *plcCreateND_int32_T( plcShort numDimensions, plcSize *size);
extern plcArray_real_T *plcCreateND_real_T( plcShort numDimensions, plcSize *size);
extern plcArray_uint8_T *plcCreateND_uint8_T( plcShort numDimensions, plcSize *size);
extern plcArray_char_T *plcCreateWrapperND_char_T(char_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_int32_T *plcCreateWrapperND_int32_T(int32_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_real_T *plcCreateWrapperND_real_T(real_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_uint8_T *plcCreateWrapperND_uint8_T(uint8_T *data, plcShort numDimensions, plcSize *size);
extern plcArray_char_T *plcCreateWrapper_char_T(char_T *data, plcSize rows, plcSize cols);
extern plcArray_int32_T *plcCreateWrapper_int32_T(int32_T *data, plcSize rows, plcSize cols);
extern plcArray_real_T *plcCreateWrapper_real_T(real_T *data, plcSize rows, plcSize cols);
extern plcArray_uint8_T *plcCreateWrapper_uint8_T(uint8_T *data, plcSize rows, plcSize cols);
extern plcArray_char_T *plcCreate_char_T( plcSize rows, plcSize cols);
extern plcArray_int32_T *plcCreate_int32_T( plcSize rows, plcSize cols);
extern plcArray_real_T *plcCreate_real_T( plcSize rows, plcSize cols);
extern plcArray_uint8_T *plcCreate_uint8_T( plcSize rows, plcSize cols);
extern void plcDestroyArray_char_T(plcArray_char_T *plcArray);
extern void plcDestroyArray_int32_T(plcArray_int32_T *plcArray);
extern void plcDestroyArray_real_T(plcArray_real_T *plcArray);
extern void plcDestroyArray_uint8_T(plcArray_uint8_T *plcArray);
#endif
