#ifndef __CRS_PRODAX_H__
#define __CRS_PRODAX_H__
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "rtwtypes.h"
#include "crs_prodAx_types.h"
extern void crs_prodAx(const struct_T A, const emxArray_real_T *x, emxArray_real_T *b, int32_T nthreads);
extern void crs_prodAx_initialize(void);
extern void crs_prodAx_mpi(const struct_T A, const emxArray_real_T *x, emxArray_real_T *b, int32_T nthreads, const b_struct_T comm);
extern void crs_prodAx_mpip(const struct_T A, const emxArray_real_T *x, emxArray_real_T *b, int32_T nthreads, const b_struct_T comm, const emxArray_real_T *pbmsg);
extern void crs_prodAx_mpip1(const struct_T A, const emxArray_real_T *x, emxArray_real_T *b, int32_T nthreads, const b_struct_T comm, const emxArray_real_T *pbmsg, int32_T pbsz);
extern void crs_prodAx_ser(const struct_T A, const emxArray_real_T *x, emxArray_real_T *b);
extern void crs_prodAx_ser1(const struct_T A, const emxArray_real_T *x, emxArray_real_T *b);
extern void crs_prodAx_terminate(void);
extern emxArray_char_T *emxCreateND_char_T(int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size);
extern emxArray_uint8_T *emxCreateND_uint8_T(int32_T numDimensions, int32_T *size);
extern emxArray_char_T *emxCreateWrapperND_char_T(char_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_uint8_T *emxCreateWrapperND_uint8_T(uint8_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_char_T *emxCreateWrapper_char_T(char_T *data, int32_T rows, int32_T cols);
extern emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T cols);
extern emxArray_uint8_T *emxCreateWrapper_uint8_T(uint8_T *data, int32_T rows, int32_T cols);
extern emxArray_char_T *emxCreate_char_T(int32_T rows, int32_T cols);
extern emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols);
extern emxArray_uint8_T *emxCreate_uint8_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_char_T(emxArray_char_T *emxArray);
extern void emxDestroyArray_int32_T(emxArray_int32_T *emxArray);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxDestroyArray_uint8_T(emxArray_uint8_T *emxArray);
#endif
