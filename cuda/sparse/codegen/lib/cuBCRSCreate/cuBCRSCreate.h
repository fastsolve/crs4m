#ifndef CUBCRSCREATE_H
#define CUBCRSCREATE_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuBCRSCreate_types.h"

extern void cuBCRSCreate(int32_T mb, int32_T nb, int32_T nnzb, int32_T blkdim,
  int32_T type, MSP_CuBCRS *mat, int32_T *errCode);
extern void cuBCRSCreate_3args(int32_T mb, int32_T nb, int32_T nnzb, MSP_CuBCRS *
  mat, int32_T *errCode);
extern void cuBCRSCreate_4args(int32_T mb, int32_T nb, int32_T nnzb, int32_T
  blkdim, MSP_CuBCRS *mat, int32_T *errCode);
extern void cuBCRSCreate_initialize(void);
extern void cuBCRSCreate_terminate(void);

#endif
