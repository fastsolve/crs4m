#ifndef CUBCRSCREATE_TYPES_H
#define CUBCRSCREATE_TYPES_H
#include "rtwtypes.h"
#ifndef typedef_MSP_CuBCRS
#define typedef_MSP_CuBCRS

typedef struct {
  uint64_T rowptr;
  uint64_T colind;
  uint64_T vals;
  int32_T type;
  int32_T dimsb[2];
  int32_T nnzb;
  int32_T blkdim;
} MSP_CuBCRS;

#endif
#endif
