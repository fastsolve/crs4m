#ifndef CUBCRSDESTROY_TYPES_H
#define CUBCRSDESTROY_TYPES_H
#include "rtwtypes.h"
#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  uint64_T rowptr;
  uint64_T colind;
  uint64_T vals;
  int32_T type;
  int32_T dimsb[2];
  int32_T nnzb;
  int32_T blkdim;
} struct0_T;

#endif
#endif
