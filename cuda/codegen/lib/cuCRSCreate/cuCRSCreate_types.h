#ifndef CUCRSCREATE_TYPES_H
#define CUCRSCREATE_TYPES_H
#include "rtwtypes.h"
#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  unsigned long rowptr;
  unsigned long colind;
  unsigned long vals;
  int type;
  int dims[2];
  int nnz;
} struct0_T;

#endif
#endif
