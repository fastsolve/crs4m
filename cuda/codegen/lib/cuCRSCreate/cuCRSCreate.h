#ifndef CUCRSCREATE_H
#define CUCRSCREATE_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuCRSCreate_types.h"

extern void cuCRSCreate(int m, int n, int nnz, int type, struct0_T *mat, int
  *errCode);
extern void cuCRSCreate_3args(int m, int n, int nnz, struct0_T *mat, int
  *errCode);
extern void cuCRSCreate_initialize(void);
extern void cuCRSCreate_terminate(void);

#endif
