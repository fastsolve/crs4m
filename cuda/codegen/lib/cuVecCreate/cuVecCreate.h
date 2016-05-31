#ifndef CUVECCREATE_H
#define CUVECCREATE_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuVecCreate_types.h"

extern void cuVecCreate(int32_T n, int32_T type, struct0_T *vec, int32_T
  *errCode, boolean_T *toplevel);
extern void cuVecCreate_1arg(int32_T n, struct0_T *vec, int32_T *errCode,
  boolean_T *toplevel);
extern void cuVecCreate_initialize(void);
extern void cuVecCreate_terminate(void);

#endif
