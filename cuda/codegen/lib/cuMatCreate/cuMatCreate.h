#ifndef CUMATCREATE_H
#define CUMATCREATE_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuMatCreate_types.h"

extern void cuMatCreate(int32_T m, int32_T n, int32_T type, struct0_T *mat,
  int32_T *errCode, boolean_T *toplevel);
extern void cuMatCreate_2args(int32_T m, int32_T n, struct0_T *vec, int32_T
  *errCode, boolean_T *toplevel);
extern void cuMatCreate_initialize(void);
extern void cuMatCreate_terminate(void);

#endif
