#ifndef CUDAVECDESTROY_H
#define CUDAVECDESTROY_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cudaVecDestroy_types.h"

extern void cudaVecDestroy(const struct0_T *vec, int *errCode, boolean_T
  *toplevel);
extern void cudaVecDestroy_initialize(void);
extern void cudaVecDestroy_terminate(void);

#endif
