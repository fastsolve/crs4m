#ifndef CUDAMATDESTROY_H
#define CUDAMATDESTROY_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cudaMatDestroy_types.h"

extern void cudaMatDestroy(const struct0_T *mat, int *errCode, boolean_T
  *toplevel);
extern void cudaMatDestroy_initialize(void);
extern void cudaMatDestroy_terminate(void);

#endif
