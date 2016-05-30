#ifndef CUCRSDESTROY_H
#define CUCRSDESTROY_H
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rtwtypes.h"
#include "cuCRSDestroy_types.h"

extern int cuCRSDestroy(struct0_T *mat);
extern void cuCRSDestroy_initialize(void);
extern void cuCRSDestroy_terminate(void);

#endif
