#ifndef CUVECCREATE_TYPES_H
#define CUVECCREATE_TYPES_H
#include "rtwtypes.h"
#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T

typedef struct emxArray_int32_T emxArray_int32_T;

#endif

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  unsigned long data;
  emxArray_int32_T *type;
  int len;
} struct0_T;

#endif

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  unsigned long data;
  int type;
  int len;
} struct1_T;

#endif
#endif
