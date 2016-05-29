#ifndef CUDAVECCOPYFROMHOST_TYPES_H
#define CUDAVECCOPYFROMHOST_TYPES_H
#include "rtwtypes.h"
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  unsigned int data[2];
  int type;
  int len;
} struct0_T;

#endif
#endif