#ifndef CUMATCOPYSUBTOGPU_TYPES_H
#define CUMATCOPYSUBTOGPU_TYPES_H
#include "rtwtypes.h"
#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

struct emxArray_char_T
{
  char_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_char_T
#define typedef_emxArray_char_T

typedef struct emxArray_char_T emxArray_char_T;

#endif

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif

#ifndef struct_emxArray_uint8_T
#define struct_emxArray_uint8_T

struct emxArray_uint8_T
{
  uint8_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_uint8_T
#define typedef_emxArray_uint8_T

typedef struct emxArray_uint8_T emxArray_uint8_T;

#endif

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  uint64_T data;
  int32_T type;
  int32_T dims[2];
} struct0_T;

#endif

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  emxArray_uint8_T *data;
  emxArray_char_T *type;
  int32_T nitems;
} struct1_T;

#endif
#endif
