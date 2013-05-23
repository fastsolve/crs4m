#ifndef __CRS_ROWIND_TYPES_H__
#define __CRS_ROWIND_TYPES_H__
#include "rtwtypes.h"
#ifndef struct_m2cArray_int32_T
#define struct_m2cArray_int32_T
struct m2cArray_int32_T
{
    int32_T *data;
    m2cSize *size;
    m2cSize allocatedSize;
    m2cShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_m2cArray_int32_T*/
#ifndef typedef_m2cArray_int32_T
#define typedef_m2cArray_int32_T
typedef struct m2cArray_int32_T m2cArray_int32_T;
#endif /*typedef_m2cArray_int32_T*/

#endif
