#ifndef __CRS_ROWIND_TYPES_H__
#define __CRS_ROWIND_TYPES_H__
#include "plctypes.h"
#ifndef struct_plcArray_int32_T
#define struct_plcArray_int32_T
struct plcArray_int32_T
{
    int32_T *data;
    plcSize *size;
    plcSize allocatedSize;
    plcShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_plcArray_int32_T*/
#ifndef typedef_plcArray_int32_T
#define typedef_plcArray_int32_T
typedef struct plcArray_int32_T plcArray_int32_T;
#endif /*typedef_plcArray_int32_T*/

#endif
