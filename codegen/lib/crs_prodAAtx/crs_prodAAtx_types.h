#ifndef __CRS_PRODAATX_TYPES_H__
#define __CRS_PRODAATX_TYPES_H__
#include "plctypes.h"
#ifndef struct_plcArray_char_T
#define struct_plcArray_char_T
struct plcArray_char_T
{
    char_T *data;
    plcSize *size;
    plcSize allocatedSize;
    plcShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_plcArray_char_T*/
#ifndef typedef_plcArray_char_T
#define typedef_plcArray_char_T
typedef struct plcArray_char_T plcArray_char_T;
#endif /*typedef_plcArray_char_T*/
#ifndef struct_plcArray_uint8_T
#define struct_plcArray_uint8_T
struct plcArray_uint8_T
{
    uint8_T *data;
    plcSize *size;
    plcSize allocatedSize;
    plcShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_plcArray_uint8_T*/
#ifndef typedef_plcArray_uint8_T
#define typedef_plcArray_uint8_T
typedef struct plcArray_uint8_T plcArray_uint8_T;
#endif /*typedef_plcArray_uint8_T*/
#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct
{
    plcArray_uint8_T *data;
    plcArray_char_T *type;
    int32_T nitems;
} b_struct_T;
#endif /*typedef_b_struct_T*/
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
#ifndef struct_plcArray_real_T
#define struct_plcArray_real_T
struct plcArray_real_T
{
    real_T *data;
    plcSize *size;
    plcSize allocatedSize;
    plcShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_plcArray_real_T*/
#ifndef typedef_plcArray_real_T
#define typedef_plcArray_real_T
typedef struct plcArray_real_T plcArray_real_T;
#endif /*typedef_plcArray_real_T*/
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct
{
    plcArray_int32_T *row_ptr;
    plcArray_int32_T *col_ind;
    plcArray_real_T *val;
    int32_T nrows;
    int32_T ncols;
} struct_T;
#endif /*typedef_struct_T*/

#endif
