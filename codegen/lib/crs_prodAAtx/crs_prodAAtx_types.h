#ifndef __CRS_PRODAATX_TYPES_H__
#define __CRS_PRODAATX_TYPES_H__
#include "rtwtypes.h"
#ifndef struct_m2cArray_char_T
#define struct_m2cArray_char_T
struct m2cArray_char_T
{
    char_T *data;
    m2cSize *size;
    m2cSize allocatedSize;
    m2cShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_m2cArray_char_T*/
#ifndef typedef_m2cArray_char_T
#define typedef_m2cArray_char_T
typedef struct m2cArray_char_T m2cArray_char_T;
#endif /*typedef_m2cArray_char_T*/
#ifndef struct_m2cArray_uint8_T
#define struct_m2cArray_uint8_T
struct m2cArray_uint8_T
{
    uint8_T *data;
    m2cSize *size;
    m2cSize allocatedSize;
    m2cShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_m2cArray_uint8_T*/
#ifndef typedef_m2cArray_uint8_T
#define typedef_m2cArray_uint8_T
typedef struct m2cArray_uint8_T m2cArray_uint8_T;
#endif /*typedef_m2cArray_uint8_T*/
#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct
{
    m2cArray_uint8_T *data;
    m2cArray_char_T *type;
    int32_T nitems;
} b_struct_T;
#endif /*typedef_b_struct_T*/
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
#ifndef struct_m2cArray_real_T
#define struct_m2cArray_real_T
struct m2cArray_real_T
{
    real_T *data;
    m2cSize *size;
    m2cSize allocatedSize;
    m2cShort numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_m2cArray_real_T*/
#ifndef typedef_m2cArray_real_T
#define typedef_m2cArray_real_T
typedef struct m2cArray_real_T m2cArray_real_T;
#endif /*typedef_m2cArray_real_T*/
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct
{
    m2cArray_int32_T *row_ptr;
    m2cArray_int32_T *col_ind;
    m2cArray_real_T *val;
    int32_T nrows;
    int32_T ncols;
} struct_T;
#endif /*typedef_struct_T*/

#endif
