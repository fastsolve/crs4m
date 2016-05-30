#include "cudaMatCreate.h"
#include "mspack.h"
#include "m2c.h"

#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

struct emxArray_char_T
{
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_char_T
#define typedef_emxArray_char_T

typedef struct emxArray_char_T emxArray_char_T;

#endif

#ifndef struct_emxArray_uint8_T
#define struct_emxArray_uint8_T

struct emxArray_uint8_T
{
  unsigned char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray_uint8_T
#define typedef_emxArray_uint8_T

typedef struct emxArray_uint8_T emxArray_uint8_T;

#endif

static void b_m2c_error(const emxArray_char_T *varargin_3);
static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxFree_char_T(emxArray_char_T **pEmxArray);
static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
static void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);
static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int numDimensions);
static void m2c_error(int varargin_3);
static void b_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int i1;
  int loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i1 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i1, (int)sizeof(char));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_varargin_3->data[i1] = varargin_3->data[i1];
  }

  M2C_error("CUDA:RuntimeError", "cudaMalloc returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

static void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_int32_T(&pStruct->type);
}

static void emxFree_char_T(emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if (((*pEmxArray)->data != (char *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint8_T *)NULL) {
    if (((*pEmxArray)->data != (unsigned char *)NULL) && (*pEmxArray)
        ->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_uint8_T *)NULL;
  }
}

static void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_int32_T(&pStruct->type, 2);
}

static void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxArray_char_T *emxArray;
  int i;
  *pEmxArray = (emxArray_char_T *)malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int numDimensions)
{
  emxArray_uint8_T *emxArray;
  int i;
  *pEmxArray = (emxArray_uint8_T *)malloc(sizeof(emxArray_uint8_T));
  emxArray = *pEmxArray;
  emxArray->data = (unsigned char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void m2c_error(int varargin_3)
{
  M2C_error("cudaMatCreate:WrongType", "Unknow data type %d.\n", varargin_3);
}

void cudaMatCreate(int m, int n, int type, struct0_T *mat, int *errCode,
                   boolean_T *toplevel)
{
  void * t_mat;
  int b_errCode;
  int mat_dims[2];
  unsigned long mat_data;
  signed char mat_type[6];
  int i0;
  static const signed char iv0[6] = { 100, 111, 117, 98, 108, 101 };

  emxArray_uint8_T *msg0;
  const char * ptr;
  int len;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  *toplevel = true;
  if (type == 2) {
    b_errCode = cudaMalloc(&t_mat, m * n << 3);
    mat_dims[0] = m;
    mat_dims[1] = n;
    mat_data = *(uint64_T *)(&t_mat);
    mat->data = mat_data;
    i0 = mat->type->size[0] * mat->type->size[1];
    mat->type->size[0] = 1;
    mat->type->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)mat->type, i0, (int)sizeof(int));
    mat->type->data[0] = 2;
    for (i0 = 0; i0 < 2; i0++) {
      mat->dims[i0] = mat_dims[i0];
    }

    *errCode = b_errCode;
  } else if (type == 1) {
    b_errCode = cudaMalloc(&t_mat, m * n << 2);
    mat_dims[0] = m;
    mat_dims[1] = n;
    mat_data = *(uint64_T *)(&t_mat);
    mat->data = mat_data;
    i0 = mat->type->size[0] * mat->type->size[1];
    mat->type->size[0] = 1;
    mat->type->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)mat->type, i0, (int)sizeof(int));
    mat->type->data[0] = 1;
    for (i0 = 0; i0 < 2; i0++) {
      mat->dims[i0] = mat_dims[i0];
    }

    *errCode = b_errCode;
  } else if (type == 3) {
    b_errCode = cudaMalloc(&t_mat, m * n << 3);
    mat_dims[0] = m;
    mat_dims[1] = n;
    mat_data = *(uint64_T *)(&t_mat);
    mat->data = mat_data;
    i0 = mat->type->size[0] * mat->type->size[1];
    mat->type->size[0] = 1;
    mat->type->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)mat->type, i0, (int)sizeof(int));
    mat->type->data[0] = 3;
    for (i0 = 0; i0 < 2; i0++) {
      mat->dims[i0] = mat_dims[i0];
    }

    *errCode = b_errCode;
  } else if (type == 4) {
    b_errCode = cudaMalloc(&t_mat, m * n << 4);
    mat_dims[0] = m;
    mat_dims[1] = n;
    mat_data = *(uint64_T *)(&t_mat);
    mat->data = mat_data;
    i0 = mat->type->size[0] * mat->type->size[1];
    mat->type->size[0] = 1;
    mat->type->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)mat->type, i0, (int)sizeof(int));
    mat->type->data[0] = 4;
    for (i0 = 0; i0 < 2; i0++) {
      mat->dims[i0] = mat_dims[i0];
    }

    *errCode = b_errCode;
  } else {
    m2c_error(type);
    b_errCode = cudaMalloc(&t_mat, m * n << 3);
    for (i0 = 0; i0 < 6; i0++) {
      mat_type[i0] = iv0[i0];
    }

    mat_dims[0] = m;
    mat_dims[1] = n;
    mat_data = *(uint64_T *)(&t_mat);
    mat->data = mat_data;
    i0 = mat->type->size[0] * mat->type->size[1];
    mat->type->size[0] = 1;
    mat->type->size[1] = 6;
    emxEnsureCapacity((emxArray__common *)mat->type, i0, (int)sizeof(int));
    for (i0 = 0; i0 < 6; i0++) {
      mat->type->data[i0] = mat_type[i0];
    }

    for (i0 = 0; i0 < 2; i0++) {
      mat->dims[i0] = mat_dims[i0];
    }

    *errCode = b_errCode;
  }

  if (*errCode != 0) {
    emxInit_uint8_T(&msg0, 2);
    ptr = cudaGetErrorString(*errCode);
    len = strlen(ptr) + 1;
    i0 = msg0->size[0] * msg0->size[1];
    msg0->size[0] = 1;
    msg0->size[1] = len;
    emxEnsureCapacity((emxArray__common *)msg0, i0, (int)sizeof(unsigned char));
    for (i0 = 0; i0 < len; i0++) {
      msg0->data[i0] = 0;
    }

    memcpy(&msg0->data[0], ptr, len);
    if (1 > len) {
      b_errCode = 0;
    } else {
      b_errCode = len;
    }

    emxInit_uint8_T(&varargin_1, 2);
    i0 = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = 1;
    varargin_1->size[1] = b_errCode;
    emxEnsureCapacity((emxArray__common *)varargin_1, i0, (int)sizeof(unsigned
      char));
    for (i0 = 0; i0 < b_errCode; i0++) {
      varargin_1->data[varargin_1->size[0] * i0] = msg0->data[i0];
    }

    emxFree_uint8_T(&msg0);
    emxInit_char_T(&b_varargin_1, 2);
    i0 = b_varargin_1->size[0] * b_varargin_1->size[1];
    b_varargin_1->size[0] = 1;
    b_varargin_1->size[1] = b_errCode;
    emxEnsureCapacity((emxArray__common *)b_varargin_1, i0, (int)sizeof(char));
    for (i0 = 0; i0 < b_errCode; i0++) {
      b_varargin_1->data[i0] = (signed char)varargin_1->data[i0];
    }

    emxFree_uint8_T(&varargin_1);
    b_m2c_error(b_varargin_1);
    emxFree_char_T(&b_varargin_1);
  }
}

void cudaMatCreate_2args(int m, int n, struct1_T *vec, int *errCode, boolean_T
  *toplevel)
{
  void * t_mat;
  emxArray_uint8_T *msg0;
  const char * ptr;
  int len;
  int i2;
  int loop_ub;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  *errCode = cudaMalloc(&t_mat, m * n << 3);
  vec->type = 2;
  vec->dims[0] = m;
  vec->dims[1] = n;
  vec->data = *(uint64_T *)(&t_mat);
  if (*errCode != 0) {
    emxInit_uint8_T(&msg0, 2);
    ptr = cudaGetErrorString(*errCode);
    len = strlen(ptr) + 1;
    i2 = msg0->size[0] * msg0->size[1];
    msg0->size[0] = 1;
    msg0->size[1] = len;
    emxEnsureCapacity((emxArray__common *)msg0, i2, (int)sizeof(unsigned char));
    for (i2 = 0; i2 < len; i2++) {
      msg0->data[i2] = 0;
    }

    memcpy(&msg0->data[0], ptr, len);
    if (1 > len) {
      loop_ub = 0;
    } else {
      loop_ub = len;
    }

    emxInit_uint8_T(&varargin_1, 2);
    i2 = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = 1;
    varargin_1->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)varargin_1, i2, (int)sizeof(unsigned
      char));
    for (i2 = 0; i2 < loop_ub; i2++) {
      varargin_1->data[varargin_1->size[0] * i2] = msg0->data[i2];
    }

    emxFree_uint8_T(&msg0);
    emxInit_char_T(&b_varargin_1, 2);
    i2 = b_varargin_1->size[0] * b_varargin_1->size[1];
    b_varargin_1->size[0] = 1;
    b_varargin_1->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_varargin_1, i2, (int)sizeof(char));
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_varargin_1->data[i2] = (signed char)varargin_1->data[i2];
    }

    emxFree_uint8_T(&varargin_1);
    b_m2c_error(b_varargin_1);
    emxFree_char_T(&b_varargin_1);
  }

  *toplevel = true;
}

void cudaMatCreate_initialize(void)
{
}

void cudaMatCreate_terminate(void)
{
}

void emxDestroy_struct0_T(struct0_T emxArray)
{
  emxFreeStruct_struct0_T(&emxArray);
}

void emxInit_struct0_T(struct0_T *pStruct)
{
  emxInitStruct_struct0_T(pStruct);
}
