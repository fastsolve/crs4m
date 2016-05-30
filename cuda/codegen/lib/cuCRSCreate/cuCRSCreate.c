#include "cuCRSCreate.h"
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
static void cuVecCreate(int n, unsigned long *vec_data, int *vec_type, int
  *vec_len, int *errCode);
static void emxFree_char_T(emxArray_char_T **pEmxArray);
static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray);
static void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);
static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int numDimensions);
static void m2c_error(int varargin_3);
static void b_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int i2;
  int loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i2 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i2, (int)sizeof(char));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    b_varargin_3->data[i2] = varargin_3->data[i2];
  }

  M2C_error("CUDA:RuntimeError", "cudaMalloc returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

static void cuVecCreate(int n, unsigned long *vec_data, int *vec_type, int
  *vec_len, int *errCode)
{
  void * t_vec;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int flag;
  const char * ptr;
  int len;
  int i1;
  *errCode = cudaMalloc(&t_vec, n << 2);
  *vec_type = 13;
  *vec_len = n;
  *vec_data = *(uint64_T *)(&t_vec);
  emxInit_uint8_T(&msg0, 2);
  emxInit_uint8_T(&varargin_1, 2);
  emxInit_char_T(&b_varargin_1, 2);
  if (*errCode != 0) {
    flag = (M2C_DEBUG);
    if (flag != 0) {
      ptr = cudaGetErrorString(*errCode);
      len = strlen(ptr) + 1;
      i1 = msg0->size[0] * msg0->size[1];
      msg0->size[0] = 1;
      msg0->size[1] = len;
      emxEnsureCapacity((emxArray__common *)msg0, i1, (int)sizeof(unsigned char));
      for (i1 = 0; i1 < len; i1++) {
        msg0->data[i1] = 0;
      }

      memcpy(&msg0->data[0], ptr, len);
      if (1 > len) {
        flag = 0;
      } else {
        flag = len;
      }

      i1 = varargin_1->size[0] * varargin_1->size[1];
      varargin_1->size[0] = 1;
      varargin_1->size[1] = flag;
      emxEnsureCapacity((emxArray__common *)varargin_1, i1, (int)sizeof(unsigned
        char));
      for (i1 = 0; i1 < flag; i1++) {
        varargin_1->data[varargin_1->size[0] * i1] = msg0->data[i1];
      }

      i1 = b_varargin_1->size[0] * b_varargin_1->size[1];
      b_varargin_1->size[0] = 1;
      b_varargin_1->size[1] = flag;
      emxEnsureCapacity((emxArray__common *)b_varargin_1, i1, (int)sizeof(char));
      for (i1 = 0; i1 < flag; i1++) {
        b_varargin_1->data[i1] = (signed char)varargin_1->data[i1];
      }

      b_m2c_error(b_varargin_1);
    }
  }

  emxFree_char_T(&b_varargin_1);
  emxFree_uint8_T(&varargin_1);
  emxFree_uint8_T(&msg0);
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
  M2C_error("mspGetSizePerElement:WrongType", "Unknow data type %d.\n",
            varargin_3);
}

void cuCRSCreate(int m, int n, int nnz, int type, struct0_T *mat, int *errCode)
{
  unsigned long rowptr_data;
  int expl_temp;
  int b_expl_temp;
  unsigned long colind_data;
  unsigned long vals_data;
  int sizepe;
  void * t_vec;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int flag;
  const char * ptr;
  int len;
  int i0;
  cuVecCreate(m + 1, &rowptr_data, &expl_temp, &b_expl_temp, errCode);
  if (!(*errCode != 0)) {
    cuVecCreate(nnz, &colind_data, &expl_temp, &b_expl_temp, errCode);
  } else {
    colind_data = 0UL;
  }

  if (!(*errCode != 0)) {
    if ((type == 2) || (type == -1) || (type == 14)) {
      sizepe = 8;
    } else if ((type == 1) || (type == 13) || (type == 13)) {
      sizepe = 4;
    } else if (type == 12) {
      sizepe = 2;
    } else if (type == 11) {
      sizepe = 1;
    } else if (type == -2) {
      sizepe = 16;
    } else {
      m2c_error(type);
      sizepe = 0;
    }

    *errCode = cudaMalloc(&t_vec, nnz * sizepe);
    vals_data = *(uint64_T *)(&t_vec);
    emxInit_uint8_T(&msg0, 2);
    emxInit_uint8_T(&varargin_1, 2);
    emxInit_char_T(&b_varargin_1, 2);
    if (*errCode != 0) {
      flag = (M2C_DEBUG);
      if (flag != 0) {
        ptr = cudaGetErrorString(*errCode);
        len = strlen(ptr) + 1;
        i0 = msg0->size[0] * msg0->size[1];
        msg0->size[0] = 1;
        msg0->size[1] = len;
        emxEnsureCapacity((emxArray__common *)msg0, i0, (int)sizeof(unsigned
          char));
        for (i0 = 0; i0 < len; i0++) {
          msg0->data[i0] = 0;
        }

        memcpy(&msg0->data[0], ptr, len);
        if (1 > len) {
          flag = 0;
        } else {
          flag = len;
        }

        i0 = varargin_1->size[0] * varargin_1->size[1];
        varargin_1->size[0] = 1;
        varargin_1->size[1] = flag;
        emxEnsureCapacity((emxArray__common *)varargin_1, i0, (int)sizeof
                          (unsigned char));
        for (i0 = 0; i0 < flag; i0++) {
          varargin_1->data[varargin_1->size[0] * i0] = msg0->data[i0];
        }

        i0 = b_varargin_1->size[0] * b_varargin_1->size[1];
        b_varargin_1->size[0] = 1;
        b_varargin_1->size[1] = flag;
        emxEnsureCapacity((emxArray__common *)b_varargin_1, i0, (int)sizeof(char));
        for (i0 = 0; i0 < flag; i0++) {
          b_varargin_1->data[i0] = (signed char)varargin_1->data[i0];
        }

        b_m2c_error(b_varargin_1);
      }
    }

    emxFree_char_T(&b_varargin_1);
    emxFree_uint8_T(&varargin_1);
    emxFree_uint8_T(&msg0);
  } else {
    vals_data = 0UL;
  }

  mat->rowptr = rowptr_data;
  mat->colind = colind_data;
  mat->vals = vals_data;
  mat->type = type;
  mat->dims[0] = m;
  mat->dims[1] = n;
  mat->nnz = nnz;
}

void cuCRSCreate_3args(int m, int n, int nnz, struct0_T *mat, int *errCode)
{
  unsigned long rowptr_data;
  int expl_temp;
  int b_expl_temp;
  unsigned long colind_data;
  unsigned long vals_data;
  void * t_vec;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int flag;
  const char * ptr;
  int len;
  int i3;
  cuVecCreate(m + 1, &rowptr_data, &expl_temp, &b_expl_temp, errCode);
  if (!(*errCode != 0)) {
    cuVecCreate(nnz, &colind_data, &expl_temp, &b_expl_temp, errCode);
  } else {
    colind_data = 0UL;
  }

  if (!(*errCode != 0)) {
    *errCode = cudaMalloc(&t_vec, nnz << 3);
    vals_data = *(uint64_T *)(&t_vec);
    emxInit_uint8_T(&msg0, 2);
    emxInit_uint8_T(&varargin_1, 2);
    emxInit_char_T(&b_varargin_1, 2);
    if (*errCode != 0) {
      flag = (M2C_DEBUG);
      if (flag != 0) {
        ptr = cudaGetErrorString(*errCode);
        len = strlen(ptr) + 1;
        i3 = msg0->size[0] * msg0->size[1];
        msg0->size[0] = 1;
        msg0->size[1] = len;
        emxEnsureCapacity((emxArray__common *)msg0, i3, (int)sizeof(unsigned
          char));
        for (i3 = 0; i3 < len; i3++) {
          msg0->data[i3] = 0;
        }

        memcpy(&msg0->data[0], ptr, len);
        if (1 > len) {
          flag = 0;
        } else {
          flag = len;
        }

        i3 = varargin_1->size[0] * varargin_1->size[1];
        varargin_1->size[0] = 1;
        varargin_1->size[1] = flag;
        emxEnsureCapacity((emxArray__common *)varargin_1, i3, (int)sizeof
                          (unsigned char));
        for (i3 = 0; i3 < flag; i3++) {
          varargin_1->data[varargin_1->size[0] * i3] = msg0->data[i3];
        }

        i3 = b_varargin_1->size[0] * b_varargin_1->size[1];
        b_varargin_1->size[0] = 1;
        b_varargin_1->size[1] = flag;
        emxEnsureCapacity((emxArray__common *)b_varargin_1, i3, (int)sizeof(char));
        for (i3 = 0; i3 < flag; i3++) {
          b_varargin_1->data[i3] = (signed char)varargin_1->data[i3];
        }

        b_m2c_error(b_varargin_1);
      }
    }

    emxFree_char_T(&b_varargin_1);
    emxFree_uint8_T(&varargin_1);
    emxFree_uint8_T(&msg0);
  } else {
    vals_data = 0UL;
  }

  mat->rowptr = rowptr_data;
  mat->colind = colind_data;
  mat->vals = vals_data;
  mat->type = 2;
  mat->dims[0] = m;
  mat->dims[1] = n;
  mat->nnz = nnz;
}

void cuCRSCreate_initialize(void)
{
}

void cuCRSCreate_terminate(void)
{
}
