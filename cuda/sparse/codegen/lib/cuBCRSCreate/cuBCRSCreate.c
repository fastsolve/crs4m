#include "cuBCRSCreate.h"
#include "mspack.h"
#include "m2c.h"

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

static void b_cuVecCreate(int32_T n, uint64_T *vec_data, int32_T *vec_type,
  int32_T *vec_len, int32_T *errCode);
static void b_m2c_error(const emxArray_char_T *varargin_3);
static void cuVecCreate(int32_T n, uint64_T *vec_data, int32_T *vec_type,
  int32_T *vec_len, int32_T *errCode);
static void emxFree_char_T(emxArray_char_T **pEmxArray);
static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray);
static void emxInit_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions);
static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int32_T numDimensions);
static void m2c_error(int32_T varargin_3);
static void b_cuVecCreate(int32_T n, uint64_T *vec_data, int32_T *vec_type,
  int32_T *vec_len, int32_T *errCode)
{
  void * t_vec;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int32_T flag;
  const char * ptr;
  int32_T len;
  int32_T i3;
  *errCode = cudaMalloc(&t_vec, n << 3);
  *vec_type = 2;
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
      i3 = msg0->size[0] * msg0->size[1];
      msg0->size[0] = 1;
      msg0->size[1] = len;
      emxEnsureCapacity((emxArray__common *)msg0, i3, (int32_T)sizeof(uint8_T));
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
      emxEnsureCapacity((emxArray__common *)varargin_1, i3, (int32_T)sizeof
                        (uint8_T));
      for (i3 = 0; i3 < flag; i3++) {
        varargin_1->data[varargin_1->size[0] * i3] = msg0->data[i3];
      }

      i3 = b_varargin_1->size[0] * b_varargin_1->size[1];
      b_varargin_1->size[0] = 1;
      b_varargin_1->size[1] = flag;
      emxEnsureCapacity((emxArray__common *)b_varargin_1, i3, (int32_T)sizeof
                        (char_T));
      for (i3 = 0; i3 < flag; i3++) {
        b_varargin_1->data[i3] = (int8_T)varargin_1->data[i3];
      }

      b_m2c_error(b_varargin_1);
    }
  }

  emxFree_char_T(&b_varargin_1);
  emxFree_uint8_T(&varargin_1);
  emxFree_uint8_T(&msg0);
}

static void b_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int32_T i2;
  int32_T loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i2 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i2, (int32_T)sizeof(char_T));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    b_varargin_3->data[i2] = varargin_3->data[i2];
  }

  M2C_error("CUDA:RuntimeError", "cudaMalloc returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

static void cuVecCreate(int32_T n, uint64_T *vec_data, int32_T *vec_type,
  int32_T *vec_len, int32_T *errCode)
{
  void * t_vec;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int32_T flag;
  const char * ptr;
  int32_T len;
  int32_T i1;
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
      emxEnsureCapacity((emxArray__common *)msg0, i1, (int32_T)sizeof(uint8_T));
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
      emxEnsureCapacity((emxArray__common *)varargin_1, i1, (int32_T)sizeof
                        (uint8_T));
      for (i1 = 0; i1 < flag; i1++) {
        varargin_1->data[varargin_1->size[0] * i1] = msg0->data[i1];
      }

      i1 = b_varargin_1->size[0] * b_varargin_1->size[1];
      b_varargin_1->size[0] = 1;
      b_varargin_1->size[1] = flag;
      emxEnsureCapacity((emxArray__common *)b_varargin_1, i1, (int32_T)sizeof
                        (char_T));
      for (i1 = 0; i1 < flag; i1++) {
        b_varargin_1->data[i1] = (int8_T)varargin_1->data[i1];
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
    if (((*pEmxArray)->data != (char_T *)NULL) && (*pEmxArray)->canFreeData) {
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
    if (((*pEmxArray)->data != (uint8_T *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_uint8_T *)NULL;
  }
}

static void emxInit_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions)
{
  emxArray_char_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_char_T *)malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int32_T numDimensions)
{
  emxArray_uint8_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_uint8_T *)malloc(sizeof(emxArray_uint8_T));
  emxArray = *pEmxArray;
  emxArray->data = (uint8_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void m2c_error(int32_T varargin_3)
{
  M2C_error("mspGetSizePerElement:WrongType", "Unknow data type %d.\n",
            varargin_3);
}

void cuBCRSCreate(int32_T mb, int32_T nb, int32_T nnzb, int32_T blkdim, int32_T
                  type, MSP_CuBCRS *mat, int32_T *errCode)
{
  uint64_T rowptr_data;
  int32_T expl_temp;
  int32_T b_expl_temp;
  uint64_T colind_data;
  uint64_T vals_data;
  int32_T sizepe;
  void * t_vec;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int32_T flag;
  const char * ptr;
  int32_T len;
  int32_T i0;
  cuVecCreate(mb + 1, &rowptr_data, &expl_temp, &b_expl_temp, errCode);
  if (!(*errCode != 0)) {
    cuVecCreate(nnzb, &colind_data, &expl_temp, &b_expl_temp, errCode);
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

    *errCode = cudaMalloc(&t_vec, nnzb * blkdim * blkdim * sizepe);
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
        emxEnsureCapacity((emxArray__common *)msg0, i0, (int32_T)sizeof(uint8_T));
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
        emxEnsureCapacity((emxArray__common *)varargin_1, i0, (int32_T)sizeof
                          (uint8_T));
        for (i0 = 0; i0 < flag; i0++) {
          varargin_1->data[varargin_1->size[0] * i0] = msg0->data[i0];
        }

        i0 = b_varargin_1->size[0] * b_varargin_1->size[1];
        b_varargin_1->size[0] = 1;
        b_varargin_1->size[1] = flag;
        emxEnsureCapacity((emxArray__common *)b_varargin_1, i0, (int32_T)sizeof
                          (char_T));
        for (i0 = 0; i0 < flag; i0++) {
          b_varargin_1->data[i0] = (int8_T)varargin_1->data[i0];
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
  mat->dimsb[0] = mb;
  mat->dimsb[1] = nb;
  mat->nnzb = nnzb;
  mat->blkdim = blkdim;
}

void cuBCRSCreate_3args(int32_T mb, int32_T nb, int32_T nnzb, MSP_CuBCRS *mat,
  int32_T *errCode)
{
  uint64_T rowptr_data;
  int32_T expl_temp;
  int32_T b_expl_temp;
  uint64_T colind_data;
  uint64_T vals_data;
  cuVecCreate(mb + 1, &rowptr_data, &expl_temp, &b_expl_temp, errCode);
  if (!(*errCode != 0)) {
    cuVecCreate(nnzb, &colind_data, &expl_temp, &b_expl_temp, errCode);
  } else {
    colind_data = 0UL;
  }

  if (!(*errCode != 0)) {
    b_cuVecCreate(nnzb, &vals_data, &expl_temp, &b_expl_temp, errCode);
  } else {
    vals_data = 0UL;
  }

  mat->rowptr = rowptr_data;
  mat->colind = colind_data;
  mat->vals = vals_data;
  mat->type = 2;
  mat->dimsb[0] = mb;
  mat->dimsb[1] = nb;
  mat->nnzb = nnzb;
  mat->blkdim = 1;
}

void cuBCRSCreate_4args(int32_T mb, int32_T nb, int32_T nnzb, int32_T blkdim,
  MSP_CuBCRS *mat, int32_T *errCode)
{
  uint64_T rowptr_data;
  int32_T expl_temp;
  int32_T b_expl_temp;
  uint64_T colind_data;
  uint64_T vals_data;
  cuVecCreate(mb + 1, &rowptr_data, &expl_temp, &b_expl_temp, errCode);
  if (!(*errCode != 0)) {
    cuVecCreate(nnzb, &colind_data, &expl_temp, &b_expl_temp, errCode);
  } else {
    colind_data = 0UL;
  }

  if (!(*errCode != 0)) {
    b_cuVecCreate(nnzb * blkdim * blkdim, &vals_data, &expl_temp, &b_expl_temp,
                  errCode);
  } else {
    vals_data = 0UL;
  }

  mat->rowptr = rowptr_data;
  mat->colind = colind_data;
  mat->vals = vals_data;
  mat->type = 2;
  mat->dimsb[0] = mb;
  mat->dimsb[1] = nb;
  mat->nnzb = nnzb;
  mat->blkdim = blkdim;
}

void cuBCRSCreate_initialize(void)
{
}

void cuBCRSCreate_terminate(void)
{
}
