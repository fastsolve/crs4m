#include "cuBCRSDestroy.h"
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

static int32_T cuVecDestroy(uint64_T *vec);
static void emxFree_char_T(emxArray_char_T **pEmxArray);
static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray);
static void emxInit_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions);
static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int32_T numDimensions);
static void m2c_error(const emxArray_char_T *varargin_3);
static int32_T cuVecDestroy(uint64_T *vec)
{
  int32_T errCode;
  void * output;
  emxArray_uint8_T *msg0;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  int32_T flag;
  const char * ptr;
  int32_T len;
  int32_T i1;
  output = *(void **)(vec);
  errCode = cudaFree(output);
  emxInit_uint8_T(&msg0, 2);
  emxInit_uint8_T(&varargin_1, 2);
  emxInit_char_T(&b_varargin_1, 2);
  if (errCode != 0) {
    flag = (M2C_DEBUG);
    if (flag != 0) {
      ptr = cudaGetErrorString(errCode);
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

      m2c_error(b_varargin_1);
    }
  }

  emxFree_char_T(&b_varargin_1);
  emxFree_uint8_T(&varargin_1);
  emxFree_uint8_T(&msg0);
  flag = (M2C_DEBUG);
  if (flag != 0) {
    *vec = 0UL;
  }

  return errCode;
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

static void m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int32_T i0;
  int32_T loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i0 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i0, (int32_T)sizeof(char_T));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_varargin_3->data[i0] = varargin_3->data[i0];
  }

  M2C_error("CUDA:RuntimeError", "cudaFree returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

int32_T cuBCRSDestroy(struct0_T *mat)
{
  int32_T errCode;
  errCode = cuVecDestroy(&mat->vals);
  if (!(errCode != 0)) {
    errCode = cuVecDestroy(&mat->colind);
  }

  if (!(errCode != 0)) {
    errCode = cuVecDestroy(&mat->rowptr);
  }

  return errCode;
}

void cuBCRSDestroy_initialize(void)
{
}

void cuBCRSDestroy_terminate(void)
{
}
