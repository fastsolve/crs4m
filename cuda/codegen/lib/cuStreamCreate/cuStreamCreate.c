#include "cuStreamCreate.h"
#include "mspack.h"
#include "m2c.h"

static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
static void emxInit_uint8_T1(emxArray_uint8_T **pEmxArray, int numDimensions);
static void m2c_error(const emxArray_char_T *varargin_3);

static void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_uint8_T(&pStruct->data);
  emxFree_char_T(&pStruct->type);
}

static void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_uint8_T(&pStruct->data, 1);
  emxInit_char_T(&pStruct->type, 2);
}

static void emxInit_uint8_T1(emxArray_uint8_T **pEmxArray, int numDimensions)
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

static void m2c_error(const emxArray_char_T *varargin_3)
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

  M2C_error("CUDA:RuntimeError", "cudaStreamCreate returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

void cuStreamCreate(struct0_T *stm, int *errCode, boolean_T *toplevel)
{
  emxArray_uint8_T *data0;
  cudaStream_t t_stm;
  int sizepe;
  int i0;
  char t0_type[12];
  static const char cv0[12] = { 'c', 'u', 'd', 'a', 'S', 't', 'r', 'e', 'a', 'm',
    '_', 't' };

  int loop_ub;
  char * ptr;
  int i;
  emxArray_uint8_T *msg0;
  const char * b_ptr;
  int len;
  emxArray_uint8_T *varargin_1;
  emxArray_char_T *b_varargin_1;
  emxInit_uint8_T(&data0, 1);
  *errCode = cudaStreamCreate(&t_stm);
  *toplevel = true;
  sizepe = sizeof(cudaStream_t);
  i0 = data0->size[0];
  data0->size[0] = sizepe;
  emxEnsureCapacity((emxArray__common *)data0, i0, (int)sizeof(unsigned char));
  for (i0 = 0; i0 < 12; i0++) {
    t0_type[i0] = cv0[i0];
  }

  i0 = stm->data->size[0];
  stm->data->size[0] = data0->size[0];
  emxEnsureCapacity((emxArray__common *)stm->data, i0, (int)sizeof(unsigned char));
  loop_ub = data0->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    stm->data->data[i0] = data0->data[i0];
  }

  emxFree_uint8_T(&data0);
  i0 = stm->type->size[0] * stm->type->size[1];
  stm->type->size[0] = 1;
  stm->type->size[1] = 12;
  emxEnsureCapacity((emxArray__common *)stm->type, i0, (int)sizeof(char));
  for (i0 = 0; i0 < 12; i0++) {
    stm->type->data[i0] = t0_type[i0];
  }

  stm->nitems = 1;
  ptr = (char *)(&t_stm);
  for (i = 1; i <= sizepe; i++) {
    stm->data->data[i - 1] = *(ptr);
    ptr = M2C_OFFSET_PTR(ptr, 1);
  }

  if (*errCode != 0) {
    emxInit_uint8_T1(&msg0, 2);
    b_ptr = cudaGetErrorString(*errCode);
    len = strlen(b_ptr) + 1;
    i0 = msg0->size[0] * msg0->size[1];
    msg0->size[0] = 1;
    msg0->size[1] = len;
    emxEnsureCapacity((emxArray__common *)msg0, i0, (int)sizeof(unsigned char));
    for (i0 = 0; i0 < len; i0++) {
      msg0->data[i0] = 0;
    }

    memcpy(&msg0->data[0], b_ptr, len);
    if (1 > len) {
      loop_ub = 0;
    } else {
      loop_ub = len;
    }

    emxInit_uint8_T1(&varargin_1, 2);
    i0 = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = 1;
    varargin_1->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)varargin_1, i0, (int)sizeof(unsigned
      char));
    for (i0 = 0; i0 < loop_ub; i0++) {
      varargin_1->data[varargin_1->size[0] * i0] = msg0->data[i0];
    }

    emxFree_uint8_T(&msg0);
    emxInit_char_T(&b_varargin_1, 2);
    i0 = b_varargin_1->size[0] * b_varargin_1->size[1];
    b_varargin_1->size[0] = 1;
    b_varargin_1->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_varargin_1, i0, (int)sizeof(char));
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_varargin_1->data[i0] = (signed char)varargin_1->data[i0];
    }

    emxFree_uint8_T(&varargin_1);
    m2c_error(b_varargin_1);
    emxFree_char_T(&b_varargin_1);
  }
}

void cuStreamCreate_initialize(void)
{
}

void cuStreamCreate_terminate(void)
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
