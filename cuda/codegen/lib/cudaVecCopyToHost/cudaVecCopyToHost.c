#include "cudaVecCopyToHost.h"
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

static void b_m2c_error(void);
static void c_m2c_error(void);
static void d_m2c_error(const emxArray_char_T *varargin_3);
static void emxFree_char_T(emxArray_char_T **pEmxArray);
static void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);
static void m2c_error(void);
static void b_m2c_error(void)
{
  M2C_error("cudaVecCopyFromHost:SizeMismatch", "Target array is too small.");
}

static void c_m2c_error(void)
{
  M2C_error("cudaVecCopyFromHost:TypeMismatch", "Expected real numbers.");
}

static void d_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int i0;
  int loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i0 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i0, (int)sizeof(char));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_varargin_3->data[i0] = varargin_3->data[i0];
  }

  M2C_error("CUDA:RuntimeError", "cublasGetVector returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
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

static void m2c_error(void)
{
  M2C_error("cudaVecCopyFromHost:TypeMismatch",
            "Real and complex numbers mismatch.");
}

void cudaVecCopyToHost(const struct0_T *cuVec, const emxArray_real_T *vec, int
  *errCode, boolean_T *toplevel)
{
  int i;
  int sizepe;
  unsigned int data[2];
  void * output;
  char * ptr;
  int varargin_2;
  int varargin_3;
  int varargin_4;
  int varargin_5;
  int varargin_6;
  int varargin_7;
  int varargin_8;
  emxArray_char_T *cstr;
  static const char cv0[14] = { 'U', 'n', 'k', 'n', 'o', 'w', 'n', ' ', 'e', 'r',
    'r', 'o', 'r', '\x00' };

  static const char cv1[22] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'S', 'U', 'C', 'C', 'E', 'S', 'S', '\x00' };

  static const char cv2[30] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'N', 'O', 'T', '_', 'I', 'N', 'I', 'T', 'I', 'A', 'L',
    'I', 'Z', 'E', 'D', '\x00' };

  static const char cv3[27] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'A', 'L', 'L', 'O', 'C', '_', 'F', 'A', 'I', 'L', 'E',
    'D', '\x00' };

  static const char cv4[28] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'I', 'N', 'V', 'A', 'L', 'I', 'D', '_', 'V', 'A', 'L',
    'U', 'E', '\x00' };

  static const char cv5[28] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'A', 'R', 'C', 'H', '_', 'M', 'I', 'S', 'M', 'A', 'T',
    'C', 'H', '\x00' };

  static const char cv6[28] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'M', 'A', 'P', 'P', 'I', 'N', 'G', '_', 'E', 'R', 'R',
    'O', 'R', '\x00' };

  static const char cv7[31] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'E', 'X', 'E', 'C', 'U', 'T', 'I', 'O', 'N', '_', 'F',
    'A', 'I', 'L', 'E', 'D', '\x00' };

  static const char cv8[29] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'I', 'N', 'T', 'E', 'R', 'N', 'A', 'L', '_', 'E', 'R',
    'R', 'O', 'R', '\x00' };

  if ((cuVec->type != 2) && (cuVec->type != 1)) {
    m2c_error();
  } else {
    i = M2C_INTDIV(vec->size[0], 1);
    if (cuVec->len > i) {
      b_m2c_error();
    }
  }

  if ((cuVec->type == 2) || (cuVec->type == 3)) {
    sizepe = 8;
  } else if (cuVec->type == 1) {
    sizepe = 4;
  } else if (cuVec->type == 4) {
    sizepe = 16;
  } else {
    sizepe = 0;
    c_m2c_error();
  }

  for (i = 0; i < 2; i++) {
    data[i] = cuVec->data[i];
  }

  output = *(void **)(data);
  ptr = (char *)(&vec->data[0]);
  *errCode = cublasGetVector(cuVec->len, sizepe, output, 1, ptr, 1);
  if (*errCode != 0) {
    i = (CUBLAS_STATUS_SUCCESS);
    varargin_2 = (CUBLAS_STATUS_NOT_INITIALIZED);
    varargin_3 = (CUBLAS_STATUS_ALLOC_FAILED);
    varargin_4 = (CUBLAS_STATUS_INVALID_VALUE);
    varargin_5 = (CUBLAS_STATUS_ARCH_MISMATCH);
    varargin_6 = (CUBLAS_STATUS_MAPPING_ERROR);
    varargin_7 = (CUBLAS_STATUS_EXECUTION_FAILED);
    varargin_8 = (CUBLAS_STATUS_INTERNAL_ERROR);
    if (i == *errCode) {
      i = 0;
    } else if (varargin_2 == *errCode) {
      i = 1;
    } else if (varargin_3 == *errCode) {
      i = 2;
    } else if (varargin_4 == *errCode) {
      i = 3;
    } else if (varargin_5 == *errCode) {
      i = 4;
    } else if (varargin_6 == *errCode) {
      i = 5;
    } else if (varargin_7 == *errCode) {
      i = 6;
    } else if (varargin_8 == *errCode) {
      i = 7;
    } else {
      i = -1;
    }

    emxInit_char_T(&cstr, 2);
    switch (i) {
     case 0:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 22;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 22; i++) {
        cstr->data[i] = cv1[i];
      }
      break;

     case 1:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 30; i++) {
        cstr->data[i] = cv2[i];
      }
      break;

     case 2:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 27;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 27; i++) {
        cstr->data[i] = cv3[i];
      }
      break;

     case 3:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 28;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 28; i++) {
        cstr->data[i] = cv4[i];
      }
      break;

     case 4:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 28;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 28; i++) {
        cstr->data[i] = cv5[i];
      }
      break;

     case 5:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 28;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 28; i++) {
        cstr->data[i] = cv6[i];
      }
      break;

     case 6:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 31;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 31; i++) {
        cstr->data[i] = cv7[i];
      }
      break;

     case 7:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 29;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 29; i++) {
        cstr->data[i] = cv8[i];
      }
      break;

     default:
      i = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 14;
      emxEnsureCapacity((emxArray__common *)cstr, i, (int)sizeof(char));
      for (i = 0; i < 14; i++) {
        cstr->data[i] = cv0[i];
      }
      break;
    }

    d_m2c_error(cstr);
    emxFree_char_T(&cstr);
  }

  *toplevel = true;
}

void cudaVecCopyToHost_initialize(void)
{
}

void cudaVecCopyToHost_terminate(void)
{
}

