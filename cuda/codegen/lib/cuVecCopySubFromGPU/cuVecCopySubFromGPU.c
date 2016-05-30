#include "cuVecCopySubFromGPU.h"
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
  M2C_error("cuVecCopyToGPU:SizeMismatch", "Target array is too small.");
}

static void c_m2c_error(void)
{
  M2C_error("cuVecCopyToGPU:TypeMismatch", "Expected real numbers.");
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
  M2C_error("cuVecCopyToGPU:TypeMismatch", "Real and complex numbers mismatch.");
}

void cuVecCopySubFromGPU(int n, const struct0_T *cuVec, int incCuVec, const
  emxArray_real_T *vec, int istart, int inc, int *errCode, boolean_T *toplevel)
{
  int quotient;
  int sizepe;
  unsigned long data;
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

  *toplevel = true;
  if ((cuVec->type != 2) && (cuVec->type != 1)) {
    m2c_error();
  } else {
    quotient = M2C_INTDIV(vec->size[0], inc);
    if (n > quotient) {
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

  data = cuVec->data;
  output = *(void **)(&data);
  quotient = (istart - 1) * sizepe;
  ptr = (char *)(&vec->data[0]);
  if (quotient != 0) {
    ptr = M2C_OFFSET_PTR(ptr, quotient);
  }

  *errCode = cublasGetVector(n, sizepe, output, incCuVec, ptr, inc);
  if (*errCode != 0) {
    quotient = (CUBLAS_STATUS_SUCCESS);
    varargin_2 = (CUBLAS_STATUS_NOT_INITIALIZED);
    varargin_3 = (CUBLAS_STATUS_ALLOC_FAILED);
    varargin_4 = (CUBLAS_STATUS_INVALID_VALUE);
    varargin_5 = (CUBLAS_STATUS_ARCH_MISMATCH);
    varargin_6 = (CUBLAS_STATUS_MAPPING_ERROR);
    varargin_7 = (CUBLAS_STATUS_EXECUTION_FAILED);
    varargin_8 = (CUBLAS_STATUS_INTERNAL_ERROR);
    if (quotient == *errCode) {
      quotient = 0;
    } else if (varargin_2 == *errCode) {
      quotient = 1;
    } else if (varargin_3 == *errCode) {
      quotient = 2;
    } else if (varargin_4 == *errCode) {
      quotient = 3;
    } else if (varargin_5 == *errCode) {
      quotient = 4;
    } else if (varargin_6 == *errCode) {
      quotient = 5;
    } else if (varargin_7 == *errCode) {
      quotient = 6;
    } else if (varargin_8 == *errCode) {
      quotient = 7;
    } else {
      quotient = -1;
    }

    emxInit_char_T(&cstr, 2);
    switch (quotient) {
     case 0:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 22;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 22; quotient++) {
        cstr->data[quotient] = cv1[quotient];
      }
      break;

     case 1:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 30; quotient++) {
        cstr->data[quotient] = cv2[quotient];
      }
      break;

     case 2:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 27;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 27; quotient++) {
        cstr->data[quotient] = cv3[quotient];
      }
      break;

     case 3:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 28;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 28; quotient++) {
        cstr->data[quotient] = cv4[quotient];
      }
      break;

     case 4:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 28;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 28; quotient++) {
        cstr->data[quotient] = cv5[quotient];
      }
      break;

     case 5:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 28;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 28; quotient++) {
        cstr->data[quotient] = cv6[quotient];
      }
      break;

     case 6:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 31;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 31; quotient++) {
        cstr->data[quotient] = cv7[quotient];
      }
      break;

     case 7:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 29;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 29; quotient++) {
        cstr->data[quotient] = cv8[quotient];
      }
      break;

     default:
      quotient = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 14;
      emxEnsureCapacity((emxArray__common *)cstr, quotient, (int)sizeof(char));
      for (quotient = 0; quotient < 14; quotient++) {
        cstr->data[quotient] = cv0[quotient];
      }
      break;
    }

    d_m2c_error(cstr);
    emxFree_char_T(&cstr);
  }
}

void cuVecCopySubFromGPU_initialize(void)
{
}

void cuVecCopySubFromGPU_terminate(void)
{
}

