#include "cuSparseDestroy.h"
#include "mspack.h"
#include "m2c.h"

static void b_m2c_error(const emxArray_char_T *varargin_3);
static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
static void m2c_error(const emxArray_char_T *varargin_3);
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

  M2C_error("CUDA:RuntimeError", "cusparseDestroy returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

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

static void m2c_error(const emxArray_char_T *varargin_3)
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

  M2C_error("m2c_opaque_obj:WrongInput",
            "Incorrect data type %s. Expected cusparseHandle_t.\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

void cuSparseDestroy(const struct0_T *hdl, int *errCode, boolean_T *toplevel)
{
  boolean_T p;
  boolean_T b_p;
  int k;
  int exitg2;
  int varargin_2;
  boolean_T exitg1;
  emxArray_char_T *b_hdl;
  static const char cv0[16] = { 'c', 'u', 's', 'p', 'a', 'r', 's', 'e', 'H', 'a',
    'n', 'd', 'l', 'e', '_', 't' };

  emxArray_uint8_T *data;
  cusparseHandle_t c_hdl;
  int varargin_3;
  int varargin_4;
  int varargin_5;
  int varargin_6;
  int varargin_7;
  int varargin_8;
  int varargin_9;
  emxArray_char_T *cstr;
  static const char cv1[14] = { 'U', 'n', 'k', 'n', 'o', 'w', 'n', ' ', 'e', 'r',
    'r', 'o', 'r', '\x00' };

  static const char cv2[24] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'S', 'U', 'C', 'C', 'E', 'S', 'S', '\x00' };

  static const char cv3[32] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'N', 'O', 'T', '_', 'I', 'N', 'I', 'T', 'I',
    'A', 'L', 'I', 'Z', 'E', 'D', '\x00' };

  static const char cv4[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'A', 'L', 'L', 'O', 'C', '_', 'F', 'A', 'I',
    'L', 'E', 'D', '\x00' };

  static const char cv5[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'V', 'A', 'L', 'I', 'D', '_', 'V',
    'A', 'L', 'U', 'E', '\x00' };

  static const char cv6[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'A', 'R', 'C', 'H', '_', 'M', 'I', 'S', 'M',
    'A', 'T', 'C', 'H', '\x00' };

  static const char cv7[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'P', 'P', 'I', 'N', 'G', '_', 'E',
    'R', 'R', 'O', 'R', '\x00' };

  static const char cv8[33] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'E', 'X', 'E', 'C', 'U', 'T', 'I', 'O', 'N',
    '_', 'F', 'A', 'I', 'L', 'E', 'D', '\x00' };

  static const char cv9[31] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'T', 'E', 'R', 'N', 'A', 'L', '_',
    'E', 'R', 'R', 'O', 'R', '\x00' };

  static const char cv10[42] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'T', 'R', 'I', 'X', '_', 'T',
    'Y', 'P', 'E', '_', 'N', 'O', 'T', '_', 'S', 'U', 'P', 'P', 'O', 'R', 'T',
    'E', 'D', '\x00' };

  p = false;
  b_p = false;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      varargin_2 = hdl->type->size[k];
      if (varargin_2 != 15 * k + 1) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      b_p = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(hdl->type->size[1] == 0))) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 16)) {
      if (!(hdl->type->data[k] == cv0[k])) {
        b_p = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (!b_p) {
  } else {
    p = true;
  }

  if (!p) {
    emxInit_char_T(&b_hdl, 2);
    varargin_2 = b_hdl->size[0] * b_hdl->size[1];
    b_hdl->size[0] = 1;
    b_hdl->size[1] = hdl->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_hdl, varargin_2, (int)sizeof(char));
    k = hdl->type->size[1];
    for (varargin_2 = 0; varargin_2 < k; varargin_2++) {
      b_hdl->data[b_hdl->size[0] * varargin_2] = hdl->type->data[hdl->type->
        size[0] * varargin_2];
    }

    b_hdl->data[b_hdl->size[0] * hdl->type->size[1]] = '\x00';
    m2c_error(b_hdl);
    emxFree_char_T(&b_hdl);
  }

  emxInit_uint8_T(&data, 1);
  varargin_2 = data->size[0];
  data->size[0] = hdl->data->size[0];
  emxEnsureCapacity((emxArray__common *)data, varargin_2, (int)sizeof(unsigned
    char));
  k = hdl->data->size[0];
  for (varargin_2 = 0; varargin_2 < k; varargin_2++) {
    data->data[varargin_2] = hdl->data->data[varargin_2];
  }

  c_hdl = *(cusparseHandle_t*)(&data->data[0]);
  *errCode = cusparseDestroy(c_hdl);
  *toplevel = true;
  emxFree_uint8_T(&data);
  if (*errCode != 0) {
    k = (CUSPARSE_STATUS_SUCCESS);
    varargin_2 = (CUSPARSE_STATUS_NOT_INITIALIZED);
    varargin_3 = (CUSPARSE_STATUS_ALLOC_FAILED);
    varargin_4 = (CUSPARSE_STATUS_INVALID_VALUE);
    varargin_5 = (CUSPARSE_STATUS_ARCH_MISMATCH);
    varargin_6 = (CUSPARSE_STATUS_MAPPING_ERROR);
    varargin_7 = (CUSPARSE_STATUS_EXECUTION_FAILED);
    varargin_8 = (CUSPARSE_STATUS_INTERNAL_ERROR);
    varargin_9 = (CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED);
    if (k == *errCode) {
      k = 0;
    } else if (varargin_2 == *errCode) {
      k = 1;
    } else if (varargin_3 == *errCode) {
      k = 2;
    } else if (varargin_4 == *errCode) {
      k = 3;
    } else if (varargin_5 == *errCode) {
      k = 4;
    } else if (varargin_6 == *errCode) {
      k = 5;
    } else if (varargin_7 == *errCode) {
      k = 6;
    } else if (varargin_8 == *errCode) {
      k = 7;
    } else if (varargin_9 == *errCode) {
      k = 8;
    } else {
      k = -1;
    }

    emxInit_char_T(&cstr, 2);
    switch (k) {
     case 0:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 24;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 24; varargin_2++) {
        cstr->data[varargin_2] = cv2[varargin_2];
      }
      break;

     case 1:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 32;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 32; varargin_2++) {
        cstr->data[varargin_2] = cv3[varargin_2];
      }
      break;

     case 2:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 29;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 29; varargin_2++) {
        cstr->data[varargin_2] = cv4[varargin_2];
      }
      break;

     case 3:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 30; varargin_2++) {
        cstr->data[varargin_2] = cv5[varargin_2];
      }
      break;

     case 4:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 30; varargin_2++) {
        cstr->data[varargin_2] = cv6[varargin_2];
      }
      break;

     case 5:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 30; varargin_2++) {
        cstr->data[varargin_2] = cv7[varargin_2];
      }
      break;

     case 6:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 33;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 33; varargin_2++) {
        cstr->data[varargin_2] = cv8[varargin_2];
      }
      break;

     case 7:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 31;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 31; varargin_2++) {
        cstr->data[varargin_2] = cv9[varargin_2];
      }
      break;

     case 8:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 42;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 42; varargin_2++) {
        cstr->data[varargin_2] = cv10[varargin_2];
      }
      break;

     default:
      varargin_2 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 14;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_2, (int)sizeof(char));
      for (varargin_2 = 0; varargin_2 < 14; varargin_2++) {
        cstr->data[varargin_2] = cv1[varargin_2];
      }
      break;
    }

    b_m2c_error(cstr);
    emxFree_char_T(&cstr);
  }
}

void cuSparseDestroy_initialize(void)
{
}

void cuSparseDestroy_terminate(void)
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
