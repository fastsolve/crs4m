#include "cuSparseCreate.h"
#include "mspack.h"
#include "m2c.h"

static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
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

  M2C_error("CUDA:RuntimeError", "cusparseCreate returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

void cuSparseCreate(struct0_T *hdl, int32_T *errCode, boolean_T *toplevel)
{
  emxArray_uint8_T *data0;
  cusparseHandle_t t_hdl;
  int32_T sizepe;
  int32_T varargin_1;
  char_T t0_type[16];
  static const char_T cv0[16] = { 'c', 'u', 's', 'p', 'a', 'r', 's', 'e', 'H',
    'a', 'n', 'd', 'l', 'e', '_', 't' };

  int32_T loop_ub;
  char * ptr;
  int32_T i;
  int32_T varargin_3;
  int32_T varargin_4;
  int32_T varargin_5;
  int32_T varargin_6;
  int32_T varargin_7;
  int32_T varargin_8;
  int32_T varargin_9;
  boolean_T result;
  emxArray_char_T *cstr;
  static const char_T cv1[14] = { 'U', 'n', 'k', 'n', 'o', 'w', 'n', ' ', 'e',
    'r', 'r', 'o', 'r', '\x00' };

  static const char_T cv2[24] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'S', 'U', 'C', 'C', 'E', 'S', 'S', '\x00'
  };

  static const char_T cv3[32] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'N', 'O', 'T', '_', 'I', 'N', 'I', 'T',
    'I', 'A', 'L', 'I', 'Z', 'E', 'D', '\x00' };

  static const char_T cv4[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'A', 'L', 'L', 'O', 'C', '_', 'F', 'A',
    'I', 'L', 'E', 'D', '\x00' };

  static const char_T cv5[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'V', 'A', 'L', 'I', 'D', '_',
    'V', 'A', 'L', 'U', 'E', '\x00' };

  static const char_T cv6[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'A', 'R', 'C', 'H', '_', 'M', 'I', 'S',
    'M', 'A', 'T', 'C', 'H', '\x00' };

  static const char_T cv7[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'P', 'P', 'I', 'N', 'G', '_',
    'E', 'R', 'R', 'O', 'R', '\x00' };

  static const char_T cv8[33] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'E', 'X', 'E', 'C', 'U', 'T', 'I', 'O',
    'N', '_', 'F', 'A', 'I', 'L', 'E', 'D', '\x00' };

  static const char_T cv9[31] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'T', 'E', 'R', 'N', 'A', 'L',
    '_', 'E', 'R', 'R', 'O', 'R', '\x00' };

  static const char_T cv10[42] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'T', 'R', 'I', 'X', '_', 'T',
    'Y', 'P', 'E', '_', 'N', 'O', 'T', '_', 'S', 'U', 'P', 'P', 'O', 'R', 'T',
    'E', 'D', '\x00' };

  emxInit_uint8_T(&data0, 1);
  *errCode = cusparseCreate(&t_hdl);
  *toplevel = true;
  sizepe = sizeof(cusparseHandle_t);
  varargin_1 = data0->size[0];
  data0->size[0] = sizepe;
  emxEnsureCapacity((emxArray__common *)data0, varargin_1, (int32_T)sizeof
                    (uint8_T));
  for (varargin_1 = 0; varargin_1 < 16; varargin_1++) {
    t0_type[varargin_1] = cv0[varargin_1];
  }

  varargin_1 = hdl->data->size[0];
  hdl->data->size[0] = data0->size[0];
  emxEnsureCapacity((emxArray__common *)hdl->data, varargin_1, (int32_T)sizeof
                    (uint8_T));
  loop_ub = data0->size[0];
  for (varargin_1 = 0; varargin_1 < loop_ub; varargin_1++) {
    hdl->data->data[varargin_1] = data0->data[varargin_1];
  }

  emxFree_uint8_T(&data0);
  varargin_1 = hdl->type->size[0] * hdl->type->size[1];
  hdl->type->size[0] = 1;
  hdl->type->size[1] = 16;
  emxEnsureCapacity((emxArray__common *)hdl->type, varargin_1, (int32_T)sizeof
                    (char_T));
  for (varargin_1 = 0; varargin_1 < 16; varargin_1++) {
    hdl->type->data[varargin_1] = t0_type[varargin_1];
  }

  hdl->nitems = 1;
  ptr = (char *)(&t_hdl);
  for (i = 1; i <= sizepe; i++) {
    hdl->data->data[i - 1] = *(ptr);
    ptr = M2C_OFFSET_PTR(ptr, 1);
  }

  if (*errCode != 0) {
    varargin_1 = (CUSPARSE_STATUS_SUCCESS);
    loop_ub = (CUSPARSE_STATUS_NOT_INITIALIZED);
    varargin_3 = (CUSPARSE_STATUS_ALLOC_FAILED);
    varargin_4 = (CUSPARSE_STATUS_INVALID_VALUE);
    varargin_5 = (CUSPARSE_STATUS_ARCH_MISMATCH);
    varargin_6 = (CUSPARSE_STATUS_MAPPING_ERROR);
    varargin_7 = (CUSPARSE_STATUS_EXECUTION_FAILED);
    varargin_8 = (CUSPARSE_STATUS_INTERNAL_ERROR);
    varargin_9 = (CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED);
    result = (varargin_1 == *errCode);
    if (result) {
      varargin_1 = 0;
    } else {
      result = (loop_ub == *errCode);
      if (result) {
        varargin_1 = 1;
      } else {
        result = (varargin_3 == *errCode);
        if (result) {
          varargin_1 = 2;
        } else {
          result = (varargin_4 == *errCode);
          if (result) {
            varargin_1 = 3;
          } else {
            result = (varargin_5 == *errCode);
            if (result) {
              varargin_1 = 4;
            } else {
              result = (varargin_6 == *errCode);
              if (result) {
                varargin_1 = 5;
              } else {
                result = (varargin_7 == *errCode);
                if (result) {
                  varargin_1 = 6;
                } else {
                  result = (varargin_8 == *errCode);
                  if (result) {
                    varargin_1 = 7;
                  } else {
                    result = (varargin_9 == *errCode);
                    if (result) {
                      varargin_1 = 8;
                    } else {
                      varargin_1 = -1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    emxInit_char_T(&cstr, 2);
    switch (varargin_1) {
     case 0:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 24;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 24; varargin_1++) {
        cstr->data[varargin_1] = cv2[varargin_1];
      }
      break;

     case 1:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 32;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 32; varargin_1++) {
        cstr->data[varargin_1] = cv3[varargin_1];
      }
      break;

     case 2:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 29;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 29; varargin_1++) {
        cstr->data[varargin_1] = cv4[varargin_1];
      }
      break;

     case 3:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 30; varargin_1++) {
        cstr->data[varargin_1] = cv5[varargin_1];
      }
      break;

     case 4:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 30; varargin_1++) {
        cstr->data[varargin_1] = cv6[varargin_1];
      }
      break;

     case 5:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 30;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 30; varargin_1++) {
        cstr->data[varargin_1] = cv7[varargin_1];
      }
      break;

     case 6:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 33;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 33; varargin_1++) {
        cstr->data[varargin_1] = cv8[varargin_1];
      }
      break;

     case 7:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 31;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 31; varargin_1++) {
        cstr->data[varargin_1] = cv9[varargin_1];
      }
      break;

     case 8:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 42;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 42; varargin_1++) {
        cstr->data[varargin_1] = cv10[varargin_1];
      }
      break;

     default:
      varargin_1 = cstr->size[0] * cstr->size[1];
      cstr->size[0] = 1;
      cstr->size[1] = 14;
      emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int32_T)sizeof
                        (char_T));
      for (varargin_1 = 0; varargin_1 < 14; varargin_1++) {
        cstr->data[varargin_1] = cv1[varargin_1];
      }
      break;
    }

    m2c_error(cstr);
    emxFree_char_T(&cstr);
  }
}

void cuSparseCreate_initialize(void)
{
}

void cuSparseCreate_terminate(void)
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
