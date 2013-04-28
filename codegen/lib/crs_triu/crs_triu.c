#include "crs_triu.h"
#include "stdio.h"
#ifdef BUILD_MEX

#include "mex.h"
#define malloc                         mxMalloc
#define calloc                         mxCalloc
#define realloc                        mxRealloc
#define free                           mxFree
#define emlrtIsMATLABThread(s)         1
#else
#define emlrtIsMATLABThread(s)         0
#define mexErrMsgIdAndTxt(a,b)
#endif

#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray__common
#define typedef_emxArray__common

typedef struct emxArray__common emxArray__common;

#endif

static void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
                     emxArray_real_T *val);
static void emxCopyStruct_struct_T(struct_T *dst, const struct_T *src);
static void emxCopy_int32_T(emxArray_int32_T **dst, emxArray_int32_T * const
  *src);
static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
static void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
                     emxArray_real_T *val)
{
  int32_T i2;
  int32_T i;
  emxArray_real_T *buf_val;
  emxArray_int32_T *buf_indx;
  boolean_T ascend;
  int32_T j;
  int32_T n;
  boolean_T exitg3;
  uint32_T ind;
  int32_T l;
  int32_T exitg1;
  boolean_T guard1 = FALSE;
  int32_T r0;
  real_T t0;
  int32_T exitg2;
  int32_T b_i;
  boolean_T guard2 = FALSE;
  i2 = row_ptr->size[0] - 1;
  i = 1;
  emxInit_real_T(&buf_val, 1);
  emxInit_int32_T(&buf_indx, 1);
  while (i <= i2) {
    ascend = TRUE;
    j = row_ptr->data[i - 1];
    n = row_ptr->data[i] - 1;
    exitg3 = FALSE;
    while ((exitg3 == FALSE) && (j + 1 <= n)) {
      if (col_ind->data[j] < col_ind->data[j - 1]) {
        ascend = FALSE;
        exitg3 = TRUE;
      } else {
        j++;
      }
    }

    if (!ascend) {
      n = buf_indx->size[0];
      buf_indx->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)buf_indx, n, (int32_T)sizeof(int32_T));
      n = buf_val->size[0];
      buf_val->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)buf_val, n, (int32_T)sizeof(real_T));
      ind = 1U;
      n = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= n; j++) {
        buf_indx->data[(int32_T)ind - 1] = col_ind->data[j - 1];
        buf_val->data[(int32_T)ind - 1] = val->data[j - 1];
        ind++;
      }

      n = buf_indx->size[0];
      if (n <= 1) {
      } else {
        l = (int32_T)((uint32_T)n >> 1U);
        do {
          exitg1 = 0;
          guard1 = FALSE;
          if (l + 1 <= 1) {
            r0 = buf_indx->data[n - 1];
            t0 = buf_val->data[n - 1];
            buf_indx->data[n - 1] = buf_indx->data[0];
            buf_val->data[n - 1] = buf_val->data[0];
            n--;
            if (n == 1) {
              exitg1 = 1;
            } else {
              guard1 = TRUE;
            }
          } else {
            l--;
            r0 = buf_indx->data[l];
            t0 = buf_val->data[l];
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            j = l;
            do {
              exitg2 = 0;
              b_i = j;
              j = ((j + 1) << 1) - 1;
              ascend = FALSE;
              guard2 = FALSE;
              if (j + 1 >= n) {
                if (j + 1 == n) {
                  ascend = TRUE;
                  guard2 = TRUE;
                } else if (j + 1 > n) {
                  exitg2 = 1;
                } else {
                  guard2 = TRUE;
                }
              } else {
                guard2 = TRUE;
              }

              if (guard2 == TRUE) {
                if ((!ascend) && (buf_indx->data[j] < buf_indx->data[j + 1])) {
                  j++;
                }

                if (r0 >= buf_indx->data[j]) {
                  exitg2 = 1;
                } else {
                  buf_indx->data[b_i] = buf_indx->data[j];
                  buf_val->data[b_i] = buf_val->data[j];
                }
              }
            } while (exitg2 == 0);

            buf_indx->data[b_i] = r0;
            buf_val->data[b_i] = t0;
          }
        } while (exitg1 == 0);

        buf_indx->data[0] = r0;
        buf_val->data[0] = t0;
      }

      ind = 1U;
      n = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= n; j++) {
        col_ind->data[j - 1] = buf_indx->data[(int32_T)ind - 1];
        val->data[j - 1] = buf_val->data[(int32_T)ind - 1];
        ind++;
      }
    }

    i++;
  }

  emxFree_int32_T(&buf_indx);
  emxFree_real_T(&buf_val);
}

static void emxCopyStruct_struct_T(struct_T *dst, const struct_T *src)
{
  emxCopy_int32_T(&dst->row_ptr, &src->row_ptr);
  emxCopy_int32_T(&dst->col_ind, &src->col_ind);
  emxCopy_real_T(&dst->val, &src->val);
  dst->nrows = src->nrows;
  dst->ncols = src->ncols;
}

static void emxCopy_int32_T(emxArray_int32_T **dst, emxArray_int32_T * const
  *src)
{
  int32_T numElDst;
  int32_T numElSrc;
  int32_T i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity((emxArray__common *)*dst, numElDst, (int32_T)sizeof(int32_T));
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src)
{
  int32_T numElDst;
  int32_T numElSrc;
  int32_T i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity((emxArray__common *)*dst, numElDst, (int32_T)sizeof(real_T));
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize)
{
  int32_T newNumel;
  int32_T i;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    if (emxArray->allocatedSize==0)
      i = newNumel;
    else {
      i = emxArray->allocatedSize;
      if (i < 16) {
        i = 16;
      }

      while (i < newNumel) {
        i <<= 1;
      }
    }

    newData = calloc((uint32_T)i, (uint32_T)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (uint32_T)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = TRUE;
  }
}


static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions)
{
  emxArray_int32_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int32_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

void crs_triu(const struct_T *A, struct_T *U)
{
  int32_T i0;
  int32_T offset;
  int32_T start;
  int32_T i;
  emxArray_int32_T *b_A;
  emxArray_real_T *c_A;
  emxCopyStruct_struct_T(U, A);
  i0 = U->col_ind->size[0];
  U->col_ind->size[0] = A->col_ind->size[0];
  emxEnsureCapacity((emxArray__common *)U->col_ind, i0, (int32_T)sizeof(int32_T));
  offset = A->col_ind->size[0];
  for (i0 = 0; i0 < offset; i0++) {
    U->col_ind->data[i0] = A->col_ind->data[i0];
  }

  i0 = U->val->size[0];
  U->val->size[0] = A->val->size[0];
  emxEnsureCapacity((emxArray__common *)U->val, i0, (int32_T)sizeof(real_T));
  offset = A->val->size[0];
  for (i0 = 0; i0 < offset; i0++) {
    U->val->data[i0] = A->val->data[i0];
  }

  crs_sort(A->row_ptr, U->col_ind, U->val);
  offset = 0;
  start = 0;
  for (i = 1; i <= A->nrows; i++) {
    i0 = U->row_ptr->data[i] - 1;
    while (start + 1 <= i0) {
      if (U->col_ind->data[start] < i) {
        offset++;
      } else {
        if (offset != 0) {
          U->col_ind->data[start - offset] = U->col_ind->data[start];
          U->val->data[start - offset] = U->val->data[start];
        }
      }

      start++;
    }

    start = U->row_ptr->data[i] - 1;
    U->row_ptr->data[i] -= offset;
  }

  if (offset != 0) {
    emxInit_int32_T(&b_A, 1);
    start = U->col_ind->size[0] - offset;
    i0 = b_A->size[0];
    b_A->size[0] = U->col_ind->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, i0, (int32_T)sizeof(int32_T));
    offset = U->col_ind->size[0];
    for (i0 = 0; i0 < offset; i0++) {
      b_A->data[i0] = U->col_ind->data[i0];
    }

    if (start < 1) {
      i0 = U->col_ind->size[0];
      U->col_ind->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)U->col_ind, i0, (int32_T)sizeof
                        (int32_T));
    } else {
      i0 = U->col_ind->size[0];
      U->col_ind->size[0] = start;
      emxEnsureCapacity((emxArray__common *)U->col_ind, i0, (int32_T)sizeof
                        (int32_T));
      for (i = 0; i < start; i++) {
        U->col_ind->data[i] = b_A->data[i];
      }
    }

    emxFree_int32_T(&b_A);
    emxInit_real_T(&c_A, 1);
    i0 = c_A->size[0];
    c_A->size[0] = U->val->size[0];
    emxEnsureCapacity((emxArray__common *)c_A, i0, (int32_T)sizeof(real_T));
    offset = U->val->size[0];
    for (i0 = 0; i0 < offset; i0++) {
      c_A->data[i0] = U->val->data[i0];
    }

    if (start < 1) {
      i0 = U->val->size[0];
      U->val->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)U->val, i0, (int32_T)sizeof(real_T));
    } else {
      i0 = U->val->size[0];
      U->val->size[0] = start;
      emxEnsureCapacity((emxArray__common *)U->val, i0, (int32_T)sizeof(real_T));
      for (i = 0; i < start; i++) {
        U->val->data[i] = c_A->data[i];
      }
    }

    emxFree_real_T(&c_A);
  }
}

void crs_triu1(const struct_T *A, int32_T k, struct_T *U)
{
  int32_T i1;
  int32_T offset;
  int32_T start;
  int32_T i;
  emxArray_int32_T *b_A;
  emxArray_real_T *c_A;
  emxCopyStruct_struct_T(U, A);
  i1 = U->col_ind->size[0];
  U->col_ind->size[0] = A->col_ind->size[0];
  emxEnsureCapacity((emxArray__common *)U->col_ind, i1, (int32_T)sizeof(int32_T));
  offset = A->col_ind->size[0];
  for (i1 = 0; i1 < offset; i1++) {
    U->col_ind->data[i1] = A->col_ind->data[i1];
  }

  i1 = U->val->size[0];
  U->val->size[0] = A->val->size[0];
  emxEnsureCapacity((emxArray__common *)U->val, i1, (int32_T)sizeof(real_T));
  offset = A->val->size[0];
  for (i1 = 0; i1 < offset; i1++) {
    U->val->data[i1] = A->val->data[i1];
  }

  crs_sort(A->row_ptr, U->col_ind, U->val);
  offset = 0;
  start = 0;
  for (i = 1; i <= A->nrows; i++) {
    i1 = U->row_ptr->data[i] - 1;
    while (start + 1 <= i1) {
      if (U->col_ind->data[start] < i + k) {
        offset++;
      } else {
        if (offset != 0) {
          U->col_ind->data[start - offset] = U->col_ind->data[start];
          U->val->data[start - offset] = U->val->data[start];
        }
      }

      start++;
    }

    start = U->row_ptr->data[i] - 1;
    U->row_ptr->data[i] -= offset;
  }

  if (offset != 0) {
    emxInit_int32_T(&b_A, 1);
    start = U->col_ind->size[0] - offset;
    i1 = b_A->size[0];
    b_A->size[0] = U->col_ind->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, i1, (int32_T)sizeof(int32_T));
    offset = U->col_ind->size[0];
    for (i1 = 0; i1 < offset; i1++) {
      b_A->data[i1] = U->col_ind->data[i1];
    }

    if (start < 1) {
      i1 = U->col_ind->size[0];
      U->col_ind->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)U->col_ind, i1, (int32_T)sizeof
                        (int32_T));
    } else {
      i1 = U->col_ind->size[0];
      U->col_ind->size[0] = start;
      emxEnsureCapacity((emxArray__common *)U->col_ind, i1, (int32_T)sizeof
                        (int32_T));
      for (i = 0; i < start; i++) {
        U->col_ind->data[i] = b_A->data[i];
      }
    }

    emxFree_int32_T(&b_A);
    emxInit_real_T(&c_A, 1);
    i1 = c_A->size[0];
    c_A->size[0] = U->val->size[0];
    emxEnsureCapacity((emxArray__common *)c_A, i1, (int32_T)sizeof(real_T));
    offset = U->val->size[0];
    for (i1 = 0; i1 < offset; i1++) {
      c_A->data[i1] = U->val->data[i1];
    }

    if (start < 1) {
      i1 = U->val->size[0];
      U->val->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)U->val, i1, (int32_T)sizeof(real_T));
    } else {
      i1 = U->val->size[0];
      U->val->size[0] = start;
      emxEnsureCapacity((emxArray__common *)U->val, i1, (int32_T)sizeof(real_T));
      for (i = 0; i < start; i++) {
        U->val->data[i] = c_A->data[i];
      }
    }

    emxFree_real_T(&c_A);
  }
}

void crs_triu_initialize(void)
{
}

void crs_triu_terminate(void)
{
}

emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions,
  int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_int32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_int32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols)
{
  emxArray_int32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_int32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_int32_T(emxArray_int32_T *emxArray)
{
  emxFree_int32_T(&emxArray);
}

void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}
