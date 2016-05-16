#include "crs_tril.h"
#include "omp.h"
#include "m2c.h"

static void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
                     emxArray_real_T *val);
static void emxCopyStruct_struct0_T(struct0_T *dst, const struct0_T *src);
static void emxCopy_int32_T(emxArray_int32_T **dst, emxArray_int32_T * const
  *src);
static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src);
static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
static void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
                     emxArray_real_T *val)
{
  int i2;
  int i;
  emxArray_real_T *buf_val;
  emxArray_int32_T *buf_indx;
  emxArray_int32_T *A;
  emxArray_real_T *b_A;
  boolean_T ascend;
  int i3;
  int j;
  boolean_T exitg3;
  int loop_ub;
  unsigned int ind;
  int n;
  int l;
  int ir;
  int exitg1;
  boolean_T guard1 = false;
  int r0;
  double t0;
  int exitg2;
  int b_i;
  boolean_T gt;
  boolean_T guard2 = false;
  i2 = row_ptr->size[0] - 1;
  i = 1;
  emxInit_real_T(&buf_val, 1);
  emxInit_int32_T(&buf_indx, 1);
  emxInit_int32_T(&A, 1);
  emxInit_real_T(&b_A, 1);
  while (i <= i2) {
    ascend = true;
    i3 = row_ptr->data[i] - 1;
    j = row_ptr->data[i - 1];
    exitg3 = false;
    while ((!exitg3) && (j + 1 <= i3)) {
      if (col_ind->data[j] < col_ind->data[j - 1]) {
        ascend = false;
        exitg3 = true;
      } else {
        j++;
      }
    }

    if (!ascend) {
      i3 = A->size[0];
      A->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)A, i3, (int)sizeof(int));
      loop_ub = row_ptr->data[i] - row_ptr->data[i - 1];
      for (i3 = 0; i3 < loop_ub; i3++) {
        A->data[i3] = 0;
      }

      i3 = buf_indx->size[0];
      buf_indx->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)buf_indx, i3, (int)sizeof(int));
      i3 = b_A->size[0];
      b_A->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)b_A, i3, (int)sizeof(double));
      loop_ub = row_ptr->data[i] - row_ptr->data[i - 1];
      for (i3 = 0; i3 < loop_ub; i3++) {
        b_A->data[i3] = 0.0;
      }

      i3 = buf_val->size[0];
      buf_val->size[0] = b_A->size[0];
      emxEnsureCapacity((emxArray__common *)buf_val, i3, (int)sizeof(double));
      ind = 1U;
      i3 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i3; j++) {
        buf_indx->data[(int)ind - 1] = col_ind->data[j - 1];
        buf_val->data[(int)ind - 1] = val->data[j - 1];
        ind++;
      }

      n = buf_indx->size[0];
      if (n <= 1) {
      } else {
        l = (int)((unsigned int)n >> 1U);
        ir = n;
        do {
          exitg1 = 0;
          guard1 = false;
          if (l + 1 <= 1) {
            r0 = buf_indx->data[ir - 1];
            t0 = buf_val->data[ir - 1];
            buf_indx->data[ir - 1] = buf_indx->data[0];
            buf_val->data[ir - 1] = buf_val->data[0];
            ir--;
            if (ir == 1) {
              exitg1 = 1;
            } else {
              guard1 = true;
            }
          } else {
            l--;
            r0 = buf_indx->data[l];
            t0 = buf_val->data[l];
            guard1 = true;
          }

          if (guard1) {
            j = l;
            do {
              exitg2 = 0;
              b_i = j;
              j = ((j + 1) << 1) - 1;
              gt = false;
              guard2 = false;
              if (j + 1 >= ir) {
                if (j + 1 == ir) {
                  gt = true;
                  guard2 = true;
                } else if (j + 1 > ir) {
                  exitg2 = 1;
                } else {
                  guard2 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                if ((!gt) && (buf_indx->data[j] < buf_indx->data[j + 1])) {
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
      i3 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i3; j++) {
        col_ind->data[j - 1] = buf_indx->data[(int)ind - 1];
        val->data[j - 1] = buf_val->data[(int)ind - 1];
        ind++;
      }
    }

    i++;
  }

  emxFree_real_T(&b_A);
  emxFree_int32_T(&A);
  emxFree_int32_T(&buf_indx);
  emxFree_real_T(&buf_val);
}

static void emxCopyStruct_struct0_T(struct0_T *dst, const struct0_T *src)
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
  int numElDst;
  int numElSrc;
  int i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity((emxArray__common *)*dst, numElDst, (int)sizeof(int));
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src)
{
  int numElDst;
  int numElSrc;
  int i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity((emxArray__common *)*dst, numElDst, (int)sizeof(double));
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

static void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_int32_T(&pStruct->row_ptr);
  emxFree_int32_T(&pStruct->col_ind);
  emxFree_real_T(&pStruct->val);
}

static void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_int32_T(&pStruct->row_ptr, 1);
  emxInit_int32_T(&pStruct->col_ind, 1);
  emxInit_real_T(&pStruct->val, 1);
}

void crs_tril(const struct0_T *A, struct0_T *L)
{
  int offset;
  int start;
  int i;
  int i0;
  emxArray_int32_T *b_A;
  int j;
  int newlen;
  int loop_ub;
  emxArray_real_T *c_A;
  emxCopyStruct_struct0_T(L, A);
  crs_sort(L->row_ptr, L->col_ind, L->val);
  offset = 0;
  start = 0;
  for (i = 1; i <= L->nrows; i++) {
    i0 = L->row_ptr->data[i] - 1;
    for (j = start; j + 1 <= i0; j++) {
      if (L->col_ind->data[j] > i) {
        offset++;
      } else {
        if (offset != 0) {
          L->col_ind->data[j - offset] = L->col_ind->data[j];
          L->val->data[j - offset] = L->val->data[j];
        }
      }
    }

    start = L->row_ptr->data[i] - 1;
    L->row_ptr->data[i] -= offset;
  }

  if (offset != 0) {
    emxInit_int32_T(&b_A, 1);
    newlen = L->col_ind->size[0] - offset;
    i0 = b_A->size[0];
    b_A->size[0] = L->col_ind->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, i0, (int)sizeof(int));
    loop_ub = L->col_ind->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_A->data[i0] = L->col_ind->data[i0];
    }

    if (newlen < 1) {
      i0 = L->col_ind->size[0];
      L->col_ind->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)L->col_ind, i0, (int)sizeof(int));
    } else {
      i0 = L->col_ind->size[0];
      L->col_ind->size[0] = newlen;
      emxEnsureCapacity((emxArray__common *)L->col_ind, i0, (int)sizeof(int));
      for (i = 0; i < newlen; i++) {
        L->col_ind->data[i] = b_A->data[i];
      }
    }

    emxFree_int32_T(&b_A);
    emxInit_real_T(&c_A, 1);
    i0 = c_A->size[0];
    c_A->size[0] = L->val->size[0];
    emxEnsureCapacity((emxArray__common *)c_A, i0, (int)sizeof(double));
    loop_ub = L->val->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_A->data[i0] = L->val->data[i0];
    }

    if (newlen < 1) {
      i0 = L->val->size[0];
      L->val->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)L->val, i0, (int)sizeof(double));
    } else {
      i0 = L->val->size[0];
      L->val->size[0] = newlen;
      emxEnsureCapacity((emxArray__common *)L->val, i0, (int)sizeof(double));
      for (i = 0; i < newlen; i++) {
        L->val->data[i] = c_A->data[i];
      }
    }

    emxFree_real_T(&c_A);
  }
}

void crs_tril1(const struct0_T *A, int k, struct0_T *L)
{
  int offset;
  int start;
  int i;
  int i1;
  emxArray_int32_T *b_A;
  int j;
  int newlen;
  int loop_ub;
  emxArray_real_T *c_A;
  emxCopyStruct_struct0_T(L, A);
  crs_sort(L->row_ptr, L->col_ind, L->val);
  offset = 0;
  start = 0;
  for (i = 1; i <= L->nrows; i++) {
    i1 = L->row_ptr->data[i] - 1;
    for (j = start; j + 1 <= i1; j++) {
      if (L->col_ind->data[j] > i + k) {
        offset++;
      } else {
        if (offset != 0) {
          L->col_ind->data[j - offset] = L->col_ind->data[j];
          L->val->data[j - offset] = L->val->data[j];
        }
      }
    }

    start = L->row_ptr->data[i] - 1;
    L->row_ptr->data[i] -= offset;
  }

  if (offset != 0) {
    emxInit_int32_T(&b_A, 1);
    newlen = L->col_ind->size[0] - offset;
    i1 = b_A->size[0];
    b_A->size[0] = L->col_ind->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, i1, (int)sizeof(int));
    loop_ub = L->col_ind->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_A->data[i1] = L->col_ind->data[i1];
    }

    if (newlen < 1) {
      i1 = L->col_ind->size[0];
      L->col_ind->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)L->col_ind, i1, (int)sizeof(int));
    } else {
      i1 = L->col_ind->size[0];
      L->col_ind->size[0] = newlen;
      emxEnsureCapacity((emxArray__common *)L->col_ind, i1, (int)sizeof(int));
      for (i = 0; i < newlen; i++) {
        L->col_ind->data[i] = b_A->data[i];
      }
    }

    emxFree_int32_T(&b_A);
    emxInit_real_T(&c_A, 1);
    i1 = c_A->size[0];
    c_A->size[0] = L->val->size[0];
    emxEnsureCapacity((emxArray__common *)c_A, i1, (int)sizeof(double));
    loop_ub = L->val->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_A->data[i1] = L->val->data[i1];
    }

    if (newlen < 1) {
      i1 = L->val->size[0];
      L->val->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)L->val, i1, (int)sizeof(double));
    } else {
      i1 = L->val->size[0];
      L->val->size[0] = newlen;
      emxEnsureCapacity((emxArray__common *)L->val, i1, (int)sizeof(double));
      for (i = 0; i < newlen; i++) {
        L->val->data[i] = c_A->data[i];
      }
    }

    emxFree_real_T(&c_A);
  }
}

void crs_tril_initialize(void)
{
}

void crs_tril_terminate(void)
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
