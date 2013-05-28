#include "crs_tril.h"
#include "palc.h"

static void crs_sort(const plcArray_int32_T *row_ptr, plcArray_int32_T *col_ind,
                     plcArray_real_T *val);
static void plcCopyStruct_struct_T(struct_T *dst, const struct_T *src);
static define_plcCopy(plcCopy_int32_T, int32_T);
static define_plcCopy(plcCopy_real_T, real_T);
static void crs_sort(const plcArray_int32_T *row_ptr, plcArray_int32_T *col_ind,
                     plcArray_real_T *val)
{
  int32_T i2;
  int32_T i;
  plcArray_real_T *buf_val;
  plcArray_int32_T *buf_indx;
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
  plcInit_real_T(&buf_val, 1);
  plcInit_int32_T(&buf_indx, 1);
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
      plcEnsureCapacity((plcArray__common *)buf_indx, n, (int32_T)sizeof(int32_T));
      n = buf_val->size[0];
      buf_val->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      plcEnsureCapacity((plcArray__common *)buf_val, n, (int32_T)sizeof(real_T));
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

  plcFree_int32_T(&buf_indx);
  plcFree_real_T(&buf_val);
}

static void plcCopyStruct_struct_T(struct_T *dst, const struct_T *src)
{
  plcCopy_int32_T(&dst->row_ptr, &src->row_ptr);
  plcCopy_int32_T(&dst->col_ind, &src->col_ind);
  plcCopy_real_T(&dst->val, &src->val);
  dst->nrows = src->nrows;
  dst->ncols = src->ncols;
}

void crs_tril(const struct_T *A, struct_T *L)
{
  int32_T i0;
  int32_T offset;
  int32_T start;
  int32_T i;
  plcArray_int32_T *b_A;
  plcArray_real_T *c_A;
  plcCopyStruct_struct_T(L, A);
  i0 = L->col_ind->size[0];
  L->col_ind->size[0] = A->col_ind->size[0];
  plcEnsureCapacity((plcArray__common *)L->col_ind, i0, (int32_T)sizeof(int32_T));
  offset = A->col_ind->size[0];
  for (i0 = 0; i0 < offset; i0++) {
    L->col_ind->data[i0] = A->col_ind->data[i0];
  }

  i0 = L->val->size[0];
  L->val->size[0] = A->val->size[0];
  plcEnsureCapacity((plcArray__common *)L->val, i0, (int32_T)sizeof(real_T));
  offset = A->val->size[0];
  for (i0 = 0; i0 < offset; i0++) {
    L->val->data[i0] = A->val->data[i0];
  }

  crs_sort(A->row_ptr, L->col_ind, L->val);
  offset = 0;
  start = 0;
  for (i = 1; i <= A->nrows; i++) {
    i0 = L->row_ptr->data[i] - 1;
    while (start + 1 <= i0) {
      if (L->col_ind->data[start] > i) {
        offset++;
      } else {
        if (offset != 0) {
          L->col_ind->data[start - offset] = L->col_ind->data[start];
          L->val->data[start - offset] = L->val->data[start];
        }
      }

      start++;
    }

    start = L->row_ptr->data[i] - 1;
    L->row_ptr->data[i] -= offset;
  }

  if (offset != 0) {
    plcInit_int32_T(&b_A, 1);
    start = L->col_ind->size[0] - offset;
    i0 = b_A->size[0];
    b_A->size[0] = L->col_ind->size[0];
    plcEnsureCapacity((plcArray__common *)b_A, i0, (int32_T)sizeof(int32_T));
    offset = L->col_ind->size[0];
    for (i0 = 0; i0 < offset; i0++) {
      b_A->data[i0] = L->col_ind->data[i0];
    }

    if (start < 1) {
      i0 = L->col_ind->size[0];
      L->col_ind->size[0] = 0;
      plcEnsureCapacity((plcArray__common *)L->col_ind, i0, (int32_T)sizeof
                        (int32_T));
    } else {
      i0 = L->col_ind->size[0];
      L->col_ind->size[0] = start;
      plcEnsureCapacity((plcArray__common *)L->col_ind, i0, (int32_T)sizeof
                        (int32_T));
      for (i = 0; i < start; i++) {
        L->col_ind->data[i] = b_A->data[i];
      }
    }

    plcFree_int32_T(&b_A);
    plcInit_real_T(&c_A, 1);
    i0 = c_A->size[0];
    c_A->size[0] = L->val->size[0];
    plcEnsureCapacity((plcArray__common *)c_A, i0, (int32_T)sizeof(real_T));
    offset = L->val->size[0];
    for (i0 = 0; i0 < offset; i0++) {
      c_A->data[i0] = L->val->data[i0];
    }

    if (start < 1) {
      i0 = L->val->size[0];
      L->val->size[0] = 0;
      plcEnsureCapacity((plcArray__common *)L->val, i0, (int32_T)sizeof(real_T));
    } else {
      i0 = L->val->size[0];
      L->val->size[0] = start;
      plcEnsureCapacity((plcArray__common *)L->val, i0, (int32_T)sizeof(real_T));
      for (i = 0; i < start; i++) {
        L->val->data[i] = c_A->data[i];
      }
    }

    plcFree_real_T(&c_A);
  }
}

void crs_tril1(const struct_T *A, int32_T k, struct_T *L)
{
  int32_T i1;
  int32_T offset;
  int32_T start;
  int32_T i;
  plcArray_int32_T *b_A;
  plcArray_real_T *c_A;
  plcCopyStruct_struct_T(L, A);
  i1 = L->col_ind->size[0];
  L->col_ind->size[0] = A->col_ind->size[0];
  plcEnsureCapacity((plcArray__common *)L->col_ind, i1, (int32_T)sizeof(int32_T));
  offset = A->col_ind->size[0];
  for (i1 = 0; i1 < offset; i1++) {
    L->col_ind->data[i1] = A->col_ind->data[i1];
  }

  i1 = L->val->size[0];
  L->val->size[0] = A->val->size[0];
  plcEnsureCapacity((plcArray__common *)L->val, i1, (int32_T)sizeof(real_T));
  offset = A->val->size[0];
  for (i1 = 0; i1 < offset; i1++) {
    L->val->data[i1] = A->val->data[i1];
  }

  crs_sort(A->row_ptr, L->col_ind, L->val);
  offset = 0;
  start = 0;
  for (i = 1; i <= A->nrows; i++) {
    i1 = L->row_ptr->data[i] - 1;
    while (start + 1 <= i1) {
      if (L->col_ind->data[start] > i + k) {
        offset++;
      } else {
        if (offset != 0) {
          L->col_ind->data[start - offset] = L->col_ind->data[start];
          L->val->data[start - offset] = L->val->data[start];
        }
      }

      start++;
    }

    start = L->row_ptr->data[i] - 1;
    L->row_ptr->data[i] -= offset;
  }

  if (offset != 0) {
    plcInit_int32_T(&b_A, 1);
    start = L->col_ind->size[0] - offset;
    i1 = b_A->size[0];
    b_A->size[0] = L->col_ind->size[0];
    plcEnsureCapacity((plcArray__common *)b_A, i1, (int32_T)sizeof(int32_T));
    offset = L->col_ind->size[0];
    for (i1 = 0; i1 < offset; i1++) {
      b_A->data[i1] = L->col_ind->data[i1];
    }

    if (start < 1) {
      i1 = L->col_ind->size[0];
      L->col_ind->size[0] = 0;
      plcEnsureCapacity((plcArray__common *)L->col_ind, i1, (int32_T)sizeof
                        (int32_T));
    } else {
      i1 = L->col_ind->size[0];
      L->col_ind->size[0] = start;
      plcEnsureCapacity((plcArray__common *)L->col_ind, i1, (int32_T)sizeof
                        (int32_T));
      for (i = 0; i < start; i++) {
        L->col_ind->data[i] = b_A->data[i];
      }
    }

    plcFree_int32_T(&b_A);
    plcInit_real_T(&c_A, 1);
    i1 = c_A->size[0];
    c_A->size[0] = L->val->size[0];
    plcEnsureCapacity((plcArray__common *)c_A, i1, (int32_T)sizeof(real_T));
    offset = L->val->size[0];
    for (i1 = 0; i1 < offset; i1++) {
      c_A->data[i1] = L->val->data[i1];
    }

    if (start < 1) {
      i1 = L->val->size[0];
      L->val->size[0] = 0;
      plcEnsureCapacity((plcArray__common *)L->val, i1, (int32_T)sizeof(real_T));
    } else {
      i1 = L->val->size[0];
      L->val->size[0] = start;
      plcEnsureCapacity((plcArray__common *)L->val, i1, (int32_T)sizeof(real_T));
      for (i = 0; i < start; i++) {
        L->val->data[i] = c_A->data[i];
      }
    }

    plcFree_real_T(&c_A);
  }
}

void crs_tril_initialize(void)
{
}

void crs_tril_terminate(void)
{
}

