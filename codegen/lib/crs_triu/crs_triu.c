#include "crs_triu.h"
#include "m2c.h"

static void crs_sort(const m2cArray_int32_T *row_ptr, m2cArray_int32_T *col_ind,
                     m2cArray_real_T *val);
static void m2cCopyStruct_struct_T(struct_T *dst, const struct_T *src);
static define_m2cCopy(m2cCopy_int32_T, int32_T);
static define_m2cCopy(m2cCopy_real_T, real_T);
static void crs_sort(const m2cArray_int32_T *row_ptr, m2cArray_int32_T *col_ind,
                     m2cArray_real_T *val)
{
  int32_T i2;
  int32_T i;
  m2cArray_real_T *buf_val;
  m2cArray_int32_T *buf_indx;
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
  m2cInit_real_T(&buf_val, 1);
  m2cInit_int32_T(&buf_indx, 1);
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
      m2cEnsureCapacity((m2cArray__common *)buf_indx, n, (int32_T)sizeof(int32_T));
      n = buf_val->size[0];
      buf_val->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      m2cEnsureCapacity((m2cArray__common *)buf_val, n, (int32_T)sizeof(real_T));
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

  m2cFree_int32_T(&buf_indx);
  m2cFree_real_T(&buf_val);
}

static void m2cCopyStruct_struct_T(struct_T *dst, const struct_T *src)
{
  m2cCopy_int32_T(&dst->row_ptr, &src->row_ptr);
  m2cCopy_int32_T(&dst->col_ind, &src->col_ind);
  m2cCopy_real_T(&dst->val, &src->val);
  dst->nrows = src->nrows;
  dst->ncols = src->ncols;
}

void crs_triu(const struct_T *A, struct_T *U)
{
  int32_T i0;
  int32_T offset;
  int32_T start;
  int32_T i;
  m2cArray_int32_T *b_A;
  m2cArray_real_T *c_A;
  m2cCopyStruct_struct_T(U, A);
  i0 = U->col_ind->size[0];
  U->col_ind->size[0] = A->col_ind->size[0];
  m2cEnsureCapacity((m2cArray__common *)U->col_ind, i0, (int32_T)sizeof(int32_T));
  offset = A->col_ind->size[0];
  for (i0 = 0; i0 < offset; i0++) {
    U->col_ind->data[i0] = A->col_ind->data[i0];
  }

  i0 = U->val->size[0];
  U->val->size[0] = A->val->size[0];
  m2cEnsureCapacity((m2cArray__common *)U->val, i0, (int32_T)sizeof(real_T));
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
    m2cInit_int32_T(&b_A, 1);
    start = U->col_ind->size[0] - offset;
    i0 = b_A->size[0];
    b_A->size[0] = U->col_ind->size[0];
    m2cEnsureCapacity((m2cArray__common *)b_A, i0, (int32_T)sizeof(int32_T));
    offset = U->col_ind->size[0];
    for (i0 = 0; i0 < offset; i0++) {
      b_A->data[i0] = U->col_ind->data[i0];
    }

    if (start < 1) {
      i0 = U->col_ind->size[0];
      U->col_ind->size[0] = 0;
      m2cEnsureCapacity((m2cArray__common *)U->col_ind, i0, (int32_T)sizeof
                        (int32_T));
    } else {
      i0 = U->col_ind->size[0];
      U->col_ind->size[0] = start;
      m2cEnsureCapacity((m2cArray__common *)U->col_ind, i0, (int32_T)sizeof
                        (int32_T));
      for (i = 0; i < start; i++) {
        U->col_ind->data[i] = b_A->data[i];
      }
    }

    m2cFree_int32_T(&b_A);
    m2cInit_real_T(&c_A, 1);
    i0 = c_A->size[0];
    c_A->size[0] = U->val->size[0];
    m2cEnsureCapacity((m2cArray__common *)c_A, i0, (int32_T)sizeof(real_T));
    offset = U->val->size[0];
    for (i0 = 0; i0 < offset; i0++) {
      c_A->data[i0] = U->val->data[i0];
    }

    if (start < 1) {
      i0 = U->val->size[0];
      U->val->size[0] = 0;
      m2cEnsureCapacity((m2cArray__common *)U->val, i0, (int32_T)sizeof(real_T));
    } else {
      i0 = U->val->size[0];
      U->val->size[0] = start;
      m2cEnsureCapacity((m2cArray__common *)U->val, i0, (int32_T)sizeof(real_T));
      for (i = 0; i < start; i++) {
        U->val->data[i] = c_A->data[i];
      }
    }

    m2cFree_real_T(&c_A);
  }
}

void crs_triu1(const struct_T *A, int32_T k, struct_T *U)
{
  int32_T i1;
  int32_T offset;
  int32_T start;
  int32_T i;
  m2cArray_int32_T *b_A;
  m2cArray_real_T *c_A;
  m2cCopyStruct_struct_T(U, A);
  i1 = U->col_ind->size[0];
  U->col_ind->size[0] = A->col_ind->size[0];
  m2cEnsureCapacity((m2cArray__common *)U->col_ind, i1, (int32_T)sizeof(int32_T));
  offset = A->col_ind->size[0];
  for (i1 = 0; i1 < offset; i1++) {
    U->col_ind->data[i1] = A->col_ind->data[i1];
  }

  i1 = U->val->size[0];
  U->val->size[0] = A->val->size[0];
  m2cEnsureCapacity((m2cArray__common *)U->val, i1, (int32_T)sizeof(real_T));
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
    m2cInit_int32_T(&b_A, 1);
    start = U->col_ind->size[0] - offset;
    i1 = b_A->size[0];
    b_A->size[0] = U->col_ind->size[0];
    m2cEnsureCapacity((m2cArray__common *)b_A, i1, (int32_T)sizeof(int32_T));
    offset = U->col_ind->size[0];
    for (i1 = 0; i1 < offset; i1++) {
      b_A->data[i1] = U->col_ind->data[i1];
    }

    if (start < 1) {
      i1 = U->col_ind->size[0];
      U->col_ind->size[0] = 0;
      m2cEnsureCapacity((m2cArray__common *)U->col_ind, i1, (int32_T)sizeof
                        (int32_T));
    } else {
      i1 = U->col_ind->size[0];
      U->col_ind->size[0] = start;
      m2cEnsureCapacity((m2cArray__common *)U->col_ind, i1, (int32_T)sizeof
                        (int32_T));
      for (i = 0; i < start; i++) {
        U->col_ind->data[i] = b_A->data[i];
      }
    }

    m2cFree_int32_T(&b_A);
    m2cInit_real_T(&c_A, 1);
    i1 = c_A->size[0];
    c_A->size[0] = U->val->size[0];
    m2cEnsureCapacity((m2cArray__common *)c_A, i1, (int32_T)sizeof(real_T));
    offset = U->val->size[0];
    for (i1 = 0; i1 < offset; i1++) {
      c_A->data[i1] = U->val->data[i1];
    }

    if (start < 1) {
      i1 = U->val->size[0];
      U->val->size[0] = 0;
      m2cEnsureCapacity((m2cArray__common *)U->val, i1, (int32_T)sizeof(real_T));
    } else {
      i1 = U->val->size[0];
      U->val->size[0] = start;
      m2cEnsureCapacity((m2cArray__common *)U->val, i1, (int32_T)sizeof(real_T));
      for (i = 0; i < start; i++) {
        U->val->data[i] = c_A->data[i];
      }
    }

    m2cFree_real_T(&c_A);
  }
}

void crs_triu_initialize(void)
{
}

void crs_triu_terminate(void)
{
}

