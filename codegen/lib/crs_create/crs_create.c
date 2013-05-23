#include "crs_create.h"
#include "m2c.h"

static void crs_sort(const m2cArray_int32_T *row_ptr, m2cArray_int32_T *col_ind,
                     m2cArray_real_T *val);
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

void crs_create(const m2cArray_int32_T *rows, const m2cArray_int32_T *cols,
                const m2cArray_real_T *vs, struct_T *A)
{
  m2cArray_int32_T *x;
  int32_T i0;
  int32_T j;
  int32_T mtmp;
  int32_T b_mtmp;
  boolean_T ascend;
  boolean_T exitg1;
  m2cInit_int32_T(&x, 1);
  if ((rows->size[0] == 1) && (cols->size[0] == 1)) {
    i0 = A->row_ptr->size[0];
    A->row_ptr->size[0] = 0;
    m2cEnsureCapacity((m2cArray__common *)A->row_ptr, i0, (int32_T)sizeof
                      (int32_T));
    i0 = A->col_ind->size[0];
    A->col_ind->size[0] = 0;
    m2cEnsureCapacity((m2cArray__common *)A->col_ind, i0, (int32_T)sizeof
                      (int32_T));
    i0 = A->val->size[0];
    A->val->size[0] = vs->size[0];
    m2cEnsureCapacity((m2cArray__common *)A->val, i0, (int32_T)sizeof(real_T));
    j = vs->size[0];
    for (i0 = 0; i0 < j; i0++) {
      A->val->data[i0] = vs->data[i0];
    }

    A->nrows = rows->data[0];
    A->ncols = cols->data[0];
  } else {
    mtmp = rows->data[0];
    if ((rows->size[0] > 1) && (1 < rows->size[0])) {
      for (j = 1; j + 1 <= rows->size[0]; j++) {
        if (rows->data[j] > mtmp) {
          mtmp = rows->data[j];
        }
      }
    }

    b_mtmp = cols->data[0];
    if ((cols->size[0] > 1) && (1 < cols->size[0])) {
      for (j = 1; j + 1 <= cols->size[0]; j++) {
        if (cols->data[j] > b_mtmp) {
          b_mtmp = cols->data[j];
        }
      }
    }

    i0 = A->row_ptr->size[0];
    A->row_ptr->size[0] = mtmp + 1;
    m2cEnsureCapacity((m2cArray__common *)A->row_ptr, i0, (int32_T)sizeof
                      (int32_T));
    for (i0 = 0; i0 <= mtmp; i0++) {
      A->row_ptr->data[i0] = 0;
    }

    A->nrows = mtmp;
    A->ncols = b_mtmp;
    i0 = rows->size[0];
    for (b_mtmp = 0; b_mtmp + 1 <= i0; b_mtmp++) {
      A->row_ptr->data[rows->data[b_mtmp]]++;
    }

    A->row_ptr->data[0] = 1;
    for (b_mtmp = 1; b_mtmp <= mtmp; b_mtmp++) {
      A->row_ptr->data[b_mtmp] += A->row_ptr->data[b_mtmp - 1];
    }

    ascend = TRUE;
    b_mtmp = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (b_mtmp <= rows->size[0] - 2)) {
      if (rows->data[b_mtmp + 1] < rows->data[b_mtmp]) {
        ascend = FALSE;
        exitg1 = TRUE;
      } else {
        b_mtmp++;
      }
    }

    if (ascend) {
      i0 = A->col_ind->size[0];
      A->col_ind->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A->col_ind, i0, (int32_T)sizeof
                        (int32_T));
      j = cols->size[0];
      for (i0 = 0; i0 < j; i0++) {
        A->col_ind->data[i0] = cols->data[i0];
      }

      i0 = A->val->size[0];
      A->val->size[0] = vs->size[0];
      m2cEnsureCapacity((m2cArray__common *)A->val, i0, (int32_T)sizeof(real_T));
      j = vs->size[0];
      for (i0 = 0; i0 < j; i0++) {
        A->val->data[i0] = vs->data[i0];
      }
    } else {
      i0 = A->col_ind->size[0];
      A->col_ind->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A->col_ind, i0, (int32_T)sizeof
                        (int32_T));
      i0 = A->val->size[0];
      A->val->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A->val, i0, (int32_T)sizeof(real_T));
      for (b_mtmp = 0; b_mtmp < rows->size[0]; b_mtmp++) {
        j = A->row_ptr->data[rows->data[b_mtmp] - 1];
        A->val->data[A->row_ptr->data[rows->data[b_mtmp] - 1] - 1] = vs->
          data[b_mtmp];
        A->col_ind->data[j - 1] = cols->data[b_mtmp];
        A->row_ptr->data[rows->data[b_mtmp] - 1]++;
      }

      i0 = x->size[0];
      x->size[0] = A->row_ptr->size[0];
      m2cEnsureCapacity((m2cArray__common *)x, i0, (int32_T)sizeof(int32_T));
      j = A->row_ptr->size[0];
      for (i0 = 0; i0 < j; i0++) {
        x->data[i0] = A->row_ptr->data[i0];
      }

      i0 = A->row_ptr->size[0];
      for (b_mtmp = 0; b_mtmp <= i0 - 2; b_mtmp++) {
        j = (x->size[0] - b_mtmp) - 2;
        A->row_ptr->data[j + 1] = A->row_ptr->data[j];
      }

      A->row_ptr->data[0] = 1;
    }

    crs_sort(A->row_ptr, A->col_ind, A->val);
  }

  m2cFree_int32_T(&x);
}

b_struct_T crs_create0(int32_T ni, int32_T nj)
{
  b_struct_T A;
  A.nrows = ni;
  A.ncols = nj;
  return A;
}

void crs_create1(const m2cArray_int32_T *is, const m2cArray_int32_T *js, const
                 m2cArray_real_T *vs, int32_T ni, int32_T nj, struct_T *A)
{
  int32_T i1;
  int32_T i;
  boolean_T ascend;
  boolean_T exitg1;
  int32_T j;
  m2cArray_int32_T *x;
  i1 = A->row_ptr->size[0];
  A->row_ptr->size[0] = ni + 1;
  m2cEnsureCapacity((m2cArray__common *)A->row_ptr, i1, (int32_T)sizeof(int32_T));
  for (i1 = 0; i1 <= ni; i1++) {
    A->row_ptr->data[i1] = 0;
  }

  A->nrows = ni;
  A->ncols = nj;
  i1 = is->size[0];
  for (i = 0; i + 1 <= i1; i++) {
    A->row_ptr->data[is->data[i]]++;
  }

  A->row_ptr->data[0] = 1;
  for (i = 1; i <= ni; i++) {
    A->row_ptr->data[i] += A->row_ptr->data[i - 1];
  }

  ascend = TRUE;
  i = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i <= is->size[0] - 2)) {
    if (is->data[i + 1] < is->data[i]) {
      ascend = FALSE;
      exitg1 = TRUE;
    } else {
      i++;
    }
  }

  if (ascend) {
    i1 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    m2cEnsureCapacity((m2cArray__common *)A->col_ind, i1, (int32_T)sizeof
                      (int32_T));
    j = js->size[0];
    for (i1 = 0; i1 < j; i1++) {
      A->col_ind->data[i1] = js->data[i1];
    }

    i1 = A->val->size[0];
    A->val->size[0] = vs->size[0];
    m2cEnsureCapacity((m2cArray__common *)A->val, i1, (int32_T)sizeof(real_T));
    j = vs->size[0];
    for (i1 = 0; i1 < j; i1++) {
      A->val->data[i1] = vs->data[i1];
    }
  } else {
    i1 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    m2cEnsureCapacity((m2cArray__common *)A->col_ind, i1, (int32_T)sizeof
                      (int32_T));
    i1 = A->val->size[0];
    A->val->size[0] = js->size[0];
    m2cEnsureCapacity((m2cArray__common *)A->val, i1, (int32_T)sizeof(real_T));
    for (i = 0; i < is->size[0]; i++) {
      j = A->row_ptr->data[is->data[i] - 1];
      A->val->data[A->row_ptr->data[is->data[i] - 1] - 1] = vs->data[i];
      A->col_ind->data[j - 1] = js->data[i];
      A->row_ptr->data[is->data[i] - 1]++;
    }

    m2cInit_int32_T(&x, 1);
    i1 = x->size[0];
    x->size[0] = A->row_ptr->size[0];
    m2cEnsureCapacity((m2cArray__common *)x, i1, (int32_T)sizeof(int32_T));
    j = A->row_ptr->size[0];
    for (i1 = 0; i1 < j; i1++) {
      x->data[i1] = A->row_ptr->data[i1];
    }

    i1 = A->row_ptr->size[0];
    for (i = 0; i <= i1 - 2; i++) {
      j = (x->size[0] - i) - 2;
      A->row_ptr->data[j + 1] = A->row_ptr->data[j];
    }

    m2cFree_int32_T(&x);
    A->row_ptr->data[0] = 1;
  }

  crs_sort(A->row_ptr, A->col_ind, A->val);
}

void crs_create_initialize(void)
{
}

void crs_create_terminate(void)
{
}

