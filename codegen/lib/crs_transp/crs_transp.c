#include "crs_transp.h"
#include "m2c.h"

static void crs_create(const m2cArray_int32_T *rows, const m2cArray_int32_T
  *cols, const m2cArray_real_T *vs, m2cArray_int32_T *A_row_ptr,
  m2cArray_int32_T *A_col_ind, m2cArray_real_T *A_val, int32_T *A_nrows, int32_T
  *A_ncols);
static void crs_create(const m2cArray_int32_T *rows, const m2cArray_int32_T
  *cols, const m2cArray_real_T *vs, m2cArray_int32_T *A_row_ptr,
  m2cArray_int32_T *A_col_ind, m2cArray_real_T *A_val, int32_T *A_nrows, int32_T
  *A_ncols)
{
  m2cArray_real_T *r0;
  m2cArray_int32_T *x;
  m2cArray_real_T *buf_val;
  m2cArray_int32_T *buf_indx;
  int32_T i1;
  int32_T ix;
  int32_T l;
  int32_T n;
  int32_T i;
  boolean_T ascend;
  boolean_T exitg4;
  int32_T j;
  boolean_T exitg3;
  uint32_T ind;
  int32_T exitg1;
  boolean_T guard1 = FALSE;
  int32_T b_r0;
  real_T t0;
  int32_T exitg2;
  boolean_T guard2 = FALSE;
  m2cInit_real_T(&r0, 1);
  m2cInit_int32_T(&x, 1);
  m2cInit_real_T(&buf_val, 1);
  m2cInit_int32_T(&buf_indx, 1);
  if ((rows->size[0] == 1) && (cols->size[0] == 1)) {
    i1 = A_row_ptr->size[0];
    A_row_ptr->size[0] = 0;
    m2cEnsureCapacity((m2cArray__common *)A_row_ptr, i1, (int32_T)sizeof(int32_T));
    i1 = A_col_ind->size[0];
    A_col_ind->size[0] = 0;
    m2cEnsureCapacity((m2cArray__common *)A_col_ind, i1, (int32_T)sizeof(int32_T));
    i1 = A_val->size[0];
    A_val->size[0] = vs->size[0];
    m2cEnsureCapacity((m2cArray__common *)A_val, i1, (int32_T)sizeof(real_T));
    ix = vs->size[0];
    for (i1 = 0; i1 < ix; i1++) {
      A_val->data[i1] = vs->data[i1];
    }

    *A_nrows = rows->data[0];
    *A_ncols = cols->data[0];
  } else {
    l = rows->data[0];
    if ((rows->size[0] > 1) && (1 < rows->size[0])) {
      for (ix = 1; ix + 1 <= rows->size[0]; ix++) {
        if (rows->data[ix] > l) {
          l = rows->data[ix];
        }
      }
    }

    n = cols->data[0];
    if ((cols->size[0] > 1) && (1 < cols->size[0])) {
      for (ix = 1; ix + 1 <= cols->size[0]; ix++) {
        if (cols->data[ix] > n) {
          n = cols->data[ix];
        }
      }
    }

    i1 = A_row_ptr->size[0];
    A_row_ptr->size[0] = l + 1;
    m2cEnsureCapacity((m2cArray__common *)A_row_ptr, i1, (int32_T)sizeof(int32_T));
    for (i1 = 0; i1 <= l; i1++) {
      A_row_ptr->data[i1] = 0;
    }

    *A_nrows = l;
    *A_ncols = n;
    i1 = rows->size[0];
    for (i = 0; i + 1 <= i1; i++) {
      A_row_ptr->data[rows->data[i]]++;
    }

    A_row_ptr->data[0] = 1;
    for (i = 1; i <= l; i++) {
      A_row_ptr->data[i] += A_row_ptr->data[i - 1];
    }

    ascend = TRUE;
    i = 0;
    exitg4 = FALSE;
    while ((exitg4 == FALSE) && (i <= rows->size[0] - 2)) {
      if (rows->data[i + 1] < rows->data[i]) {
        ascend = FALSE;
        exitg4 = TRUE;
      } else {
        i++;
      }
    }

    if (ascend) {
      i1 = A_col_ind->size[0];
      A_col_ind->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_col_ind, i1, (int32_T)sizeof
                        (int32_T));
      ix = cols->size[0];
      for (i1 = 0; i1 < ix; i1++) {
        A_col_ind->data[i1] = cols->data[i1];
      }

      i1 = A_val->size[0];
      A_val->size[0] = vs->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_val, i1, (int32_T)sizeof(real_T));
      ix = vs->size[0];
      for (i1 = 0; i1 < ix; i1++) {
        A_val->data[i1] = vs->data[i1];
      }
    } else {
      i1 = A_col_ind->size[0];
      A_col_ind->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_col_ind, i1, (int32_T)sizeof
                        (int32_T));
      i1 = A_val->size[0];
      A_val->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_val, i1, (int32_T)sizeof(real_T));
      for (i = 0; i < rows->size[0]; i++) {
        j = A_row_ptr->data[rows->data[i] - 1];
        A_val->data[A_row_ptr->data[rows->data[i] - 1] - 1] = vs->data[i];
        A_col_ind->data[j - 1] = cols->data[i];
        A_row_ptr->data[rows->data[i] - 1]++;
      }

      i1 = x->size[0];
      x->size[0] = A_row_ptr->size[0];
      m2cEnsureCapacity((m2cArray__common *)x, i1, (int32_T)sizeof(int32_T));
      ix = A_row_ptr->size[0];
      for (i1 = 0; i1 < ix; i1++) {
        x->data[i1] = A_row_ptr->data[i1];
      }

      i1 = A_row_ptr->size[0];
      for (i = 0; i <= i1 - 2; i++) {
        ix = (x->size[0] - i) - 2;
        A_row_ptr->data[ix + 1] = A_row_ptr->data[ix];
      }

      A_row_ptr->data[0] = 1;
    }

    i1 = r0->size[0];
    r0->size[0] = A_val->size[0];
    m2cEnsureCapacity((m2cArray__common *)r0, i1, (int32_T)sizeof(real_T));
    ix = A_val->size[0];
    for (i1 = 0; i1 < ix; i1++) {
      r0->data[i1] = A_val->data[i1];
    }

    i1 = x->size[0];
    x->size[0] = A_col_ind->size[0];
    m2cEnsureCapacity((m2cArray__common *)x, i1, (int32_T)sizeof(int32_T));
    ix = A_col_ind->size[0];
    for (i1 = 0; i1 < ix; i1++) {
      x->data[i1] = A_col_ind->data[i1];
    }

    i1 = A_row_ptr->size[0] - 1;
    for (i = 1; i <= i1; i++) {
      ascend = TRUE;
      j = A_row_ptr->data[i - 1];
      ix = A_row_ptr->data[i] - 1;
      exitg3 = FALSE;
      while ((exitg3 == FALSE) && (j + 1 <= ix)) {
        if (x->data[j] < x->data[j - 1]) {
          ascend = FALSE;
          exitg3 = TRUE;
        } else {
          j++;
        }
      }

      if (!ascend) {
        ix = buf_indx->size[0];
        buf_indx->size[0] = A_row_ptr->data[i] - A_row_ptr->data[i - 1];
        m2cEnsureCapacity((m2cArray__common *)buf_indx, ix, (int32_T)sizeof
                          (int32_T));
        ix = buf_val->size[0];
        buf_val->size[0] = A_row_ptr->data[i] - A_row_ptr->data[i - 1];
        m2cEnsureCapacity((m2cArray__common *)buf_val, ix, (int32_T)sizeof
                          (real_T));
        ind = 1U;
        ix = A_row_ptr->data[i] - 1;
        for (j = A_row_ptr->data[i - 1]; j <= ix; j++) {
          buf_indx->data[(int32_T)ind - 1] = x->data[j - 1];
          buf_val->data[(int32_T)ind - 1] = r0->data[j - 1];
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
              b_r0 = buf_indx->data[n - 1];
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
              b_r0 = buf_indx->data[l];
              t0 = buf_val->data[l];
              guard1 = TRUE;
            }

            if (guard1 == TRUE) {
              j = l;
              do {
                exitg2 = 0;
                ix = j;
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
                  if ((!ascend) && (buf_indx->data[j] < buf_indx->data[j + 1]))
                  {
                    j++;
                  }

                  if (b_r0 >= buf_indx->data[j]) {
                    exitg2 = 1;
                  } else {
                    buf_indx->data[ix] = buf_indx->data[j];
                    buf_val->data[ix] = buf_val->data[j];
                  }
                }
              } while (exitg2 == 0);

              buf_indx->data[ix] = b_r0;
              buf_val->data[ix] = t0;
            }
          } while (exitg1 == 0);

          buf_indx->data[0] = b_r0;
          buf_val->data[0] = t0;
        }

        ind = 1U;
        ix = A_row_ptr->data[i] - 1;
        for (j = A_row_ptr->data[i - 1]; j <= ix; j++) {
          x->data[j - 1] = buf_indx->data[(int32_T)ind - 1];
          r0->data[j - 1] = buf_val->data[(int32_T)ind - 1];
          ind++;
        }
      }
    }

    i1 = A_col_ind->size[0];
    A_col_ind->size[0] = x->size[0];
    m2cEnsureCapacity((m2cArray__common *)A_col_ind, i1, (int32_T)sizeof(int32_T));
    ix = x->size[0];
    for (i1 = 0; i1 < ix; i1++) {
      A_col_ind->data[i1] = x->data[i1];
    }

    i1 = A_val->size[0];
    A_val->size[0] = r0->size[0];
    m2cEnsureCapacity((m2cArray__common *)A_val, i1, (int32_T)sizeof(real_T));
    ix = r0->size[0];
    for (i1 = 0; i1 < ix; i1++) {
      A_val->data[i1] = r0->data[i1];
    }
  }

  m2cFree_int32_T(&buf_indx);
  m2cFree_real_T(&buf_val);
  m2cFree_int32_T(&x);
  m2cFree_real_T(&r0);
}

void crs_transp(const struct_T *A, struct_T *At)
{
  m2cArray_int32_T *js;
  uint32_T unnamed_idx_0;
  int32_T i0;
  int32_T nrows;
  int32_T i;
  int32_T j;
  m2cInit_int32_T(&js, 1);
  unnamed_idx_0 = (uint32_T)A->col_ind->size[0];
  i0 = js->size[0];
  js->size[0] = (int32_T)unnamed_idx_0;
  m2cEnsureCapacity((m2cArray__common *)js, i0, (int32_T)sizeof(int32_T));
  nrows = A->row_ptr->size[0] - 1;
  for (i = 1; i <= nrows; i++) {
    i0 = A->row_ptr->data[i] - 1;
    for (j = A->row_ptr->data[i - 1]; j <= i0; j++) {
      js->data[j - 1] = i;
    }
  }

  crs_create(A->col_ind, js, A->val, At->row_ptr, At->col_ind, At->val,
             &At->nrows, &At->ncols);
  m2cFree_int32_T(&js);
}

void crs_transp_initialize(void)
{
}

void crs_transp_terminate(void)
{
}

