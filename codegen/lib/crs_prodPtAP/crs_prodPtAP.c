#include "crs_prodPtAP.h"
#include "m2c.h"

static void crs_create(const m2cArray_int32_T *rows, const m2cArray_int32_T
  *cols, const m2cArray_real_T *vs, m2cArray_int32_T *A_row_ptr,
  m2cArray_int32_T *A_col_ind, m2cArray_real_T *A_val, int32_T *A_nrows, int32_T
  *A_ncols);
static void crs_prodAB(const m2cArray_int32_T *A_row_ptr, const m2cArray_int32_T
  *A_col_ind, const m2cArray_real_T *A_val, int32_T A_nrows, int32_T A_ncols,
  const m2cArray_int32_T *B_row_ptr, const m2cArray_int32_T *B_col_ind, const
  m2cArray_real_T *B_val, int32_T B_nrows, int32_T B_ncols, m2cArray_int32_T
  *C_row_ptr, m2cArray_int32_T *C_col_ind, m2cArray_real_T *C_val, int32_T
  *C_nrows, int32_T *C_ncols);
static void msg_error(void);
static void crs_create(const m2cArray_int32_T *rows, const m2cArray_int32_T
  *cols, const m2cArray_real_T *vs, m2cArray_int32_T *A_row_ptr,
  m2cArray_int32_T *A_col_ind, m2cArray_real_T *A_val, int32_T *A_nrows, int32_T
  *A_ncols)
{
  m2cArray_real_T *r0;
  m2cArray_int32_T *x;
  m2cArray_real_T *buf_val;
  m2cArray_int32_T *buf_indx;
  int32_T i3;
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
    i3 = A_row_ptr->size[0];
    A_row_ptr->size[0] = 0;
    m2cEnsureCapacity((m2cArray__common *)A_row_ptr, i3, (int32_T)sizeof(int32_T));
    i3 = A_col_ind->size[0];
    A_col_ind->size[0] = 0;
    m2cEnsureCapacity((m2cArray__common *)A_col_ind, i3, (int32_T)sizeof(int32_T));
    i3 = A_val->size[0];
    A_val->size[0] = vs->size[0];
    m2cEnsureCapacity((m2cArray__common *)A_val, i3, (int32_T)sizeof(real_T));
    ix = vs->size[0];
    for (i3 = 0; i3 < ix; i3++) {
      A_val->data[i3] = vs->data[i3];
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

    i3 = A_row_ptr->size[0];
    A_row_ptr->size[0] = l + 1;
    m2cEnsureCapacity((m2cArray__common *)A_row_ptr, i3, (int32_T)sizeof(int32_T));
    for (i3 = 0; i3 <= l; i3++) {
      A_row_ptr->data[i3] = 0;
    }

    *A_nrows = l;
    *A_ncols = n;
    i3 = rows->size[0];
    for (i = 0; i + 1 <= i3; i++) {
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
      i3 = A_col_ind->size[0];
      A_col_ind->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_col_ind, i3, (int32_T)sizeof
                        (int32_T));
      ix = cols->size[0];
      for (i3 = 0; i3 < ix; i3++) {
        A_col_ind->data[i3] = cols->data[i3];
      }

      i3 = A_val->size[0];
      A_val->size[0] = vs->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_val, i3, (int32_T)sizeof(real_T));
      ix = vs->size[0];
      for (i3 = 0; i3 < ix; i3++) {
        A_val->data[i3] = vs->data[i3];
      }
    } else {
      i3 = A_col_ind->size[0];
      A_col_ind->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_col_ind, i3, (int32_T)sizeof
                        (int32_T));
      i3 = A_val->size[0];
      A_val->size[0] = cols->size[0];
      m2cEnsureCapacity((m2cArray__common *)A_val, i3, (int32_T)sizeof(real_T));
      for (i = 0; i < rows->size[0]; i++) {
        j = A_row_ptr->data[rows->data[i] - 1];
        A_val->data[A_row_ptr->data[rows->data[i] - 1] - 1] = vs->data[i];
        A_col_ind->data[j - 1] = cols->data[i];
        A_row_ptr->data[rows->data[i] - 1]++;
      }

      i3 = x->size[0];
      x->size[0] = A_row_ptr->size[0];
      m2cEnsureCapacity((m2cArray__common *)x, i3, (int32_T)sizeof(int32_T));
      ix = A_row_ptr->size[0];
      for (i3 = 0; i3 < ix; i3++) {
        x->data[i3] = A_row_ptr->data[i3];
      }

      i3 = A_row_ptr->size[0];
      for (i = 0; i <= i3 - 2; i++) {
        ix = (x->size[0] - i) - 2;
        A_row_ptr->data[ix + 1] = A_row_ptr->data[ix];
      }

      A_row_ptr->data[0] = 1;
    }

    i3 = r0->size[0];
    r0->size[0] = A_val->size[0];
    m2cEnsureCapacity((m2cArray__common *)r0, i3, (int32_T)sizeof(real_T));
    ix = A_val->size[0];
    for (i3 = 0; i3 < ix; i3++) {
      r0->data[i3] = A_val->data[i3];
    }

    i3 = x->size[0];
    x->size[0] = A_col_ind->size[0];
    m2cEnsureCapacity((m2cArray__common *)x, i3, (int32_T)sizeof(int32_T));
    ix = A_col_ind->size[0];
    for (i3 = 0; i3 < ix; i3++) {
      x->data[i3] = A_col_ind->data[i3];
    }

    i3 = A_row_ptr->size[0] - 1;
    for (i = 1; i <= i3; i++) {
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

    i3 = A_col_ind->size[0];
    A_col_ind->size[0] = x->size[0];
    m2cEnsureCapacity((m2cArray__common *)A_col_ind, i3, (int32_T)sizeof(int32_T));
    ix = x->size[0];
    for (i3 = 0; i3 < ix; i3++) {
      A_col_ind->data[i3] = x->data[i3];
    }

    i3 = A_val->size[0];
    A_val->size[0] = r0->size[0];
    m2cEnsureCapacity((m2cArray__common *)A_val, i3, (int32_T)sizeof(real_T));
    ix = r0->size[0];
    for (i3 = 0; i3 < ix; i3++) {
      A_val->data[i3] = r0->data[i3];
    }
  }

  m2cFree_int32_T(&buf_indx);
  m2cFree_real_T(&buf_val);
  m2cFree_int32_T(&x);
  m2cFree_real_T(&r0);
}

static void crs_prodAB(const m2cArray_int32_T *A_row_ptr, const m2cArray_int32_T
  *A_col_ind, const m2cArray_real_T *A_val, int32_T A_nrows, int32_T A_ncols,
  const m2cArray_int32_T *B_row_ptr, const m2cArray_int32_T *B_col_ind, const
  m2cArray_real_T *B_val, int32_T B_nrows, int32_T B_ncols, m2cArray_int32_T
  *C_row_ptr, m2cArray_int32_T *C_col_ind, m2cArray_real_T *C_val, int32_T
  *C_nrows, int32_T *C_ncols)
{
  m2cArray_int32_T *b_index;
  int32_T i1;
  int32_T maxval;
  int32_T i;
  int32_T istart;
  int32_T clength;
  int32_T i2;
  int32_T k;
  m2cArray_real_T *temp;
  if (A_ncols != B_nrows) {
    msg_error();
  }

  m2cInit_int32_T(&b_index, 1);
  i1 = C_row_ptr->size[0];
  C_row_ptr->size[0] = 0;
  m2cEnsureCapacity((m2cArray__common *)C_row_ptr, i1, (int32_T)sizeof(int32_T));
  i1 = C_col_ind->size[0];
  C_col_ind->size[0] = 0;
  m2cEnsureCapacity((m2cArray__common *)C_col_ind, i1, (int32_T)sizeof(int32_T));
  i1 = C_val->size[0];
  C_val->size[0] = 0;
  m2cEnsureCapacity((m2cArray__common *)C_val, i1, (int32_T)sizeof(real_T));
  *C_nrows = A_nrows;
  *C_ncols = B_ncols;
  i1 = C_row_ptr->size[0];
  C_row_ptr->size[0] = A_row_ptr->size[0];
  m2cEnsureCapacity((m2cArray__common *)C_row_ptr, i1, (int32_T)sizeof(int32_T));
  C_row_ptr->data[0] = 1;
  if (A_ncols >= B_ncols) {
    maxval = A_ncols;
  } else {
    maxval = B_ncols;
  }

  i1 = b_index->size[0];
  b_index->size[0] = maxval;
  m2cEnsureCapacity((m2cArray__common *)b_index, i1, (int32_T)sizeof(int32_T));
  for (i1 = 0; i1 < maxval; i1++) {
    b_index->data[i1] = 0;
  }

  for (i = 1; i <= A_nrows; i++) {
    istart = -1;
    clength = 0;
    i1 = A_row_ptr->data[i] - 1;
    for (maxval = A_row_ptr->data[i - 1]; maxval <= i1; maxval++) {
      i2 = B_row_ptr->data[A_col_ind->data[maxval - 1]] - 1;
      for (k = B_row_ptr->data[A_col_ind->data[maxval - 1] - 1] - 1; k + 1 <= i2;
           k++) {
        if (b_index->data[B_col_ind->data[k] - 1] == 0) {
          b_index->data[B_col_ind->data[k] - 1] = istart;
          istart = B_col_ind->data[k];
          clength++;
        }
      }
    }

    C_row_ptr->data[i] = C_row_ptr->data[i - 1] + clength;
    i1 = C_row_ptr->data[i] - 1;
    for (maxval = C_row_ptr->data[i - 1]; maxval <= i1; maxval++) {
      k = istart;
      istart = b_index->data[istart - 1];
      b_index->data[k - 1] = 0;
    }

    b_index->data[i - 1] = 0;
  }

  i1 = C_col_ind->size[0];
  C_col_ind->size[0] = C_row_ptr->data[A_nrows] - 1;
  m2cEnsureCapacity((m2cArray__common *)C_col_ind, i1, (int32_T)sizeof(int32_T));
  for (i = 1; i <= A_nrows; i++) {
    istart = -1;
    clength = 0;
    i1 = A_row_ptr->data[i] - 1;
    for (maxval = A_row_ptr->data[i - 1]; maxval <= i1; maxval++) {
      i2 = B_row_ptr->data[A_col_ind->data[maxval - 1]] - 1;
      for (k = B_row_ptr->data[A_col_ind->data[maxval - 1] - 1] - 1; k + 1 <= i2;
           k++) {
        if (b_index->data[B_col_ind->data[k] - 1] == 0) {
          b_index->data[B_col_ind->data[k] - 1] = istart;
          istart = B_col_ind->data[k];
          clength++;
        }
      }
    }

    C_row_ptr->data[i] = C_row_ptr->data[i - 1] + clength;
    i1 = C_row_ptr->data[i] - 1;
    for (maxval = C_row_ptr->data[i - 1]; maxval <= i1; maxval++) {
      C_col_ind->data[maxval - 1] = istart;
      istart = b_index->data[istart - 1];
      b_index->data[C_col_ind->data[maxval - 1] - 1] = 0;
    }

    b_index->data[i - 1] = 0;
  }

  m2cInit_real_T(&temp, 1);
  i1 = C_val->size[0];
  C_val->size[0] = C_row_ptr->data[A_nrows] - 1;
  m2cEnsureCapacity((m2cArray__common *)C_val, i1, (int32_T)sizeof(real_T));
  i1 = temp->size[0];
  temp->size[0] = b_index->size[0];
  m2cEnsureCapacity((m2cArray__common *)temp, i1, (int32_T)sizeof(real_T));
  maxval = b_index->size[0];
  m2cFree_int32_T(&b_index);
  for (i1 = 0; i1 < maxval; i1++) {
    temp->data[i1] = 0.0;
  }

  for (i = 1; i <= A_nrows; i++) {
    i1 = A_row_ptr->data[i] - 1;
    for (maxval = A_row_ptr->data[i - 1] - 1; maxval + 1 <= i1; maxval++) {
      i2 = B_row_ptr->data[A_col_ind->data[maxval]] - 1;
      for (k = B_row_ptr->data[A_col_ind->data[maxval] - 1] - 1; k + 1 <= i2; k
           ++) {
        temp->data[B_col_ind->data[k] - 1] += A_val->data[maxval] * B_val->
          data[k];
      }
    }

    i1 = C_row_ptr->data[i] - 1;
    for (maxval = C_row_ptr->data[i - 1] - 1; maxval + 1 <= i1; maxval++) {
      C_val->data[maxval] = temp->data[C_col_ind->data[maxval] - 1];
      temp->data[C_col_ind->data[maxval] - 1] = 0.0;
    }
  }

  m2cFree_real_T(&temp);
}

static void msg_error(void)
{
  const char * msgid;
  msgid = "crs_prodMatMat:WrongSizes";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid,
                      "Number of columns of A must be equal to number of rows in B.");
  } else {
    printf("Error %s\nNumber of columns of A must be equal to number of rows in B.",
           msgid);
  }
}

void crs_prodPtAP(const struct_T *A, const struct_T *P, struct_T *B)
{
  m2cArray_int32_T *C_row_ptr;
  m2cArray_int32_T *C_col_ind;
  m2cArray_real_T *C_val;
  m2cArray_int32_T *js;
  int32_T C_ncols;
  int32_T C_nrows;
  uint32_T unnamed_idx_0;
  int32_T i0;
  int32_T nrows;
  int32_T i;
  int32_T j;
  m2cArray_int32_T *Pt_row_ptr;
  m2cArray_int32_T *Pt_col_ind;
  m2cArray_real_T *Pt_val;
  int32_T Pt_ncols;
  int32_T Pt_nrows;
  m2cInit_int32_T(&C_row_ptr, 1);
  m2cInit_int32_T(&C_col_ind, 1);
  m2cInit_real_T(&C_val, 1);
  m2cInit_int32_T(&js, 1);
  crs_prodAB(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, P->row_ptr,
             P->col_ind, P->val, P->nrows, P->ncols, C_row_ptr, C_col_ind, C_val,
             &C_nrows, &C_ncols);
  unnamed_idx_0 = (uint32_T)P->col_ind->size[0];
  i0 = js->size[0];
  js->size[0] = (int32_T)unnamed_idx_0;
  m2cEnsureCapacity((m2cArray__common *)js, i0, (int32_T)sizeof(int32_T));
  nrows = P->row_ptr->size[0] - 1;
  for (i = 1; i <= nrows; i++) {
    i0 = P->row_ptr->data[i] - 1;
    for (j = P->row_ptr->data[i - 1]; j <= i0; j++) {
      js->data[j - 1] = i;
    }
  }

  m2cInit_int32_T(&Pt_row_ptr, 1);
  m2cInit_int32_T(&Pt_col_ind, 1);
  m2cInit_real_T(&Pt_val, 1);
  crs_create(P->col_ind, js, P->val, Pt_row_ptr, Pt_col_ind, Pt_val, &Pt_nrows,
             &Pt_ncols);
  crs_prodAB(Pt_row_ptr, Pt_col_ind, Pt_val, Pt_nrows, Pt_ncols, C_row_ptr,
             C_col_ind, C_val, C_nrows, C_ncols, B->row_ptr, B->col_ind, B->val,
             &B->nrows, &B->ncols);
  m2cFree_int32_T(&js);
  m2cFree_real_T(&Pt_val);
  m2cFree_int32_T(&Pt_col_ind);
  m2cFree_int32_T(&Pt_row_ptr);
  m2cFree_real_T(&C_val);
  m2cFree_int32_T(&C_col_ind);
  m2cFree_int32_T(&C_row_ptr);
}

void crs_prodPtAP_initialize(void)
{
}

void crs_prodPtAP_terminate(void)
{
}

