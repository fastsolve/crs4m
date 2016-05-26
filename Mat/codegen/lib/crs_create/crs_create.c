#include "crs_create.h"
#include "m2c.h"

static void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
                     emxArray_real_T *val);
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

void crs_create(const emxArray_int32_T *rows, const emxArray_int32_T *cols,
                const emxArray_real_T *vs, struct0_T *A)
{
  int mtmp;
  int i0;
  int ix;
  int b_mtmp;
  int i;
  boolean_T ascend;
  boolean_T exitg1;
  int j;
  int b_i;
  if ((rows->size[0] == 1) && (cols->size[0] == 1)) {
    i0 = A->row_ptr->size[0];
    A->row_ptr->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)A->row_ptr, i0, (int)sizeof(int));
    i0 = A->col_ind->size[0];
    A->col_ind->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)A->col_ind, i0, (int)sizeof(int));
    i0 = A->val->size[0];
    A->val->size[0] = vs->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i0, (int)sizeof(double));
    ix = vs->size[0];
    for (i0 = 0; i0 < ix; i0++) {
      A->val->data[i0] = vs->data[i0];
    }

    A->nrows = rows->data[0];
    A->ncols = cols->data[0];
  } else {
    mtmp = rows->data[0];
    if (rows->size[0] > 1) {
      for (ix = 1; ix + 1 <= rows->size[0]; ix++) {
        if (rows->data[ix] > mtmp) {
          mtmp = rows->data[ix];
        }
      }
    }

    b_mtmp = cols->data[0];
    if (cols->size[0] > 1) {
      for (ix = 1; ix + 1 <= cols->size[0]; ix++) {
        if (cols->data[ix] > b_mtmp) {
          b_mtmp = cols->data[ix];
        }
      }
    }

    i0 = A->row_ptr->size[0];
    A->row_ptr->size[0] = mtmp + 1;
    emxEnsureCapacity((emxArray__common *)A->row_ptr, i0, (int)sizeof(int));
    for (i0 = 0; i0 <= mtmp; i0++) {
      A->row_ptr->data[i0] = 0;
    }

    A->nrows = mtmp;
    A->ncols = b_mtmp;
    i0 = rows->size[0];
    for (i = 0; i + 1 <= i0; i++) {
      A->row_ptr->data[rows->data[i]]++;
    }

    A->row_ptr->data[0] = 1;
    for (i = 1; i <= mtmp; i++) {
      A->row_ptr->data[i] += A->row_ptr->data[i - 1];
    }

    ascend = true;
    i = 0;
    exitg1 = false;
    while ((!exitg1) && (i <= rows->size[0] - 2)) {
      if (rows->data[1 + i] < rows->data[i]) {
        ascend = false;
        exitg1 = true;
      } else {
        i++;
      }
    }

    if (ascend) {
      i0 = A->col_ind->size[0];
      A->col_ind->size[0] = cols->size[0];
      emxEnsureCapacity((emxArray__common *)A->col_ind, i0, (int)sizeof(int));
      ix = cols->size[0];
      for (i0 = 0; i0 < ix; i0++) {
        A->col_ind->data[i0] = cols->data[i0];
      }

      i0 = A->val->size[0];
      A->val->size[0] = vs->size[0];
      emxEnsureCapacity((emxArray__common *)A->val, i0, (int)sizeof(double));
      ix = vs->size[0];
      for (i0 = 0; i0 < ix; i0++) {
        A->val->data[i0] = vs->data[i0];
      }
    } else {
      i0 = A->col_ind->size[0];
      A->col_ind->size[0] = cols->size[0];
      emxEnsureCapacity((emxArray__common *)A->col_ind, i0, (int)sizeof(int));
      i0 = A->val->size[0];
      A->val->size[0] = cols->size[0];
      emxEnsureCapacity((emxArray__common *)A->val, i0, (int)sizeof(double));
      for (i = 0; i < rows->size[0]; i++) {
        j = A->row_ptr->data[rows->data[i] - 1];
        A->val->data[A->row_ptr->data[rows->data[i] - 1] - 1] = vs->data[i];
        A->col_ind->data[j - 1] = cols->data[i];
        A->row_ptr->data[rows->data[i] - 1]++;
      }

      ix = A->row_ptr->size[0] - 2;
      i0 = A->row_ptr->size[0];
      for (i = 0; i <= i0 - 2; i++) {
        b_i = ix - i;
        A->row_ptr->data[b_i + 1] = A->row_ptr->data[b_i];
      }

      A->row_ptr->data[0] = 1;
    }

    crs_sort(A->row_ptr, A->col_ind, A->val);
  }
}

void crs_create0(int ni, int nj, struct1_T *A)
{
  A->nrows = ni;
  A->ncols = nj;
}

void crs_create1(const emxArray_int32_T *is, const emxArray_int32_T *js, const
                 emxArray_real_T *vs, int ni, int nj, struct0_T *A)
{
  int i1;
  int i;
  boolean_T ascend;
  boolean_T exitg1;
  int varargin_1;
  int j;
  int b_i;
  i1 = A->row_ptr->size[0];
  A->row_ptr->size[0] = ni + 1;
  emxEnsureCapacity((emxArray__common *)A->row_ptr, i1, (int)sizeof(int));
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

  ascend = true;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i <= is->size[0] - 2)) {
    if (is->data[1 + i] < is->data[i]) {
      ascend = false;
      exitg1 = true;
    } else {
      i++;
    }
  }

  if (ascend) {
    i1 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->col_ind, i1, (int)sizeof(int));
    varargin_1 = js->size[0];
    for (i1 = 0; i1 < varargin_1; i1++) {
      A->col_ind->data[i1] = js->data[i1];
    }

    i1 = A->val->size[0];
    A->val->size[0] = vs->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i1, (int)sizeof(double));
    varargin_1 = vs->size[0];
    for (i1 = 0; i1 < varargin_1; i1++) {
      A->val->data[i1] = vs->data[i1];
    }
  } else {
    i1 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->col_ind, i1, (int)sizeof(int));
    i1 = A->val->size[0];
    A->val->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i1, (int)sizeof(double));
    for (i = 0; i < is->size[0]; i++) {
      j = A->row_ptr->data[is->data[i] - 1];
      A->val->data[A->row_ptr->data[is->data[i] - 1] - 1] = vs->data[i];
      A->col_ind->data[j - 1] = js->data[i];
      A->row_ptr->data[is->data[i] - 1]++;
    }

    varargin_1 = A->row_ptr->size[0];
    i1 = A->row_ptr->size[0];
    for (i = 1; i - 1 <= i1 - 2; i++) {
      b_i = varargin_1 - i;
      A->row_ptr->data[b_i] = A->row_ptr->data[b_i - 1];
    }

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

void emxDestroy_struct0_T(struct0_T emxArray)
{
  emxFreeStruct_struct0_T(&emxArray);
}

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}

void emxInit_struct0_T(struct0_T *pStruct)
{
  emxInitStruct_struct0_T(pStruct);
}
