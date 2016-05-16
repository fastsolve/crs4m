#include "crs_sort.h"
#include "omp.h"
#include "m2c.h"

void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
              emxArray_real_T *val)
{
  int i0;
  int i;
  emxArray_real_T *buf_val;
  emxArray_int32_T *buf_indx;
  emxArray_int32_T *A;
  emxArray_real_T *b_A;
  boolean_T ascend;
  int i1;
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
  i0 = row_ptr->size[0] - 1;
  i = 1;
  emxInit_real_T(&buf_val, 1);
  emxInit_int32_T(&buf_indx, 1);
  emxInit_int32_T(&A, 1);
  emxInit_real_T(&b_A, 1);
  while (i <= i0) {
    ascend = true;
    i1 = row_ptr->data[i] - 1;
    j = row_ptr->data[i - 1];
    exitg3 = false;
    while ((!exitg3) && (j + 1 <= i1)) {
      if (col_ind->data[j] < col_ind->data[j - 1]) {
        ascend = false;
        exitg3 = true;
      } else {
        j++;
      }
    }

    if (!ascend) {
      i1 = A->size[0];
      A->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)A, i1, (int)sizeof(int));
      loop_ub = row_ptr->data[i] - row_ptr->data[i - 1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        A->data[i1] = 0;
      }

      i1 = buf_indx->size[0];
      buf_indx->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)buf_indx, i1, (int)sizeof(int));
      i1 = b_A->size[0];
      b_A->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)b_A, i1, (int)sizeof(double));
      loop_ub = row_ptr->data[i] - row_ptr->data[i - 1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_A->data[i1] = 0.0;
      }

      i1 = buf_val->size[0];
      buf_val->size[0] = b_A->size[0];
      emxEnsureCapacity((emxArray__common *)buf_val, i1, (int)sizeof(double));
      ind = 1U;
      i1 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i1; j++) {
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
      i1 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i1; j++) {
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

void crs_sort0(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind)
{
  int i2;
  int i;
  emxArray_int32_T *buf_indx;
  emxArray_int32_T *A;
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
  int exitg2;
  int b_i;
  boolean_T gt;
  boolean_T guard2 = false;
  i2 = row_ptr->size[0] - 1;
  i = 1;
  emxInit_int32_T(&buf_indx, 1);
  emxInit_int32_T(&A, 1);
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
      ind = 1U;
      i3 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i3; j++) {
        buf_indx->data[(int)ind - 1] = col_ind->data[j - 1];
        ind++;
      }

      n = buf_indx->size[0];
      if (n <= 1) {
      } else {
        l = (int)((unsigned int)n >> 1U) + 1;
        ir = n;
        do {
          exitg1 = 0;
          guard1 = false;
          if (l <= 1) {
            r0 = buf_indx->data[ir - 1];
            buf_indx->data[ir - 1] = buf_indx->data[0];
            ir--;
            if (ir == 1) {
              exitg1 = 1;
            } else {
              guard1 = true;
            }
          } else {
            l--;
            r0 = buf_indx->data[l - 1];
            guard1 = true;
          }

          if (guard1) {
            j = l;
            do {
              exitg2 = 0;
              b_i = j - 1;
              j <<= 1;
              gt = false;
              guard2 = false;
              if (j >= ir) {
                if (j == ir) {
                  gt = true;
                  guard2 = true;
                } else if (j > ir) {
                  exitg2 = 1;
                } else {
                  guard2 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                if ((!gt) && (buf_indx->data[j - 1] < buf_indx->data[j])) {
                  j++;
                }

                if (r0 >= buf_indx->data[j - 1]) {
                  exitg2 = 1;
                } else {
                  buf_indx->data[b_i] = buf_indx->data[j - 1];
                }
              }
            } while (exitg2 == 0);

            buf_indx->data[b_i] = r0;
          }
        } while (exitg1 == 0);

        buf_indx->data[0] = r0;
      }

      ind = 1U;
      i3 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i3; j++) {
        col_ind->data[j - 1] = buf_indx->data[(int)ind - 1];
        ind++;
      }
    }

    i++;
  }

  emxFree_int32_T(&A);
  emxFree_int32_T(&buf_indx);
}

void crs_sort_initialize(void)
{
}

void crs_sort_terminate(void)
{
}

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}
