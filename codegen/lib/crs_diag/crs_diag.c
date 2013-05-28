#include "crs_diag.h"
#include "palc.h"

void crs_diag(const struct_T *A, plcArray_real_T *D)
{
  int32_T i0;
  int32_T loop_ub;
  int32_T j;
  boolean_T exitg1;
  i0 = D->size[0];
  D->size[0] = A->nrows;
  plcEnsureCapacity((plcArray__common *)D, i0, (int32_T)sizeof(real_T));
  loop_ub = A->nrows;
  for (i0 = 0; i0 < loop_ub; i0++) {
    D->data[i0] = 0.0;
  }

  for (loop_ub = 1; loop_ub <= A->nrows; loop_ub++) {
    i0 = A->row_ptr->data[loop_ub] - 1;
    j = A->row_ptr->data[loop_ub - 1];
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (j <= i0)) {
      if (A->col_ind->data[j - 1] == loop_ub) {
        D->data[loop_ub - 1] = A->val->data[j - 1];
        exitg1 = TRUE;
      } else {
        j++;
      }
    }
  }
}

void crs_diag1(const struct_T *A, int32_T k, plcArray_real_T *D)
{
  int32_T y;
  int32_T i1;
  int32_T j;
  boolean_T exitg2;
  boolean_T exitg1;
  if (k < 0) {
    y = -k;
  } else {
    y = k;
  }

  i1 = D->size[0];
  D->size[0] = A->nrows - y;
  plcEnsureCapacity((plcArray__common *)D, i1, (int32_T)sizeof(real_T));
  y = A->nrows - y;
  for (i1 = 0; i1 < y; i1++) {
    D->data[i1] = 0.0;
  }

  if (k >= 0) {
    for (y = 1; y <= A->nrows; y++) {
      i1 = A->row_ptr->data[y] - 1;
      j = A->row_ptr->data[y - 1];
      exitg2 = FALSE;
      while ((exitg2 == FALSE) && (j <= i1)) {
        if (A->col_ind->data[j - 1] == y + k) {
          D->data[y - 1] = A->val->data[j - 1];
          exitg2 = TRUE;
        } else {
          j++;
        }
      }
    }
  } else {
    for (y = 1; y <= A->nrows; y++) {
      i1 = A->row_ptr->data[y] - 1;
      j = A->row_ptr->data[y - 1];
      exitg1 = FALSE;
      while ((exitg1 == FALSE) && (j <= i1)) {
        if (A->col_ind->data[j - 1] == y + k) {
          D->data[(y + k) - 1] = A->val->data[j - 1];
          exitg1 = TRUE;
        } else {
          j++;
        }
      }
    }
  }
}

void crs_diag_initialize(void)
{
}

void crs_diag_terminate(void)
{
}

