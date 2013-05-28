#include "crs_prodAB.h"
#include "palc.h"

static void msg_error(void);

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

void crs_prodAB(const struct_T *A, const struct_T *B, struct_T *C)
{
  plcArray_int32_T *b_index;
  int32_T i0;
  int32_T u0;
  int32_T u1;
  int32_T i;
  int32_T istart;
  int32_T clength;
  int32_T k;
  plcArray_real_T *temp;
  if (A->ncols != B->nrows) {
    msg_error();
  }

  plcInit_int32_T(&b_index, 1);
  i0 = C->row_ptr->size[0];
  C->row_ptr->size[0] = 0;
  plcEnsureCapacity((plcArray__common *)C->row_ptr, i0, (int32_T)sizeof(int32_T));
  i0 = C->col_ind->size[0];
  C->col_ind->size[0] = 0;
  plcEnsureCapacity((plcArray__common *)C->col_ind, i0, (int32_T)sizeof(int32_T));
  i0 = C->val->size[0];
  C->val->size[0] = 0;
  plcEnsureCapacity((plcArray__common *)C->val, i0, (int32_T)sizeof(real_T));
  C->nrows = A->nrows;
  C->ncols = B->ncols;
  i0 = C->row_ptr->size[0];
  C->row_ptr->size[0] = A->row_ptr->size[0];
  plcEnsureCapacity((plcArray__common *)C->row_ptr, i0, (int32_T)sizeof(int32_T));
  C->row_ptr->data[0] = 1;
  u0 = A->ncols;
  u1 = B->ncols;
  if (u0 >= u1) {
  } else {
    u0 = u1;
  }

  i0 = b_index->size[0];
  b_index->size[0] = u0;
  plcEnsureCapacity((plcArray__common *)b_index, i0, (int32_T)sizeof(int32_T));
  for (i0 = 0; i0 < u0; i0++) {
    b_index->data[i0] = 0;
  }

  for (i = 1; i <= A->nrows; i++) {
    istart = -1;
    clength = 0;
    i0 = A->row_ptr->data[i] - 1;
    for (u0 = A->row_ptr->data[i - 1]; u0 <= i0; u0++) {
      u1 = B->row_ptr->data[A->col_ind->data[u0 - 1]] - 1;
      for (k = B->row_ptr->data[A->col_ind->data[u0 - 1] - 1] - 1; k + 1 <= u1;
           k++) {
        if (b_index->data[B->col_ind->data[k] - 1] == 0) {
          b_index->data[B->col_ind->data[k] - 1] = istart;
          istart = B->col_ind->data[k];
          clength++;
        }
      }
    }

    C->row_ptr->data[i] = C->row_ptr->data[i - 1] + clength;
    i0 = C->row_ptr->data[i] - 1;
    for (u0 = C->row_ptr->data[i - 1]; u0 <= i0; u0++) {
      k = istart;
      istart = b_index->data[istart - 1];
      b_index->data[k - 1] = 0;
    }

    b_index->data[i - 1] = 0;
  }

  u0 = C->row_ptr->data[A->nrows] - 1;
  i0 = C->col_ind->size[0];
  C->col_ind->size[0] = u0;
  plcEnsureCapacity((plcArray__common *)C->col_ind, i0, (int32_T)sizeof(int32_T));
  for (i = 1; i <= A->nrows; i++) {
    istart = -1;
    clength = 0;
    i0 = A->row_ptr->data[i] - 1;
    for (u0 = A->row_ptr->data[i - 1]; u0 <= i0; u0++) {
      u1 = B->row_ptr->data[A->col_ind->data[u0 - 1]] - 1;
      for (k = B->row_ptr->data[A->col_ind->data[u0 - 1] - 1] - 1; k + 1 <= u1;
           k++) {
        if (b_index->data[B->col_ind->data[k] - 1] == 0) {
          b_index->data[B->col_ind->data[k] - 1] = istart;
          istart = B->col_ind->data[k];
          clength++;
        }
      }
    }

    C->row_ptr->data[i] = C->row_ptr->data[i - 1] + clength;
    i0 = C->row_ptr->data[i] - 1;
    for (u0 = C->row_ptr->data[i - 1]; u0 <= i0; u0++) {
      C->col_ind->data[u0 - 1] = istart;
      istart = b_index->data[istart - 1];
      b_index->data[C->col_ind->data[u0 - 1] - 1] = 0;
    }

    b_index->data[i - 1] = 0;
  }

  plcInit_real_T(&temp, 1);
  u0 = C->row_ptr->data[A->nrows] - 1;
  i0 = C->val->size[0];
  C->val->size[0] = u0;
  plcEnsureCapacity((plcArray__common *)C->val, i0, (int32_T)sizeof(real_T));
  i0 = temp->size[0];
  temp->size[0] = b_index->size[0];
  plcEnsureCapacity((plcArray__common *)temp, i0, (int32_T)sizeof(real_T));
  u0 = b_index->size[0];
  plcFree_int32_T(&b_index);
  for (i0 = 0; i0 < u0; i0++) {
    temp->data[i0] = 0.0;
  }

  for (i = 1; i <= A->nrows; i++) {
    i0 = A->row_ptr->data[i] - 1;
    for (u0 = A->row_ptr->data[i - 1] - 1; u0 + 1 <= i0; u0++) {
      u1 = B->row_ptr->data[A->col_ind->data[u0]] - 1;
      for (k = B->row_ptr->data[A->col_ind->data[u0] - 1] - 1; k + 1 <= u1; k++)
      {
        temp->data[B->col_ind->data[k] - 1] += A->val->data[u0] * B->val->data[k];
      }
    }

    i0 = C->row_ptr->data[i] - 1;
    for (u0 = C->row_ptr->data[i - 1] - 1; u0 + 1 <= i0; u0++) {
      C->val->data[u0] = temp->data[C->col_ind->data[u0] - 1];
      temp->data[C->col_ind->data[u0] - 1] = 0.0;
    }
  }

  plcFree_real_T(&temp);
}

void crs_prodAB_initialize(void)
{
}

void crs_prodAB_terminate(void)
{
}

