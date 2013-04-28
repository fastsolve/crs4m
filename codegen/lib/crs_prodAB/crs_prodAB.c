#include "crs_prodAB.h"
#include "stdio.h"
#ifdef BUILD_MEX

#include "mex.h"
#define malloc                         mxMalloc
#define calloc                         mxCalloc
#define realloc                        mxRealloc
#define free                           mxFree
#define emlrtIsMATLABThread(s)         1
#else
#define emlrtIsMATLABThread(s)         0
#define mexErrMsgIdAndTxt(a,b)
#endif

#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif

#ifndef typedef_emxArray__common
#define typedef_emxArray__common

typedef struct emxArray__common emxArray__common;

#endif

static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
static void msg_error(void);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize)
{
  int32_T newNumel;
  int32_T i;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    if (emxArray->allocatedSize==0)
      i = newNumel;
    else {
      i = emxArray->allocatedSize;
      if (i < 16) {
        i = 16;
      }

      while (i < newNumel) {
        i <<= 1;
      }
    }

    newData = calloc((uint32_T)i, (uint32_T)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (uint32_T)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = TRUE;
  }
}


static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions)
{
  emxArray_int32_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int32_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
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

void crs_prodAB(const struct_T *A, const struct_T *B, struct_T *C)
{
  emxArray_int32_T *b_index;
  int32_T i0;
  int32_T u0;
  int32_T u1;
  int32_T i;
  int32_T istart;
  int32_T clength;
  int32_T k;
  emxArray_real_T *temp;
  if (A->ncols != B->nrows) {
    msg_error();
  }

  emxInit_int32_T(&b_index, 1);
  i0 = C->row_ptr->size[0];
  C->row_ptr->size[0] = 0;
  emxEnsureCapacity((emxArray__common *)C->row_ptr, i0, (int32_T)sizeof(int32_T));
  i0 = C->col_ind->size[0];
  C->col_ind->size[0] = 0;
  emxEnsureCapacity((emxArray__common *)C->col_ind, i0, (int32_T)sizeof(int32_T));
  i0 = C->val->size[0];
  C->val->size[0] = 0;
  emxEnsureCapacity((emxArray__common *)C->val, i0, (int32_T)sizeof(real_T));
  C->nrows = A->nrows;
  C->ncols = B->ncols;
  i0 = C->row_ptr->size[0];
  C->row_ptr->size[0] = A->row_ptr->size[0];
  emxEnsureCapacity((emxArray__common *)C->row_ptr, i0, (int32_T)sizeof(int32_T));
  C->row_ptr->data[0] = 1;
  u0 = A->ncols;
  u1 = B->ncols;
  if (u0 >= u1) {
  } else {
    u0 = u1;
  }

  i0 = b_index->size[0];
  b_index->size[0] = u0;
  emxEnsureCapacity((emxArray__common *)b_index, i0, (int32_T)sizeof(int32_T));
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
  emxEnsureCapacity((emxArray__common *)C->col_ind, i0, (int32_T)sizeof(int32_T));
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

  emxInit_real_T(&temp, 1);
  u0 = C->row_ptr->data[A->nrows] - 1;
  i0 = C->val->size[0];
  C->val->size[0] = u0;
  emxEnsureCapacity((emxArray__common *)C->val, i0, (int32_T)sizeof(real_T));
  i0 = temp->size[0];
  temp->size[0] = b_index->size[0];
  emxEnsureCapacity((emxArray__common *)temp, i0, (int32_T)sizeof(real_T));
  u0 = b_index->size[0];
  emxFree_int32_T(&b_index);
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

  emxFree_real_T(&temp);
}

void crs_prodAB_initialize(void)
{
}

void crs_prodAB_terminate(void)
{
}

emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions,
  int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreateWrapper_int32_T(int32_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_int32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_int32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols)
{
  emxArray_int32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_int32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_int32_T(emxArray_int32_T *emxArray)
{
  emxFree_int32_T(&emxArray);
}

void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}
