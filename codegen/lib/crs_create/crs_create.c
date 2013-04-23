#include "crs_create.h"
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

void crs_create(const emxArray_int32_T *is, const emxArray_int32_T *js, const
                emxArray_real_T *vs, struct_T *A)
{
  int32_T mtmp;
  int32_T j;
  int32_T b_mtmp;
  int32_T i0;
  boolean_T ascend;
  boolean_T exitg1;
  emxArray_int32_T *x;
  mtmp = is->data[0];
  if ((is->size[0] > 1) && (1 < is->size[0])) {
    for (j = 1; j + 1 <= is->size[0]; j++) {
      if (is->data[j] > mtmp) {
        mtmp = is->data[j];
      }
    }
  }

  b_mtmp = js->data[0];
  if ((js->size[0] > 1) && (1 < js->size[0])) {
    for (j = 1; j + 1 <= js->size[0]; j++) {
      if (js->data[j] > b_mtmp) {
        b_mtmp = js->data[j];
      }
    }
  }

  i0 = A->row_ptr->size[0];
  A->row_ptr->size[0] = mtmp + 1;
  emxEnsureCapacity((emxArray__common *)A->row_ptr, i0, (int32_T)sizeof(int32_T));
  A->nrows = mtmp;
  A->ncols = b_mtmp;
  i0 = is->size[0];
  for (b_mtmp = 0; b_mtmp + 1 <= i0; b_mtmp++) {
    A->row_ptr->data[is->data[b_mtmp]]++;
  }

  A->row_ptr->data[0] = 1;
  for (b_mtmp = 1; b_mtmp <= mtmp; b_mtmp++) {
    A->row_ptr->data[b_mtmp] += A->row_ptr->data[b_mtmp - 1];
  }

  ascend = TRUE;
  b_mtmp = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (b_mtmp <= is->size[0] - 2)) {
    if (is->data[b_mtmp + 1] < is->data[b_mtmp]) {
      ascend = FALSE;
      exitg1 = TRUE;
    } else {
      b_mtmp++;
    }
  }

  if (ascend) {
    i0 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->col_ind, i0, (int32_T)sizeof
                      (int32_T));
    j = js->size[0];
    for (i0 = 0; i0 < j; i0++) {
      A->col_ind->data[i0] = js->data[i0];
    }

    i0 = A->val->size[0];
    A->val->size[0] = vs->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i0, (int32_T)sizeof(real_T));
    j = vs->size[0];
    for (i0 = 0; i0 < j; i0++) {
      A->val->data[i0] = vs->data[i0];
    }
  } else {
    i0 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->col_ind, i0, (int32_T)sizeof
                      (int32_T));
    i0 = A->val->size[0];
    A->val->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i0, (int32_T)sizeof(real_T));
    for (b_mtmp = 0; b_mtmp < is->size[0]; b_mtmp++) {
      j = A->row_ptr->data[is->data[b_mtmp] - 1];
      A->val->data[A->row_ptr->data[is->data[b_mtmp] - 1] - 1] = vs->data[b_mtmp];
      A->col_ind->data[j - 1] = js->data[b_mtmp];
      A->row_ptr->data[is->data[b_mtmp] - 1]++;
    }

    emxInit_int32_T(&x, 1);
    i0 = x->size[0];
    x->size[0] = A->row_ptr->size[0];
    emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)sizeof(int32_T));
    j = A->row_ptr->size[0];
    for (i0 = 0; i0 < j; i0++) {
      x->data[i0] = A->row_ptr->data[i0];
    }

    i0 = A->row_ptr->size[0];
    for (b_mtmp = 0; b_mtmp <= i0 - 2; b_mtmp++) {
      j = (x->size[0] - b_mtmp) - 2;
      A->row_ptr->data[j + 1] = A->row_ptr->data[j];
    }

    emxFree_int32_T(&x);
    A->row_ptr->data[0] = 1;
  }
}

void crs_create1(const emxArray_int32_T *is, const emxArray_int32_T *js, const
                 emxArray_real_T *vs, int32_T ni, int32_T nj, struct_T *A)
{
  int32_T i1;
  int32_T i;
  boolean_T ascend;
  boolean_T exitg1;
  int32_T j;
  emxArray_int32_T *x;
  i1 = A->row_ptr->size[0];
  A->row_ptr->size[0] = ni + 1;
  emxEnsureCapacity((emxArray__common *)A->row_ptr, i1, (int32_T)sizeof(int32_T));
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
    emxEnsureCapacity((emxArray__common *)A->col_ind, i1, (int32_T)sizeof
                      (int32_T));
    j = js->size[0];
    for (i1 = 0; i1 < j; i1++) {
      A->col_ind->data[i1] = js->data[i1];
    }

    i1 = A->val->size[0];
    A->val->size[0] = vs->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i1, (int32_T)sizeof(real_T));
    j = vs->size[0];
    for (i1 = 0; i1 < j; i1++) {
      A->val->data[i1] = vs->data[i1];
    }
  } else {
    i1 = A->col_ind->size[0];
    A->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->col_ind, i1, (int32_T)sizeof
                      (int32_T));
    i1 = A->val->size[0];
    A->val->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)A->val, i1, (int32_T)sizeof(real_T));
    for (i = 0; i < is->size[0]; i++) {
      j = A->row_ptr->data[is->data[i] - 1];
      A->val->data[A->row_ptr->data[is->data[i] - 1] - 1] = vs->data[i];
      A->col_ind->data[j - 1] = js->data[i];
      A->row_ptr->data[is->data[i] - 1]++;
    }

    emxInit_int32_T(&x, 1);
    i1 = x->size[0];
    x->size[0] = A->row_ptr->size[0];
    emxEnsureCapacity((emxArray__common *)x, i1, (int32_T)sizeof(int32_T));
    j = A->row_ptr->size[0];
    for (i1 = 0; i1 < j; i1++) {
      x->data[i1] = A->row_ptr->data[i1];
    }

    i1 = A->row_ptr->size[0];
    for (i = 0; i <= i1 - 2; i++) {
      j = (x->size[0] - i) - 2;
      A->row_ptr->data[j + 1] = A->row_ptr->data[j];
    }

    emxFree_int32_T(&x);
    A->row_ptr->data[0] = 1;
  }
}

void crs_create_initialize(void)
{
}

void crs_create_terminate(void)
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
