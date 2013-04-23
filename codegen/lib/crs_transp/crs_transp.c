#include "crs_transp.h"
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

void crs_transp(const struct_T A, struct_T *At)
{
  emxArray_int32_T *js;
  uint32_T unnamed_idx_0;
  int32_T i0;
  int32_T nrows;
  int32_T i;
  int32_T j;
  int32_T mtmp;
  boolean_T ascend;
  boolean_T exitg1;
  emxInit_int32_T(&js, 1);
  unnamed_idx_0 = (uint32_T)A.col_ind->size[0];
  i0 = js->size[0];
  js->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)js, i0, (int32_T)sizeof(int32_T));
  nrows = A.row_ptr->size[0] - 1;
  for (i = 1; i <= nrows; i++) {
    i0 = A.row_ptr->data[i] - 1;
    for (j = A.row_ptr->data[i - 1]; j <= i0; j++) {
      js->data[j - 1] = i;
    }
  }

  mtmp = A.col_ind->data[0];
  if ((A.col_ind->size[0] > 1) && (1 < A.col_ind->size[0])) {
    for (nrows = 1; nrows + 1 <= A.col_ind->size[0]; nrows++) {
      if (A.col_ind->data[nrows] > mtmp) {
        mtmp = A.col_ind->data[nrows];
      }
    }
  }

  j = js->data[0];
  if ((js->size[0] > 1) && (1 < js->size[0])) {
    for (nrows = 1; nrows + 1 <= js->size[0]; nrows++) {
      if (js->data[nrows] > j) {
        j = js->data[nrows];
      }
    }
  }

  i0 = At->row_ptr->size[0];
  At->row_ptr->size[0] = mtmp + 1;
  emxEnsureCapacity((emxArray__common *)At->row_ptr, i0, (int32_T)sizeof(int32_T));
  At->nrows = mtmp;
  At->ncols = j;
  i0 = A.col_ind->size[0];
  for (i = 0; i + 1 <= i0; i++) {
    At->row_ptr->data[A.col_ind->data[i]]++;
  }

  At->row_ptr->data[0] = 1;
  for (i = 1; i <= mtmp; i++) {
    At->row_ptr->data[i] += At->row_ptr->data[i - 1];
  }

  ascend = TRUE;
  i = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i <= A.col_ind->size[0] - 2)) {
    if (A.col_ind->data[i + 1] < A.col_ind->data[i]) {
      ascend = FALSE;
      exitg1 = TRUE;
    } else {
      i++;
    }
  }

  if (ascend) {
    i0 = At->col_ind->size[0];
    At->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)At->col_ind, i0, (int32_T)sizeof
                      (int32_T));
    nrows = js->size[0];
    for (i0 = 0; i0 < nrows; i0++) {
      At->col_ind->data[i0] = js->data[i0];
    }

    i0 = At->val->size[0];
    At->val->size[0] = A.val->size[0];
    emxEnsureCapacity((emxArray__common *)At->val, i0, (int32_T)sizeof(real_T));
    nrows = A.val->size[0];
    for (i0 = 0; i0 < nrows; i0++) {
      At->val->data[i0] = A.val->data[i0];
    }
  } else {
    i0 = At->col_ind->size[0];
    At->col_ind->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)At->col_ind, i0, (int32_T)sizeof
                      (int32_T));
    i0 = At->val->size[0];
    At->val->size[0] = js->size[0];
    emxEnsureCapacity((emxArray__common *)At->val, i0, (int32_T)sizeof(real_T));
    for (i = 0; i < A.col_ind->size[0]; i++) {
      j = At->row_ptr->data[A.col_ind->data[i] - 1];
      At->val->data[At->row_ptr->data[A.col_ind->data[i] - 1] - 1] = A.val->
        data[i];
      At->col_ind->data[j - 1] = js->data[i];
      At->row_ptr->data[A.col_ind->data[i] - 1]++;
    }

    i0 = js->size[0];
    js->size[0] = At->row_ptr->size[0];
    emxEnsureCapacity((emxArray__common *)js, i0, (int32_T)sizeof(int32_T));
    nrows = At->row_ptr->size[0];
    for (i0 = 0; i0 < nrows; i0++) {
      js->data[i0] = At->row_ptr->data[i0];
    }

    i0 = At->row_ptr->size[0];
    for (i = 0; i <= i0 - 2; i++) {
      nrows = (js->size[0] - i) - 2;
      At->row_ptr->data[nrows + 1] = At->row_ptr->data[nrows];
    }

    At->row_ptr->data[0] = 1;
  }

  emxFree_int32_T(&js);
}

void crs_transp_initialize(void)
{
}

void crs_transp_terminate(void)
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
