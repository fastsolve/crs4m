#include "crs_rowind.h"
#include "m2c.h"
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

extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);



void crs_rowind(const emxArray_int32_T *row_ptr, const emxArray_int32_T *col_ind,
                emxArray_int32_T *row_ind)
{
  uint32_T unnamed_idx_0;
  int32_T i0;
  int32_T nrows;
  int32_T i;
  int32_T j;
  unnamed_idx_0 = (uint32_T)col_ind->size[0];
  i0 = row_ind->size[0];
  row_ind->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)row_ind, i0, (int32_T)sizeof(int32_T));
  nrows = row_ptr->size[0] - 1;
  for (i = 1; i <= nrows; i++) {
    i0 = row_ptr->data[i] - 1;
    for (j = row_ptr->data[i - 1]; j <= i0; j++) {
      row_ind->data[j - 1] = i;
    }
  }
}

void crs_rowind_initialize(void)
{
}

void crs_rowind_terminate(void)
{
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


