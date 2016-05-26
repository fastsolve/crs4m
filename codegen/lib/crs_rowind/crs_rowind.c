#include "crs_rowind.h"
#include "m2c.h"

void crs_rowind(const emxArray_int32_T *row_ptr, const emxArray_int32_T *col_ind,
                emxArray_int32_T *row_ind)
{
  unsigned int unnamed_idx_0;
  int i0;
  int nrows;
  int i;
  int j;
  unnamed_idx_0 = (unsigned int)col_ind->size[0];
  i0 = row_ind->size[0];
  row_ind->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)row_ind, i0, (int)sizeof(int));
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

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}
