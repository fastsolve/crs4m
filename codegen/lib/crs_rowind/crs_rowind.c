#include "crs_rowind.h"
#include "m2c.h"

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

