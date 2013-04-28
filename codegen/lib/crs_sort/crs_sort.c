#include "crs_sort.h"
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
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);





void crs_sort(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind,
              emxArray_real_T *val)
{
  int32_T i0;
  int32_T i;
  emxArray_real_T *buf_val;
  emxArray_int32_T *buf_indx;
  boolean_T ascend;
  int32_T j;
  int32_T n;
  boolean_T exitg3;
  uint32_T ind;
  int32_T l;
  int32_T exitg1;
  boolean_T guard1 = FALSE;
  int32_T r0;
  real_T t0;
  int32_T exitg2;
  int32_T b_i;
  boolean_T guard2 = FALSE;
  i0 = row_ptr->size[0] - 1;
  i = 1;
  emxInit_real_T(&buf_val, 1);
  emxInit_int32_T(&buf_indx, 1);
  while (i <= i0) {
    ascend = TRUE;
    j = row_ptr->data[i - 1];
    n = row_ptr->data[i] - 1;
    exitg3 = FALSE;
    while ((exitg3 == FALSE) && (j + 1 <= n)) {
      if (col_ind->data[j] < col_ind->data[j - 1]) {
        ascend = FALSE;
        exitg3 = TRUE;
      } else {
        j++;
      }
    }

    if (!ascend) {
      n = buf_indx->size[0];
      buf_indx->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)buf_indx, n, (int32_T)sizeof(int32_T));
      n = buf_val->size[0];
      buf_val->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)buf_val, n, (int32_T)sizeof(real_T));
      ind = 1U;
      n = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= n; j++) {
        buf_indx->data[(int32_T)ind - 1] = col_ind->data[j - 1];
        buf_val->data[(int32_T)ind - 1] = val->data[j - 1];
        ind++;
      }

      n = buf_indx->size[0];
      if (n <= 1) {
      } else {
        l = (int32_T)((uint32_T)n >> 1U);
        do {
          exitg1 = 0;
          guard1 = FALSE;
          if (l + 1 <= 1) {
            r0 = buf_indx->data[n - 1];
            t0 = buf_val->data[n - 1];
            buf_indx->data[n - 1] = buf_indx->data[0];
            buf_val->data[n - 1] = buf_val->data[0];
            n--;
            if (n == 1) {
              exitg1 = 1;
            } else {
              guard1 = TRUE;
            }
          } else {
            l--;
            r0 = buf_indx->data[l];
            t0 = buf_val->data[l];
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            j = l;
            do {
              exitg2 = 0;
              b_i = j;
              j = ((j + 1) << 1) - 1;
              ascend = FALSE;
              guard2 = FALSE;
              if (j + 1 >= n) {
                if (j + 1 == n) {
                  ascend = TRUE;
                  guard2 = TRUE;
                } else if (j + 1 > n) {
                  exitg2 = 1;
                } else {
                  guard2 = TRUE;
                }
              } else {
                guard2 = TRUE;
              }

              if (guard2 == TRUE) {
                if ((!ascend) && (buf_indx->data[j] < buf_indx->data[j + 1])) {
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
      n = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= n; j++) {
        col_ind->data[j - 1] = buf_indx->data[(int32_T)ind - 1];
        val->data[j - 1] = buf_val->data[(int32_T)ind - 1];
        ind++;
      }
    }

    i++;
  }

  emxFree_int32_T(&buf_indx);
  emxFree_real_T(&buf_val);
}

void crs_sort0(const emxArray_int32_T *row_ptr, emxArray_int32_T *col_ind)
{
  int32_T i1;
  int32_T i;
  emxArray_int32_T *buf_indx;
  boolean_T ascend;
  int32_T j;
  int32_T n;
  boolean_T exitg3;
  uint32_T ind;
  int32_T l;
  int32_T exitg1;
  boolean_T guard1 = FALSE;
  int32_T r0;
  int32_T exitg2;
  int32_T b_i;
  boolean_T guard2 = FALSE;
  i1 = row_ptr->size[0] - 1;
  i = 1;
  emxInit_int32_T(&buf_indx, 1);
  while (i <= i1) {
    ascend = TRUE;
    j = row_ptr->data[i - 1];
    n = row_ptr->data[i] - 1;
    exitg3 = FALSE;
    while ((exitg3 == FALSE) && (j + 1 <= n)) {
      if (col_ind->data[j] < col_ind->data[j - 1]) {
        ascend = FALSE;
        exitg3 = TRUE;
      } else {
        j++;
      }
    }

    if (!ascend) {
      n = buf_indx->size[0];
      buf_indx->size[0] = row_ptr->data[i] - row_ptr->data[i - 1];
      emxEnsureCapacity((emxArray__common *)buf_indx, n, (int32_T)sizeof(int32_T));
      ind = 1U;
      n = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= n; j++) {
        buf_indx->data[(int32_T)ind - 1] = col_ind->data[j - 1];
        ind++;
      }

      n = buf_indx->size[0];
      if (n <= 1) {
      } else {
        l = (int32_T)((uint32_T)n >> 1U) + 1;
        do {
          exitg1 = 0;
          guard1 = FALSE;
          if (l <= 1) {
            r0 = buf_indx->data[n - 1];
            buf_indx->data[n - 1] = buf_indx->data[0];
            n--;
            if (n == 1) {
              exitg1 = 1;
            } else {
              guard1 = TRUE;
            }
          } else {
            l--;
            r0 = buf_indx->data[l - 1];
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            j = l;
            do {
              exitg2 = 0;
              b_i = j - 1;
              j <<= 1;
              ascend = FALSE;
              guard2 = FALSE;
              if (j >= n) {
                if (j == n) {
                  ascend = TRUE;
                  guard2 = TRUE;
                } else if (j > n) {
                  exitg2 = 1;
                } else {
                  guard2 = TRUE;
                }
              } else {
                guard2 = TRUE;
              }

              if (guard2 == TRUE) {
                if ((!ascend) && (buf_indx->data[j - 1] < buf_indx->data[j])) {
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
      n = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= n; j++) {
        col_ind->data[j - 1] = buf_indx->data[(int32_T)ind - 1];
        ind++;
      }
    }

    i++;
  }

  emxFree_int32_T(&buf_indx);
}

void crs_sort_initialize(void)
{
}

void crs_sort_terminate(void)
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




