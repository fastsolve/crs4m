#include "crs_prodAx.h"
#include "omp.h"
#include "m2c.h"

static void b_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, const
  emxArray_real_T *x, emxArray_real_T *b);
static void c_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, const
  emxArray_real_T *x, emxArray_real_T *b);
static void crs_prodAx_kernel(const int32_T *row_ptr, const
  int32_T *col_ind, const real_T *val, const real_T *x, int x_m, real_T *b, int b_m, int nrows, int nrhs, boolean_T ismt);
static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static void m2c_error(void);
static void m2c_warn(void);
static double rt_roundd(double u);
static void b_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, const
  emxArray_real_T *x, emxArray_real_T *b)
{
  int i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A_nrows;
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int)sizeof(double));
  crs_prodAx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, A_nrows, A_nrows, x->size[1], false);
}

static void c_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, const
  emxArray_real_T *x, emxArray_real_T *b)
{
  int i3;
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  i3 = b->size[0];
  crs_prodAx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i3, A_nrows, x->size[1], false);
}

static void crs_prodAx_kernel(const int32_T *row_ptr, const
  int32_T *col_ind, const real_T *val, const real_T *x, int x_m, real_T *b, int b_m, int nrows, int nrhs, boolean_T ismt)
{
  int istart;
  int nthreads;
  int iend;
  int thr_id;
  int xoffset;
  double y;
  int boffset;
  int k;
  double b_y;
  int chunk;
  int i;
  int b_remainder;
  int c_remainder;
  double t;
  int i2;
  int j;

  if (ismt) {
    nthreads = omp_get_num_threads();
    if (nthreads == 1) {
      istart = 0;
      iend = nrows;
    } else {
      thr_id = omp_get_thread_num();
      y = (double)nrows / (double)nthreads;
      if (y < 0.0) {
        b_y = ceil(y);
      } else {
        b_y = floor(y);
      }

      chunk = (int)rt_roundd(b_y);
      b_remainder = nrows - nthreads * chunk;
      if (b_remainder <= thr_id) {
        c_remainder = b_remainder;
      } else {
        c_remainder = thr_id;
      }

      istart = thr_id * chunk + c_remainder;
      iend = (istart + chunk) + (thr_id < b_remainder);
    }
  } else {
    istart = 0;
    iend = nrows;
  }

  xoffset = -1;
  boffset = -1;
  for (k = 1; k <= nrhs; k++) {
    for (i = istart + 1; i <= iend; i++) {
      t = 0.0;
      i2 = row_ptr[i] - 1;
      for (j = row_ptr[i - 1]; j <= i2; j++) {
        t += val[j - 1] * x[xoffset + col_ind[j - 1]];
      }

      b[boffset + i] = t;
    }

    xoffset += x_m;
    boffset += b_m;
  }
}

static void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_int32_T(&pStruct->row_ptr);
  emxFree_int32_T(&pStruct->col_ind);
  emxFree_real_T(&pStruct->val);
}

static void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_int32_T(&pStruct->row_ptr, 1);
  emxInit_int32_T(&pStruct->col_ind, 1);
  emxInit_real_T(&pStruct->val, 1);
}

static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void m2c_error(void)
{
  M2C_error("crs_prodAx:BufferTooSmal",
            "Buffer space for output b is too small.");
}

static void m2c_warn(void)
{
  M2C_warn("crs_prodAx:NestedParallel",
           "You are trying to use nested parallel regions. Solution may be incorrect.");
}

static double rt_roundd(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

void crs_prodAx(const struct0_T *A, const emxArray_real_T *x, emxArray_real_T *b,
                const emxArray_int32_T *nthreads)
{
  int n;
  boolean_T b0;
  int i1;
  if ((b->size[0] < A->nrows) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  n = omp_get_num_threads();
  b0 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b0 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      m2c_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    n = omp_get_num_threads();
    i1 = b->size[0];
    crs_prodAx_kernel(A->row_ptr->data, A->col_ind->data, A->val->data, x->data, x->size[0], b->data, i1, A->nrows, x->size[1], n > 1);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    i1 = b->size[0];
    crs_prodAx_kernel(A->row_ptr->data, A->col_ind->data, A->val->data, x->data, x->size[0], b->data, i1, A->nrows, x->size[1], b0);
  }
}

void crs_prodAx_initialize(void)
{
}

void crs_prodAx_ser(const struct0_T *A, const emxArray_real_T *x,
                    emxArray_real_T *b)
{
  b_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b);
}

void crs_prodAx_ser1(const struct0_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  c_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b);
}

void crs_prodAx_terminate(void)
{
}

void emxDestroy_struct0_T(struct0_T emxArray)
{
  emxFreeStruct_struct0_T(&emxArray);
}

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T1(pEmxArray, numDimensions);
}

void emxInit_struct0_T(struct0_T *pStruct)
{
  emxInitStruct_struct0_T(pStruct);
}
