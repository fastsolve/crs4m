#include "crs_prodAtx.h"
#include "omp.h"
#include "m2c.h"

static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, int
  A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, int
  A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void crs_prodAtx_kernel(const int32_T *row_ptr, const
  int32_T *col_ind, const real_T *val, const real_T *x, int x_m, real_T *b, int b_m, int nrows, int ncols, int nrhs,
  boolean_T ismt);
static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static void get_local_chunk(int m, int varargin_2, int *istart, int *iend);
static void m2c_error(void);
static void m2c_warn(void);
static double rt_roundd(double u);
static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, int
  A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A_ncols;
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int)sizeof(double));
  crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, A_ncols, A_nrows, A_ncols, x->size[1], false);
}

static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int A_nrows, int
  A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int i3;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  i3 = b->size[0];
  crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i3, A_nrows, A_ncols, x->size[1], false);
}

static void crs_prodAtx_kernel(const int32_T *row_ptr, const
  int32_T *col_ind, const real_T *val, const real_T *x, int x_m, real_T *b, int b_m, int nrows, int ncols, int nrhs,
  boolean_T ismt)
{
  int nthreads;
  int varargin_1;
  int boffset;
  int varargin_2;
  int istart;
  int iend;
  int xoffset;
  int b_nthreads;
  int i2;
  int thr_id;
  double y;
  double b_y;
  int offset;
  int i;
  int chunk;
  int b_remainder;
  int c_remainder;
  int r;

  if (ismt) {
    varargin_1 = omp_get_num_threads();
    varargin_2 = (int)rt_roundd(floor((double)b_m / (double)ncols));
    if (varargin_1 <= varargin_2) {
      nthreads = varargin_1;
    } else {
      nthreads = varargin_2;
    }

    varargin_1 = omp_get_thread_num();
    boffset = varargin_1 * ncols;
    get_local_chunk(nrows, nthreads, &istart, &iend);
  } else {
    nthreads = 1;
    boffset = 0;
    istart = 1;
    iend = nrows;
  }

  if (istart <= iend) {
    xoffset = -1;
    for (varargin_1 = 1; varargin_1 <= nrhs; varargin_1++) {
      i2 = boffset + ncols;
      for (varargin_2 = boffset; varargin_2 + 1 <= i2; varargin_2++) {
        b[varargin_2] = 0.0;
      }

      for (i = istart; i <= iend; i++) {
        i2 = row_ptr[i] - 1;
        for (varargin_2 = row_ptr[i - 1]; varargin_2 <= i2; varargin_2++)
        {
          r = (boffset + col_ind[varargin_2 - 1]) - 1;
          b[r] += x[i + xoffset] * val[varargin_2 - 1];
        }
      }

      xoffset += x_m;
      boffset += b_m;
    }
  }

  if (nthreads > 1) {

#pragma omp barrier

    b_nthreads = omp_get_num_threads();
    if (b_nthreads == 1) {
      istart = 0;
      iend = ncols;
    } else {
      thr_id = omp_get_thread_num();
      y = (double)ncols / (double)b_nthreads;
      if (y < 0.0) {
        b_y = ceil(y);
      } else {
        b_y = floor(y);
      }

      chunk = (int)rt_roundd(b_y);
      b_remainder = ncols - b_nthreads * chunk;
      if (b_remainder <= thr_id) {
        c_remainder = b_remainder;
      } else {
        c_remainder = thr_id;
      }

      istart = thr_id * chunk + c_remainder;
      iend = (istart + chunk) + (thr_id < b_remainder);
    }

    offset = ncols;
    for (varargin_2 = 2; varargin_2 <= nthreads; varargin_2++) {
      boffset = 1;
      for (varargin_1 = 1; varargin_1 <= nrhs; varargin_1++) {
        i2 = (boffset + iend) - 1;
        for (i = boffset + istart; i <= i2; i++) {
          b[i - 1] += b[(offset + i) - 1];
        }

        boffset += b_m;
      }

      offset += ncols;
    }
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

static void get_local_chunk(int m, int varargin_2, int *istart, int *iend)
{
  int thr_id;
  double y;
  double b_y;
  int chunk;
  int b_remainder;
  int c_remainder;
  if (varargin_2 == 1) {
    *istart = 1;
    *iend = m;
  } else {
    thr_id = omp_get_thread_num();
    y = (double)m / (double)varargin_2;
    if (y < 0.0) {
      b_y = ceil(y);
    } else {
      b_y = floor(y);
    }

    chunk = (int)rt_roundd(b_y);
    b_remainder = m - varargin_2 * chunk;
    if (b_remainder <= thr_id) {
      c_remainder = b_remainder;
    } else {
      c_remainder = thr_id;
    }

    *istart = (thr_id * chunk + c_remainder) + 1;
    *iend = ((*istart + chunk) + (thr_id < b_remainder)) - 1;
  }
}

static void m2c_error(void)
{
  M2C_error("crs_prodAtx:BufferTooSmal",
            "Buffer space for output b is too small.");
}

static void m2c_warn(void)
{
  M2C_warn("crs_prodAtx:NestedParallel",
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

void crs_prodAtx(const struct0_T *A, const emxArray_real_T *x, emxArray_real_T
                 *b, const emxArray_int32_T *nthreads)
{
  int n;
  boolean_T b0;
  int i1;
  if ((b->size[0] < A->ncols) || (b->size[1] < x->size[1])) {
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
    crs_prodAtx_kernel(A->row_ptr->data, A->col_ind->data, A->val->data, x->data, x->size[0], b->data, i1, A->nrows, A->ncols, x->size[1], n > 1);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    i1 = b->size[0];
    crs_prodAtx_kernel(A->row_ptr->data, A->col_ind->data, A->val->data, x->data, x->size[0], b->data, i1, A->nrows, A->ncols, x->size[1], b0);
  }
}

void crs_prodAtx_initialize(void)
{
}

void crs_prodAtx_ser(const struct0_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  b_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtx_ser1(const struct0_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b)
{
  c_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtx_terminate(void)
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
