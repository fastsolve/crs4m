#include "prodAx.h"
#include "omp.h"
#include "m2c.h"

static void b_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b);
static void c_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b);
static void m2c_error(void);
static void m2c_warn(void);
static void prodAx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static double rt_roundd(double u);
static void b_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  int i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A->size[0];
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int)sizeof(double));
  prodAx_internal(A, x, b);
}

static void c_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  prodAx_internal(A, x, b);
}

static void m2c_error(void)
{
  M2C_error("prodAx:BufferTooSmal", "Buffer space for output b is too small.");
}

static void m2c_warn(void)
{
  M2C_warn("prodAx:NestedParallel",
           "You are trying to use nested parallel regions. Solution may be incorrect.");
}

static void prodAx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  int nthreads;
  int thr_id;
  int istart;
  double y;
  int iend;
  double b_y;
  int i;
  int chunk;
  int b_remainder;
  int k;
  int c_remainder;
  double t;
  int j;
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    istart = 0;
    iend = A->size[0];
  } else {
    thr_id = omp_get_thread_num();
    y = (double)A->size[0] / (double)nthreads;
    if (y < 0.0) {
      b_y = ceil(y);
    } else {
      b_y = floor(y);
    }

    chunk = (int)rt_roundd(b_y);
    b_remainder = (int)((double)A->size[0] - (double)(nthreads * chunk));
    if (b_remainder <= thr_id) {
      c_remainder = b_remainder;
    } else {
      c_remainder = thr_id;
    }

    istart = thr_id * chunk + c_remainder;
    iend = (istart + chunk) + (thr_id < b_remainder);
  }

  for (i = istart; i + 1 <= iend; i++) {
    for (k = 0; k + 1 <= x->size[1]; k++) {
      t = 0.0;
      for (j = 0; j + 1 <= A->size[1]; j++) {
        t += A->data[i + A->size[0] * j] * x->data[j + x->size[0] * k];
      }

      b->data[i + b->size[0] * k] = t;
    }
  }
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

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}

void prodAx(const emxArray_real_T *A, const emxArray_real_T *x, emxArray_real_T *
            b, const emxArray_int32_T *nthreads)
{
  int n;
  int b_n;
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  n = omp_get_num_threads();
  if (!(nthreads->size[0] == 0)) {
    b_n = omp_get_nested();
    if ((!(b_n != 0)) && (n > 1) && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      m2c_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    prodAx_internal(A, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    prodAx_internal(A, x, b);
  }
}

void prodAx_initialize(void)
{
}

void prodAx_ser(const emxArray_real_T *A, const emxArray_real_T *x,
                emxArray_real_T *b)
{
  b_prodAx(A, x, b);
}

void prodAx_ser1(const emxArray_real_T *A, const emxArray_real_T *x,
                 emxArray_real_T *b)
{
  c_prodAx(A, x, b);
}

void prodAx_terminate(void)
{
}
