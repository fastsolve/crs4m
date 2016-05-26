#include "prodAtx.h"
#include "omp.h"
#include "m2c.h"

static void b_prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b);
static void c_prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b);
static void m2c_error(void);
static void m2c_warn(void);
static void prodAtx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static double rt_roundd(double u);
static void b_prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b)
{
  int i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A->size[1];
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int)sizeof(double));
  prodAtx_internal(A, x, b);
}

static void c_prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b)
{
  if ((b->size[0] < A->size[1]) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  prodAtx_internal(A, x, b);
}

static void m2c_error(void)
{
  M2C_error("prodAtx:BufferTooSmal", "Buffer space for output b is too small.");
}

static void m2c_warn(void)
{
  M2C_warn("prodAtx:NestedParallel",
           "You are trying to use nested parallel regions. Solution may be incorrect.");
}

static void prodAtx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  int nthreads;
  int thr_id;
  int jstart;
  double y;
  int jend;
  double b_y;
  int j;
  int chunk;
  int b_remainder;
  int k;
  int c_remainder;
  double t;
  int i;
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    jstart = 0;
    jend = A->size[1];
  } else {
    thr_id = omp_get_thread_num();
    y = (double)A->size[1] / (double)nthreads;
    if (y < 0.0) {
      b_y = ceil(y);
    } else {
      b_y = floor(y);
    }

    chunk = (int)rt_roundd(b_y);
    b_remainder = (int)((double)A->size[1] - (double)(nthreads * chunk));
    if (b_remainder <= thr_id) {
      c_remainder = b_remainder;
    } else {
      c_remainder = thr_id;
    }

    jstart = thr_id * chunk + c_remainder;
    jend = (jstart + chunk) + (thr_id < b_remainder);
  }

  for (j = jstart; j + 1 <= jend; j++) {
    for (k = 0; k + 1 <= x->size[1]; k++) {
      t = 0.0;
      for (i = 0; i + 1 <= A->size[0]; i++) {
        t += A->data[i + A->size[0] * j] * x->data[i + x->size[0] * k];
      }

      b->data[j + b->size[0] * k] = t;
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

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}

void prodAtx(const emxArray_real_T *A, const emxArray_real_T *x, emxArray_real_T
             *b, int nthreads)
{
  int n;
  int b_n;
  if ((b->size[0] < A->size[1]) || (b->size[1] < x->size[1])) {
    m2c_error();
  }

  n = omp_get_num_threads();
  b_n = omp_get_nested();
  if ((!(b_n != 0)) && (n > 1) && (nthreads > 1)) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    m2c_warn();

    M2C_END_REGION(/*omp master*/)

  }

#pragma omp parallel default(shared) num_threads(nthreads)
  M2C_BEGIN_REGION(/*omp parallel*/)

  prodAtx_internal(A, x, b);

  M2C_END_REGION(/*omp parallel*/)

}

void prodAtx_initialize(void)
{
}

void prodAtx_ser(const emxArray_real_T *A, const emxArray_real_T *x,
                 emxArray_real_T *b)
{
  b_prodAtx(A, x, b);
}

void prodAtx_ser1(const emxArray_real_T *A, const emxArray_real_T *x,
                  emxArray_real_T *b)
{
  c_prodAtx(A, x, b);
}

void prodAtx_terminate(void)
{
}
