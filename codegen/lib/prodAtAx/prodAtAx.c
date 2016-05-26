#include "prodAtAx.h"
#include "omp.h"
#include "m2c.h"

static void b_m2c_error(void);
static void b_prodAtAx(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static void b_prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b);
static void b_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b);
static void c_m2c_error(void);
static void c_prodAtAx(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static void d_m2c_error(void);
static void get_local_chunk(double m, int *istart, int *iend);
static void m2c_error(void);
static void m2c_warn(void);
static void prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                    emxArray_real_T *b);
static void prodAtx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static void prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                   emxArray_real_T *b);
static void prodAx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static double rt_roundd(double u);
static void b_m2c_error(void)
{
  M2C_error("prodAtAx:IncorrectBuffer", "Buffer Ax has incorrect size.");
}

static void b_prodAtAx(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  emxArray_real_T *Ax;
  int i0;
  int loop_ub;
  emxInit_real_T(&Ax, 2);
  i0 = Ax->size[0] * Ax->size[1];
  Ax->size[0] = A->size[0];
  Ax->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)Ax, i0, (int)sizeof(double));
  loop_ub = A->size[0] * x->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    Ax->data[i0] = 0.0;
  }

  b_prodAx(A, x, Ax);
  i0 = b->size[0] * b->size[1];
  b->size[0] = A->size[1];
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int)sizeof(double));
  loop_ub = A->size[1] * x->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b->data[i0] = 0.0;
  }

  b_prodAtx(A, Ax, b);
  emxFree_real_T(&Ax);
}

static void b_prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b)
{
  if ((b->size[0] < A->size[1]) || (b->size[1] < x->size[1])) {
    d_m2c_error();
  }

  prodAtx_internal(A, x, b);
}

static void b_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
    c_m2c_error();
  }

  prodAx_internal(A, x, b);
}

static void c_m2c_error(void)
{
  M2C_error("prodAx:BufferTooSmal", "Buffer space for output b is too small.");
}

static void c_prodAtAx(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  emxArray_real_T *Ax;
  int i2;
  int loop_ub;
  if ((b->size[0] < A->size[1]) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    m2c_error();

    M2C_END_REGION(/*omp master*/)

  }

  emxInit_real_T(&Ax, 2);
  i2 = Ax->size[0] * Ax->size[1];
  Ax->size[0] = A->size[0];
  Ax->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)Ax, i2, (int)sizeof(double));
  loop_ub = A->size[0] * x->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    Ax->data[i2] = 0.0;
  }

  b_prodAx(A, x, Ax);
  b_prodAtx(A, Ax, b);
  emxFree_real_T(&Ax);
}

static void d_m2c_error(void)
{
  M2C_error("prodAtx:BufferTooSmal", "Buffer space for output b is too small.");
}

static void get_local_chunk(double m, int *istart, int *iend)
{
  int nthreads;
  int thr_id;
  double y;
  double b_y;
  int chunk;
  int b_remainder;
  int c_remainder;
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    *istart = 1;
    *iend = (int)rt_roundd(m);
  } else {
    thr_id = omp_get_thread_num();
    y = m / (double)nthreads;
    if (y < 0.0) {
      b_y = ceil(y);
    } else {
      b_y = floor(y);
    }

    chunk = (int)rt_roundd(b_y);
    b_remainder = (int)rt_roundd(m - (double)(nthreads * chunk));
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
  M2C_error("prodAtAx:IncorrectBuffer", "Buffer b has incorrect size.");
}

static void m2c_warn(void)
{
  M2C_warn("prodAtAx:NestedParallel",
           "You are trying to use nested parallel regions, but nested parallelism is not enabled.");
}

static void prodAtx(const emxArray_real_T *A, const emxArray_real_T *x,
                    emxArray_real_T *b)
{
  if ((b->size[0] < A->size[1]) || (b->size[1] < x->size[1])) {
    d_m2c_error();
  }

  omp_get_num_threads();
  prodAtx_internal(A, x, b);
}

static void prodAtx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  int j;
  int jend;
  int k;
  double t;
  int i;
  get_local_chunk(A->size[1], &j, &jend);
  while (j <= jend) {
    for (k = 0; k + 1 <= x->size[1]; k++) {
      t = 0.0;
      for (i = 0; i + 1 <= A->size[0]; i++) {
        t += A->data[i + A->size[0] * (j - 1)] * x->data[i + x->size[0] * k];
      }

      b->data[(j + b->size[0] * k) - 1] = t;
    }

    j++;
  }
}

static void prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                   emxArray_real_T *b)
{
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
    c_m2c_error();
  }

  omp_get_num_threads();
  prodAx_internal(A, x, b);
}

static void prodAx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  int i;
  int iend;
  int k;
  double t;
  int j;
  get_local_chunk(A->size[0], &i, &iend);
  while (i <= iend) {
    for (k = 0; k + 1 <= x->size[1]; k++) {
      t = 0.0;
      for (j = 0; j + 1 <= A->size[1]; j++) {
        t += A->data[(i + A->size[0] * j) - 1] * x->data[j + x->size[0] * k];
      }

      b->data[(i + b->size[0] * k) - 1] = t;
    }

    i++;
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

void prodAtAx(const emxArray_real_T *A, const emxArray_real_T *x,
              emxArray_real_T *b, emxArray_real_T *Ax, const emxArray_int32_T
              *nthreads)
{
  int n;
  boolean_T b0;
  if ((b->size[0] < A->size[1]) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    m2c_error();

    M2C_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A->size[0]) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    b_m2c_error();

    M2C_END_REGION(/*omp master*/)

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

    prodAx(A, x, Ax);

#pragma omp barrier

    prodAtx(A, Ax, b);

    M2C_END_REGION(/*omp parallel*/)

  } else if (b0) {
    prodAx(A, x, Ax);

#pragma omp barrier

    prodAtx(A, Ax, b);
  } else {
    b_prodAx(A, x, Ax);
    b_prodAtx(A, Ax, b);
  }
}

void prodAtAx_initialize(void)
{
}

void prodAtAx_ser(const emxArray_real_T *A, const emxArray_real_T *x,
                  emxArray_real_T *b)
{
  b_prodAtAx(A, x, b);
}

void prodAtAx_ser1(const emxArray_real_T *A, const emxArray_real_T *x,
                   emxArray_real_T *b)
{
  c_prodAtAx(A, x, b);
}

void prodAtAx_terminate(void)
{
}
