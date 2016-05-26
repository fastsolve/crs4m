#include "dotprod.h"
#include "omp.h"
#include "m2c.h"

static void accu_partsum(emxArray_real_T *prod);
static void b_dotprod(const emxArray_real_T *u, const emxArray_real_T *v,
                      emxArray_real_T *prod);
static void b_m2c_error(void);
static void dotprod_partial(const emxArray_real_T *u, const emxArray_real_T *v,
  emxArray_real_T *prod);
static void m2c_error(void);
static void m2c_warn(void);
static double rt_roundd(double u);
static void accu_partsum(emxArray_real_T *prod)
{
  int nthreads;
  int thr_id;
  int istart;
  double y;
  int iend;
  double b_y;
  int varargin_1;
  int chunk;
  int n;
  int b_remainder;
  int k;
  int c_remainder;
  int j;
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    istart = 0;
    iend = prod->size[1];
  } else {
    thr_id = omp_get_thread_num();
    y = (double)prod->size[1] / (double)nthreads;
    if (y < 0.0) {
      b_y = ceil(y);
    } else {
      b_y = floor(y);
    }

    chunk = (int)rt_roundd(b_y);
    b_remainder = prod->size[1] - nthreads * chunk;
    if (b_remainder <= thr_id) {
      c_remainder = b_remainder;
    } else {
      c_remainder = thr_id;
    }

    istart = thr_id * chunk + c_remainder;
    iend = (istart + chunk) + (thr_id < b_remainder);
  }

  varargin_1 = omp_get_num_threads();
  if (varargin_1 <= prod->size[0]) {
    n = varargin_1;
  } else {
    n = prod->size[0];
  }

  for (k = istart; k + 1 <= iend; k++) {
    for (j = 2; j <= n; j++) {
      prod->data[prod->size[0] * k] += prod->data[(j + prod->size[0] * k) - 1];
    }
  }
}

static void b_dotprod(const emxArray_real_T *u, const emxArray_real_T *v,
                      emxArray_real_T *prod)
{
  int i0;
  int k;
  int i;
  if ((u->size[0] != v->size[0]) || (u->size[1] != v->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    m2c_error();

    M2C_END_REGION(/*omp master*/)

  }

  i0 = prod->size[0] * prod->size[1];
  prod->size[0] = 1;
  prod->size[1] = u->size[1];
  emxEnsureCapacity((emxArray__common *)prod, i0, (int)sizeof(double));
  for (k = 0; k + 1 <= u->size[1]; k++) {
    prod->data[prod->size[0] * k] = 0.0;
    for (i = 0; i + 1 <= u->size[0]; i++) {
      prod->data[prod->size[0] * k] += u->data[i + u->size[0] * k] * v->data[i +
        v->size[0] * k];
    }
  }
}

static void b_m2c_error(void)
{
  M2C_error("dotprod:IncorrectSize",
            "prod and u must have the same number of columns.");
}

static void dotprod_partial(const emxArray_real_T *u, const emxArray_real_T *v,
  emxArray_real_T *prod)
{
  int n;
  int nthreads;
  int thr_id;
  int istart;
  double y;
  int iend;
  double b_y;
  int k;
  int chunk;
  int b_remainder;
  int i;
  int c_remainder;
  n = omp_get_thread_num();
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    istart = 0;
    iend = u->size[0];
  } else {
    thr_id = omp_get_thread_num();
    y = (double)u->size[0] / (double)nthreads;
    if (y < 0.0) {
      b_y = ceil(y);
    } else {
      b_y = floor(y);
    }

    chunk = (int)rt_roundd(b_y);
    b_remainder = u->size[0] - nthreads * chunk;
    if (b_remainder <= thr_id) {
      c_remainder = b_remainder;
    } else {
      c_remainder = thr_id;
    }

    istart = thr_id * chunk + c_remainder;
    iend = (istart + chunk) + (thr_id < b_remainder);
  }

  for (k = 0; k + 1 <= u->size[1]; k++) {
    prod->data[n + prod->size[0] * k] = 0.0;
    for (i = istart; i + 1 <= iend; i++) {
      prod->data[n + prod->size[0] * k] += u->data[i + u->size[0] * k] * v->
        data[i + v->size[0] * k];
    }
  }
}

static void m2c_error(void)
{
  M2C_error("dotprod:IncorrectSize", "Dimensions for u and v must match.");
}

static void m2c_warn(void)
{
  M2C_warn("dotprod:NestedParallel",
           "You are trying to use nested parallel regions, but nested parallelism is not enabled.");
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

void dotprod(const emxArray_real_T *u, const emxArray_real_T *v, emxArray_real_T
             *prod, const emxArray_int32_T *nthreads)
{
  int n;
  int k;
  int i;
  if ((u->size[0] != v->size[0]) || (u->size[1] != v->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    m2c_error();

    M2C_END_REGION(/*omp master*/)

  }

  if ((prod->size[0] < 1) || (prod->size[1] != u->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    b_m2c_error();

    M2C_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if (!(n != 0)) {
      n = omp_get_num_threads();
      if ((n > 1) && (nthreads->data[0] > 1)) {

#pragma omp master

        M2C_BEGIN_REGION(/*omp master*/)

        m2c_warn();

        M2C_END_REGION(/*omp master*/)

      }
    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    dotprod_partial(u, v, prod);

#pragma omp barrier

    accu_partsum(prod);

    M2C_END_REGION(/*omp parallel*/)

  } else if (n > 1) {
    dotprod_partial(u, v, prod);

#pragma omp barrier

    accu_partsum(prod);

#pragma omp barrier

  } else {
    for (k = 0; k + 1 <= u->size[1]; k++) {
      prod->data[prod->size[0] * k] = 0.0;
      for (i = 0; i + 1 <= u->size[0]; i++) {
        prod->data[prod->size[0] * k] += u->data[i + u->size[0] * k] * v->data[i
          + v->size[0] * k];
      }
    }
  }
}

void dotprod_initialize(void)
{
}

void dotprod_ser(const emxArray_real_T *u, const emxArray_real_T *v,
                 emxArray_real_T *prod)
{
  b_dotprod(u, v, prod);
}

void dotprod_terminate(void)
{
}

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}
