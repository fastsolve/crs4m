#include "sqnormf.h"
#include "omp.h"
#include "m2c.h"

static double b_sqnormf(const emxArray_real_T *A);
static void m2c_warn(void);
static double rt_roundd(double u);
static double b_sqnormf(const emxArray_real_T *A)
{
  int iend;
  double s_local;
  int i;
  iend = A->size[0] * A->size[1];
  s_local = 0.0;
  for (i = 0; i + 1 <= iend; i++) {
    s_local += A->data[i] * A->data[i];
  }

  return s_local;
}

static void m2c_warn(void)
{
  M2C_warn("sqnormf:NestedParallel",
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

void emxInitArray_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxInit_int32_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}

void sqnormf(const emxArray_real_T *A, double *sqnrm, const emxArray_int32_T
             *nthreads)
{
  int iend;
  int n;
  boolean_T b0;
  int istart;
  double s_local;
  int i;
  int b_nthreads;
  int thr_id;
  double y;
  double b_y;
  int chunk;
  int b_remainder;
  int c_remainder;
  iend = A->size[0] * A->size[1];
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

    s_local = 0.0;
    i = 0;

#pragma omp parallel for default(none) shared(A, iend) private(i) reduction(+:s_local) num_threads(nthreads->data[0]) 
    for (i = 0; i <=(iend)-( 1 ); i++) {
      s_local += A->data[i] * A->data[i];
    }

    *sqnrm = s_local;
  } else {
    if (b0) {

#pragma omp barrier
#pragma omp single

      M2C_BEGIN_REGION(/*omp single*/)

      *sqnrm = 0.0;

      M2C_END_REGION(/*omp single*/)

      b_nthreads = omp_get_num_threads();
      if (b_nthreads == 1) {
        istart = 0;
      } else {
        thr_id = omp_get_thread_num();
        y = (double)iend / (double)b_nthreads;
        if (y < 0.0) {
          b_y = ceil(y);
        } else {
          b_y = floor(y);
        }

        chunk = (int)rt_roundd(b_y);
        b_remainder = iend - b_nthreads * chunk;
        if (b_remainder <= thr_id) {
          c_remainder = b_remainder;
        } else {
          c_remainder = thr_id;
        }

        istart = thr_id * chunk + c_remainder;
        iend = (istart + chunk) + (thr_id < b_remainder);
      }
    } else {
      *sqnrm = 0.0;
      istart = 0;
    }

    s_local = 0.0;
    for (i = istart; i + 1 <= iend; i++) {
      s_local += A->data[i] * A->data[i];
    }

    if (b0) {

#pragma omp atomic

      *sqnrm += s_local;

#pragma omp barrier

    } else {
      *sqnrm = s_local;
    }
  }
}

void sqnormf_initialize(void)
{
}

double sqnormf_ser(const emxArray_real_T *A)
{
  return b_sqnormf(A);
}

void sqnormf_terminate(void)
{
}
