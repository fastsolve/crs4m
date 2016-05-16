#include "dotprod.h"
#include "omp.h"
#include "m2c.h"

static void accu_partsum(emxArray_real_T *prod);
static void b_dotprod(const emxArray_real_T *u, const emxArray_real_T *v,
                      emxArray_real_T *prod);
static void b_m2c_error(void);
static void c_dotprod(const emxArray_real_T *u, const emxArray_real_T *v,
                      emxArray_real_T *prod, const emxArray_int32_T *nthreads,
                      MPI_Comm varargin_1);
static void c_m2c_error(const emxArray_char_T *varargin_3);
static void dotprod_partial(const emxArray_real_T *u, const emxArray_real_T *v,
  emxArray_real_T *prod);
static void emxFreeStruct_struct0_T(struct0_T *pStruct);
static void emxInitStruct_struct0_T(struct0_T *pStruct);
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

    PLC_BEGIN_REGION(/*omp master*/)

    m2c_error();

    PLC_END_REGION(/*omp master*/)

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

static void c_dotprod(const emxArray_real_T *u, const emxArray_real_T *v,
                      emxArray_real_T *prod, const emxArray_int32_T *nthreads,
                      MPI_Comm varargin_1)
{
  int n;
  int k;
  int i;
  int size;
  void * rptr;
  if ((u->size[0] != v->size[0]) || (u->size[1] != v->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    m2c_error();

    PLC_END_REGION(/*omp master*/)

  }

  if ((prod->size[0] < 1) || (prod->size[1] != u->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    b_m2c_error();

    PLC_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if (!(n != 0)) {
      n = omp_get_num_threads();
      if ((n > 1) && (nthreads->data[0] > 1)) {

#pragma omp master

        PLC_BEGIN_REGION(/*omp master*/)

        m2c_warn();

        PLC_END_REGION(/*omp master*/)

      }
    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    PLC_BEGIN_REGION(/*omp parallel*/)

    dotprod_partial(u, v, prod);

#pragma omp barrier

    accu_partsum(prod);

    PLC_END_REGION(/*omp parallel*/)

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

#pragma omp single

  PLC_BEGIN_REGION(/*omp single*/)

  n = prod->size[1];
  MPI_Comm_size(varargin_1, &size);
  if (size > 1) {

    rptr = (void *)(&prod->data[0]);
    MPI_Allreduce(MPI_IN_PLACE, rptr, n, MPI_DOUBLE, MPI_SUM, varargin_1);
  }

  PLC_END_REGION(/*omp single*/)

}

static void c_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int i1;
  int loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i1 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i1, (int)sizeof(char));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_varargin_3->data[i1] = varargin_3->data[i1];
  }

  M2C_error("MPI_Comm:WrongType", "Incorrect data type %s. Expected MPI_Comm.",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
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

static void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_uint8_T(&pStruct->data);
  emxFree_char_T(&pStruct->type);
}

static void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_uint8_T(&pStruct->data, 1);
  emxInit_char_T(&pStruct->type, 2);
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

    PLC_BEGIN_REGION(/*omp master*/)

    m2c_error();

    PLC_END_REGION(/*omp master*/)

  }

  if ((prod->size[0] < 1) || (prod->size[1] != u->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    b_m2c_error();

    PLC_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if (!(n != 0)) {
      n = omp_get_num_threads();
      if ((n > 1) && (nthreads->data[0] > 1)) {

#pragma omp master

        PLC_BEGIN_REGION(/*omp master*/)

        m2c_warn();

        PLC_END_REGION(/*omp master*/)

      }
    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    PLC_BEGIN_REGION(/*omp parallel*/)

    dotprod_partial(u, v, prod);

#pragma omp barrier

    accu_partsum(prod);

    PLC_END_REGION(/*omp parallel*/)

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

void dotprod_mpi(const emxArray_real_T *u, const emxArray_real_T *v,
                 emxArray_real_T *prod, const emxArray_int32_T *n, const
                 struct0_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int k;
  int exitg2;
  int i3;
  boolean_T exitg1;
  emxArray_char_T *b_comm;
  static const char cv0[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  emxArray_uint8_T *data;
  MPI_Comm c_comm;
  p = false;
  b_p = false;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i3 = comm->type->size[k];
      if (i3 != 7 * k + 1) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      b_p = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(comm->type->size[1] == 0))) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 8)) {
      if (!(comm->type->data[k] == cv0[k])) {
        b_p = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (!b_p) {
  } else {
    p = true;
  }

  if (!p) {
    emxInit_char_T(&b_comm, 2);
    i3 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i3, (int)sizeof(char));
    k = comm->type->size[1];
    for (i3 = 0; i3 < k; i3++) {
      b_comm->data[b_comm->size[0] * i3] = comm->type->data[comm->type->size[0] *
        i3];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    c_m2c_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 1);
  i3 = data->size[0];
  data->size[0] = comm->data->size[0];
  emxEnsureCapacity((emxArray__common *)data, i3, (int)sizeof(unsigned char));
  k = comm->data->size[0];
  for (i3 = 0; i3 < k; i3++) {
    data->data[i3] = comm->data->data[i3];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  c_dotprod(u, v, prod, n, c_comm);
  emxFree_uint8_T(&data);
}

void dotprod_ser(const emxArray_real_T *u, const emxArray_real_T *v,
                 emxArray_real_T *prod)
{
  b_dotprod(u, v, prod);
}

void dotprod_terminate(void)
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
  emxInit_real_T(pEmxArray, numDimensions);
}

void emxInit_struct0_T(struct0_T *pStruct)
{
  emxInitStruct_struct0_T(pStruct);
}
