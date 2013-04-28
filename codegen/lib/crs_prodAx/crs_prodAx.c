#include "crs_prodAx.h"
#include "spalab_kernel.h"
#include "mpi.h"
#include "omp.h"
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

static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
static void allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1);
static void b_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2);
static void b_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b);
static void b_msg_error(const emxArray_char_T *varargin_3);
static void c_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3);
static void c_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b);
static void c_msg_error(void);
static void crs_prodAx_internal(int32_T nrows, const emxArray_int32_T *row_ptr,
  const emxArray_int32_T *col_ind, const emxArray_real_T *val, const
  emxArray_real_T *x, emxArray_real_T *b);
static void d_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b, const emxArray_int32_T *nthreads,
  MPI_Comm varargin_1);
static void e_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b, const emxArray_int32_T *nthreads,
  MPI_Comm varargin_1, const emxArray_real_T *varargin_2);
extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize);
extern void emxFree_char_T(emxArray_char_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxFree_uint8_T(emxArray_uint8_T **pEmxArray);
extern void emxInit_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
extern void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int32_T numDimensions);
static void f_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b, const emxArray_int32_T *nthreads,
  MPI_Comm varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3);
static void msg_error(void);
static void msg_warn(void);
static real_T rt_roundd(real_T u);
static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce(sptr, rptr, count, datatype, op, comm);
}

static void allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1)
{
  int32_T flag;
  void * ptr;
  MPI_Initialized(&flag);
  if (flag != 0) {
    MPI_Comm_size(varargin_1, &flag);
    if (flag > 1) {

      ptr = (void *)(&b->data[0]);
      MMPI_Allreduce(MPI_IN_PLACE, ptr, sz_in, MPI_DOUBLE, op, varargin_1);
    }
  }
}

static void b_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2)
{
  int32_T s1;
  int32_T i;
  void * ptr;
  s1 = sz_in + varargin_2->size[0];
  if (s1 > b->size[0] * b->size[1]) {
    c_msg_error();
  } else {
    for (i = sz_in; i + 1 <= s1; i++) {
      b->data[i] = varargin_2->data[i - sz_in];
    }
  }

  MPI_Initialized(&i);
  if (i != 0) {
    MPI_Comm_size(varargin_1, &i);
    if (i > 1) {

      ptr = (void *)(&b->data[0]);
      MMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
    }
  }
}

static void b_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A_nrows;
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
  crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
}

static void b_msg_error(const emxArray_char_T *varargin_3)
{
  const char * msgid;
  emxArray_char_T *b_varargin_3;
  int32_T i1;
  int32_T loop_ub;
  msgid = "MPI_Comm:WrongType";
  emxInit_char_T(&b_varargin_3, 2);
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    i1 = b_varargin_3->size[0] * b_varargin_3->size[1];
    b_varargin_3->size[0] = 1;
    b_varargin_3->size[1] = varargin_3->size[1];
    emxEnsureCapacity((emxArray__common *)b_varargin_3, i1, (int32_T)sizeof
                      (char_T));
    loop_ub = varargin_3->size[0] * varargin_3->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_varargin_3->data[i1] = varargin_3->data[i1];
    }

    mexErrMsgIdAndTxt(msgid, "Incorrect data type %s. Expected MPI_Datatype.",
                      &b_varargin_3->data[0]);
  } else {
    i1 = b_varargin_3->size[0] * b_varargin_3->size[1];
    b_varargin_3->size[0] = 1;
    b_varargin_3->size[1] = varargin_3->size[1];
    emxEnsureCapacity((emxArray__common *)b_varargin_3, i1, (int32_T)sizeof
                      (char_T));
    loop_ub = varargin_3->size[0] * varargin_3->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_varargin_3->data[i1] = varargin_3->data[i1];
    }

    printf("Error %s\nIncorrect data type %s. Expected MPI_Datatype.", msgid,
           &b_varargin_3->data[0]);
  }

  emxFree_char_T(&b_varargin_3);
}

static void c_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T b_varargin_2;
  int32_T c_varargin_2;
  int32_T s1;
  void * ptr;
  b_varargin_2 = varargin_2->size[0];
  if (varargin_3 <= b_varargin_2) {
    b_varargin_2 = varargin_3;
  }

  if (b_varargin_2 < 0) {
    c_varargin_2 = 0;
  } else {
    c_varargin_2 = b_varargin_2;
  }

  s1 = sz_in + c_varargin_2;
  if (s1 > b->size[0] * b->size[1]) {
    c_msg_error();
  } else {
    for (b_varargin_2 = sz_in; b_varargin_2 + 1 <= s1; b_varargin_2++) {
      b->data[b_varargin_2] = varargin_2->data[b_varargin_2 - sz_in];
    }
  }

  MPI_Initialized(&b_varargin_2);
  if (b_varargin_2 != 0) {
    MPI_Comm_size(varargin_1, &b_varargin_2);
    if (b_varargin_2 > 1) {

      ptr = (void *)(&b->data[0]);
      MMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
    }
  }
}

static void c_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b)
{
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
}

static void c_msg_error(void)
{
  const char * msgid;
  msgid = "allreduce:bufferOverflow";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid,
                      "The given buffer is too small for the piggy-back message.");
  } else {
    printf("Error %s\nThe given buffer is too small for the piggy-back message.",
           msgid);
  }
}

static void crs_prodAx_internal(int32_T nrows, const emxArray_int32_T *row_ptr,
  const emxArray_int32_T *col_ind, const emxArray_real_T *val, const
  emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T nthreads;
  int32_T thr_id;
  int32_T i2;
  int32_T k;
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    thr_id = 0;
    nthreads = nrows;
  } else {
    thr_id = omp_get_thread_num();
    if (nthreads >= nrows) {
      nthreads = thr_id + (thr_id < nrows);
    } else {
      nthreads = (int32_T)rt_roundd(ceil((real_T)nrows / (real_T)nthreads));
      thr_id *= nthreads;
      nthreads += thr_id;
      if (nrows <= nthreads) {
        nthreads = nrows;
      }
    }
  }

  while (thr_id + 1 <= nthreads) {
    i2 = x->size[1];
    for (k = 0; k + 1 <= i2; k++) {
      b->data[thr_id + b->size[0] * k] = SPL_ddot(&val->data[row_ptr->
        data[thr_id] - 1], &col_ind->data[row_ptr->data[thr_id] - 1], &x->data
        [x->size[0] * k], row_ptr->data[thr_id + 1] - row_ptr->data[thr_id]);
    }

    thr_id++;
  }
}

static void d_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b, const emxArray_int32_T *nthreads,
  MPI_Comm varargin_1)
{
  int32_T n;
  boolean_T b0;
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b0 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b0 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
  }

  if ((nthreads->size[0] == 0) && b0) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  allreduce(b, (int32_T)rt_roundd((real_T)A_nrows * (real_T)x->size[1]), MPI_SUM,
            varargin_1);

            M2C_END_REGION(/*omp single*/)

}

static void e_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b, const emxArray_int32_T *nthreads,
  MPI_Comm varargin_1, const emxArray_real_T *varargin_2)
{
  int32_T n;
  boolean_T b1;
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b1 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b1 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
  }

  if ((nthreads->size[0] == 0) && b1) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  b_allreduce(b, (int32_T)rt_roundd((real_T)A_nrows * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2);

              M2C_END_REGION(/*omp single*/)

}










static void f_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b, const emxArray_int32_T *nthreads,
  MPI_Comm varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b2;
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b2 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b2 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
  }

  if ((nthreads->size[0] == 0) && b2) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  c_allreduce(b, (int32_T)rt_roundd((real_T)A_nrows * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2, varargin_3);

              M2C_END_REGION(/*omp single*/)

}

static void msg_error(void)
{
  const char * msgid;
  msgid = "crs_prodAx:BufferTooSmal";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer space for output b is too small.");
  } else {
    printf("Error %s\nBuffer space for output b is too small.", msgid);
  }
}

static void msg_warn(void)
{
  const char * msgid;
  msgid = "crs_prodAx:NestedParallel";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexWarnMsgIdAndTxt(msgid,
                       "You are trying to use nested parallel regions. Solution may be incorrect.");
  } else {
    printf("Warning %s\nYou are trying to use nested parallel regions. Solution may be incorrect.",
           msgid);
  }
}

static real_T rt_roundd(real_T u)
{
  real_T y;
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

void crs_prodAx(const struct_T *A, const emxArray_real_T *x, emxArray_real_T *b,
                const emxArray_int32_T *nthreads)
{
  int32_T n;
  int32_T b_n;
  if ((b->size[0] < A->nrows) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  if (!(nthreads->size[0] == 0)) {
    b_n = omp_get_nested();
    if ((!(b_n != 0)) && (n > 1) && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx_internal(A->nrows, A->row_ptr, A->col_ind, A->val, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    crs_prodAx_internal(A->nrows, A->row_ptr, A->col_ind, A->val, x, b);
  }
}

void crs_prodAx_initialize(void)
{
}

void crs_prodAx_mpi(const struct_T *A, const emxArray_real_T *x, emxArray_real_T
                    *b, const emxArray_int32_T *nthreads, const b_struct_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i3;
  boolean_T exitg1;
  static const char_T cv0[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  emxArray_char_T *b_comm;
  emxArray_uint8_T *data;
  MPI_Comm c_comm;
  p = FALSE;
  b_p = FALSE;
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
      b_p = TRUE;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(comm->type->size[1] == 0))) {
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (k <= comm->type->size[1] - 1)) {
      if (!(comm->type->data[k] == cv0[k])) {
        b_p = FALSE;
        exitg1 = TRUE;
      } else {
        k++;
      }
    }
  }

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (!p) {
    emxInit_char_T(&b_comm, 2);
    i3 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i3, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i3 = 0; i3 < k; i3++) {
      b_comm->data[b_comm->size[0] * i3] = comm->type->data[comm->type->size[0] *
        i3];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    b_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i3 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i3, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i3 = 0; i3 < k; i3++) {
    data->data[i3] = comm->data->data[i3];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  d_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b, nthreads, c_comm);
  emxFree_uint8_T(&data);
}

void crs_prodAx_mpip(const struct_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads, const
                     b_struct_T *comm, const emxArray_real_T *pbmsg)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i4;
  boolean_T exitg1;
  static const char_T cv1[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  emxArray_char_T *b_comm;
  emxArray_uint8_T *data;
  MPI_Comm c_comm;
  p = FALSE;
  b_p = FALSE;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i4 = comm->type->size[k];
      if (i4 != 7 * k + 1) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      b_p = TRUE;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(comm->type->size[1] == 0))) {
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (k <= comm->type->size[1] - 1)) {
      if (!(comm->type->data[k] == cv1[k])) {
        b_p = FALSE;
        exitg1 = TRUE;
      } else {
        k++;
      }
    }
  }

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (!p) {
    emxInit_char_T(&b_comm, 2);
    i4 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i4, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i4 = 0; i4 < k; i4++) {
      b_comm->data[b_comm->size[0] * i4] = comm->type->data[comm->type->size[0] *
        i4];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    b_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i4 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i4, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i4 = 0; i4 < k; i4++) {
    data->data[i4] = comm->data->data[i4];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  e_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b, nthreads, c_comm,
               pbmsg);
  emxFree_uint8_T(&data);
}

void crs_prodAx_mpip1(const struct_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b, const emxArray_int32_T *nthreads,
                      const b_struct_T *comm, const emxArray_real_T *pbmsg,
                      int32_T pbsz)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i5;
  boolean_T exitg1;
  static const char_T cv2[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  emxArray_char_T *b_comm;
  emxArray_uint8_T *data;
  MPI_Comm c_comm;
  p = FALSE;
  b_p = FALSE;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i5 = comm->type->size[k];
      if (i5 != 7 * k + 1) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      b_p = TRUE;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(comm->type->size[1] == 0))) {
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (k <= comm->type->size[1] - 1)) {
      if (!(comm->type->data[k] == cv2[k])) {
        b_p = FALSE;
        exitg1 = TRUE;
      } else {
        k++;
      }
    }
  }

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (!p) {
    emxInit_char_T(&b_comm, 2);
    i5 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i5, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i5 = 0; i5 < k; i5++) {
      b_comm->data[b_comm->size[0] * i5] = comm->type->data[comm->type->size[0] *
        i5];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    b_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i5 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i5, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i5 = 0; i5 < k; i5++) {
    data->data[i5] = comm->data->data[i5];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  f_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b, nthreads, c_comm,
               pbmsg, pbsz);
  emxFree_uint8_T(&data);
}

void crs_prodAx_ser(const struct_T *A, const emxArray_real_T *x, emxArray_real_T
                    *b)
{
  b_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b);
}

void crs_prodAx_ser1(const struct_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  c_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, b);
}

void crs_prodAx_terminate(void)
{
}









emxArray_char_T *emxCreateWrapper_char_T(char_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_char_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_char_T(&emx, 2);
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

emxArray_uint8_T *emxCreateWrapper_uint8_T(uint8_T *data, int32_T rows, int32_T
  cols)
{
  emxArray_uint8_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  emxInit_uint8_T(&emx, 2);
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








