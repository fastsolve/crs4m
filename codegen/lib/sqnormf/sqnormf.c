#include "sqnormf.h"
#include "mpi.h"
#include "omp.h"
#include "m2c.h"

static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
static void allreduce(real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm varargin_1);
static real_T b_sqnormf(const emxArray_real_T *A);
static void c_sqnormf(const emxArray_real_T *A, real_T *sqnrm, const
                      emxArray_int32_T *nthreads, MPI_Comm varargin_1);
static void msg_error(const emxArray_char_T *varargin_3);
static void msg_warn(void);
static real_T rt_roundd(real_T u);
static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce(sptr, rptr, count, datatype, op, comm);
}

static void allreduce(real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm varargin_1)
{
  int32_T size;
  void * ptr;
  MPI_Comm_size(varargin_1, &size);
  if (size > 1) {

    ptr = (void *)(b);
    MMPI_Allreduce(MPI_IN_PLACE, ptr, sz_in, MPI_DOUBLE, op, varargin_1);
  }
}

static real_T b_sqnormf(const emxArray_real_T *A)
{
  real_T sqnrm;
  int32_T iend;
  int32_T i;
  iend = A->size[0] * A->size[1];
  sqnrm = 0.0;
  for (i = 0; i + 1 <= iend; i++) {
    sqnrm += A->data[i] * A->data[i];
  }

  return sqnrm;
}

static void c_sqnormf(const emxArray_real_T *A, real_T *sqnrm, const
                      emxArray_int32_T *nthreads, MPI_Comm varargin_1)
{
  int32_T iend;
  int32_T istart;
  boolean_T b1;
  real_T s_local;
  int32_T thr_id;
  int32_T chunk;
  iend = A->size[0] * A->size[1];
  istart = omp_get_num_threads();
  b1 = (istart > 1);
  if (!(nthreads->size[0] == 0)) {
    istart = omp_get_nested();
    if ((!(istart != 0)) && b1 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

    s_local = 0.0;
    istart = 0;

#pragma omp parallel for default(none) shared(A, iend) private(istart) reduction(+:s_local) num_threads(nthreads->data[0]) 
    for (istart = 0; istart <=(iend)-( 1 ); istart++) {
      s_local += A->data[istart] * A->data[istart];
    }

    *sqnrm = s_local;
  } else {
    if (b1) {

#pragma omp barrier
#pragma omp single

      M2C_BEGIN_REGION(/*omp single*/)

      *sqnrm = 0.0;

      M2C_END_REGION(/*omp single*/)

      istart = omp_get_num_threads();
      if (istart == 1) {
        istart = 0;
      } else {
        thr_id = omp_get_thread_num();
        s_local = (real_T)iend / (real_T)istart;
        if (s_local < 0.0) {
          s_local = ceil(s_local);
        } else {
          s_local = floor(s_local);
        }

        chunk = (int32_T)rt_roundd(s_local);
        iend -= istart * chunk;
        if (iend <= thr_id) {
          istart = iend;
        } else {
          istart = thr_id;
        }

        istart += thr_id * chunk;
        iend = (istart + chunk) + (thr_id < iend);
      }
    } else {
      *sqnrm = 0.0;
      istart = 0;
    }

    s_local = 0.0;
    while (istart + 1 <= iend) {
      s_local += A->data[istart] * A->data[istart];
      istart++;
    }

    if (b1) {

#pragma omp atomic

      *sqnrm += s_local;

#pragma omp barrier

    } else {
      *sqnrm = s_local;
    }
  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  allreduce(sqnrm, 1, MPI_SUM, varargin_1);

  M2C_END_REGION(/*omp single*/)

}

static void msg_error(const emxArray_char_T *varargin_3)
{
  const char * msgid;
  emxArray_char_T *b_varargin_3;
  int32_T i0;
  int32_T loop_ub;
  msgid = "MPI_Comm:WrongType";
  emxInit_char_T(&b_varargin_3, 2);
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    i0 = b_varargin_3->size[0] * b_varargin_3->size[1];
    b_varargin_3->size[0] = 1;
    b_varargin_3->size[1] = varargin_3->size[1];
    emxEnsureCapacity((emxArray__common *)b_varargin_3, i0, (int32_T)sizeof
                      (char_T));
    loop_ub = varargin_3->size[0] * varargin_3->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_varargin_3->data[i0] = varargin_3->data[i0];
    }

    mexErrMsgIdAndTxt(msgid, "Incorrect data type %s. Expected MPI_Datatype.",
                      &b_varargin_3->data[0]);
  } else {
    i0 = b_varargin_3->size[0] * b_varargin_3->size[1];
    b_varargin_3->size[0] = 1;
    b_varargin_3->size[1] = varargin_3->size[1];
    emxEnsureCapacity((emxArray__common *)b_varargin_3, i0, (int32_T)sizeof
                      (char_T));
    loop_ub = varargin_3->size[0] * varargin_3->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_varargin_3->data[i0] = varargin_3->data[i0];
    }

    printf("Error %s\nIncorrect data type %s. Expected MPI_Datatype.", msgid,
           &b_varargin_3->data[0]);
  }

  emxFree_char_T(&b_varargin_3);
}

static void msg_warn(void)
{
  const char * msgid;
  msgid = "sqnormf:NestedParallel";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexWarnMsgIdAndTxt(msgid,
                       "You are trying to use nested parallel regions, but nested parallelism is not enabled.");
  } else {
    printf("Warning %s\nYou are trying to use nested parallel regions, but nested parallelism is not enabled.",
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

void sqnormf(const emxArray_real_T *A, real_T *sqnrm, const emxArray_int32_T
             *nthreads)
{
  int32_T iend;
  int32_T istart;
  boolean_T b0;
  real_T s_local;
  int32_T thr_id;
  int32_T chunk;
  iend = A->size[0] * A->size[1];
  istart = omp_get_num_threads();
  b0 = (istart > 1);
  if (!(nthreads->size[0] == 0)) {
    istart = omp_get_nested();
    if ((!(istart != 0)) && b0 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

    s_local = 0.0;
    istart = 0;

#pragma omp parallel for default(none) shared(A, iend) private(istart) reduction(+:s_local) num_threads(nthreads->data[0]) 
    for (istart = 0; istart <=(iend)-( 1 ); istart++) {
      s_local += A->data[istart] * A->data[istart];
    }

    *sqnrm = s_local;
  } else {
    if (b0) {

#pragma omp barrier
#pragma omp single

      M2C_BEGIN_REGION(/*omp single*/)

      *sqnrm = 0.0;

      M2C_END_REGION(/*omp single*/)

      istart = omp_get_num_threads();
      if (istart == 1) {
        istart = 0;
      } else {
        thr_id = omp_get_thread_num();
        s_local = (real_T)iend / (real_T)istart;
        if (s_local < 0.0) {
          s_local = ceil(s_local);
        } else {
          s_local = floor(s_local);
        }

        chunk = (int32_T)rt_roundd(s_local);
        iend -= istart * chunk;
        if (iend <= thr_id) {
          istart = iend;
        } else {
          istart = thr_id;
        }

        istart += thr_id * chunk;
        iend = (istart + chunk) + (thr_id < iend);
      }
    } else {
      *sqnrm = 0.0;
      istart = 0;
    }

    s_local = 0.0;
    while (istart + 1 <= iend) {
      s_local += A->data[istart] * A->data[istart];
      istart++;
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

void sqnormf_mpi(const emxArray_real_T *A, real_T *s, const emxArray_int32_T
                 *nthreads, const struct_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i1;
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
      i1 = comm->type->size[k];
      if (i1 != 7 * k + 1) {
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
    i1 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i1, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i1 = 0; i1 < k; i1++) {
      b_comm->data[b_comm->size[0] * i1] = comm->type->data[comm->type->size[0] *
        i1];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i1 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i1, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i1 = 0; i1 < k; i1++) {
    data->data[i1] = comm->data->data[i1];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  c_sqnormf(A, s, nthreads, c_comm);
  emxFree_uint8_T(&data);
}

real_T sqnormf_ser(const emxArray_real_T *A)
{
  return b_sqnormf(A);
}

void sqnormf_terminate(void)
{
}
