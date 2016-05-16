#include "prodAx.h"
#include "mpi.h"
#include "omp.h"
#include "m2c.h"

static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
static void allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1);
static void b_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2);
static void b_msg_error(const emxArray_char_T *varargin_3);
static void b_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b);
static void c_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3);
static void c_msg_error(void);
static void c_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b);
static void d_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads,
                     MPI_Comm varargin_1);
static void e_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads,
                     MPI_Comm varargin_1, const emxArray_real_T *varargin_2);
static void f_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads,
                     MPI_Comm varargin_1, const emxArray_real_T *varargin_2,
                     int32_T varargin_3);
static void msg_error(void);
static void msg_warn(void);
static void prodAx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b);
static real_T rt_roundd(real_T u);
static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce(sptr, rptr, count, datatype, op, comm);
}

static void allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1)
{
  int32_T size;
  void * ptr;
  MPI_Comm_size(varargin_1, &size);
  if (size > 1) {

    ptr = (void *)(&b->data[0]);
    MMPI_Allreduce(MPI_IN_PLACE, ptr, sz_in, MPI_DOUBLE, op, varargin_1);
  }
}

static void b_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2)
{
  int32_T s1;
  int32_T i;
  void * ptr;
  s1 = sz_in + varargin_2->size[0] * varargin_2->size[1];
  if (s1 > b->size[0] * b->size[1]) {
    c_msg_error();
  } else {
    for (i = sz_in; i + 1 <= s1; i++) {
      b->data[i] = varargin_2->data[i - sz_in];
    }
  }

  MPI_Comm_size(varargin_1, &i);
  if (i > 1) {

    ptr = (void *)(&b->data[0]);
    MMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
  }
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

static void b_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  int32_T i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A->size[0];
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
  prodAx_internal(A, x, b);
}

static void c_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T b_varargin_2;
  int32_T c_varargin_2;
  int32_T s1;
  void * ptr;
  b_varargin_2 = varargin_2->size[0] * varargin_2->size[1];
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

  MPI_Comm_size(varargin_1, &b_varargin_2);
  if (b_varargin_2 > 1) {

    ptr = (void *)(&b->data[0]);
    MMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
  }
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

static void c_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  prodAx_internal(A, x, b);
}

static void d_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads,
                     MPI_Comm varargin_1)
{
  int32_T n;
  boolean_T b0;
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
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

    prodAx_internal(A, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    prodAx_internal(A, x, b);
  }

  if ((nthreads->size[0] == 0) && b0) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  allreduce(b, (int32_T)rt_roundd((real_T)A->size[0] * (real_T)x->size[1]),
            MPI_SUM, varargin_1);

  M2C_END_REGION(/*omp single*/)

}

static void e_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads,
                     MPI_Comm varargin_1, const emxArray_real_T *varargin_2)
{
  int32_T n;
  boolean_T b1;
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
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

    prodAx_internal(A, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    prodAx_internal(A, x, b);
  }

  if ((nthreads->size[0] == 0) && b1) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  b_allreduce(b, (int32_T)rt_roundd((real_T)A->size[0] * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2);

  M2C_END_REGION(/*omp single*/)

}

static void f_prodAx(const emxArray_real_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads,
                     MPI_Comm varargin_1, const emxArray_real_T *varargin_2,
                     int32_T varargin_3)
{
  int32_T n;
  boolean_T b2;
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
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

    prodAx_internal(A, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    prodAx_internal(A, x, b);
  }

  if ((nthreads->size[0] == 0) && b2) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  c_allreduce(b, (int32_T)rt_roundd((real_T)A->size[0] * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2, varargin_3);

  M2C_END_REGION(/*omp single*/)

}

static void msg_error(void)
{
  const char * msgid;
  msgid = "prodAx:BufferTooSmal";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer space for output b is too small.");
  } else {
    printf("Error %s\nBuffer space for output b is too small.", msgid);
  }
}

static void msg_warn(void)
{
  const char * msgid;
  msgid = "prodAx:NestedParallel";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexWarnMsgIdAndTxt(msgid,
                       "You are trying to use nested parallel regions. Solution may be incorrect.");
  } else {
    printf("Warning %s\nYou are trying to use nested parallel regions. Solution may be incorrect.",
           msgid);
  }
}

static void prodAx_internal(const emxArray_real_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  int32_T ncols;
  int32_T nrhs;
  int32_T nthreads;
  int32_T istart;
  int32_T thr_id;
  real_T t;
  int32_T chunk;
  int32_T b_remainder;
  ncols = A->size[1];
  nrhs = x->size[1];
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    istart = 0;
    nthreads = A->size[0];
  } else {
    thr_id = omp_get_thread_num();
    t = (real_T)A->size[0] / (real_T)nthreads;
    if (t < 0.0) {
      t = ceil(t);
    } else {
      t = floor(t);
    }

    chunk = (int32_T)rt_roundd(t);
    b_remainder = (int32_T)rt_roundd((real_T)A->size[0] - (real_T)(nthreads *
      chunk));
    if (b_remainder <= thr_id) {
      nthreads = b_remainder;
    } else {
      nthreads = thr_id;
    }

    istart = thr_id * chunk + nthreads;
    nthreads = (istart + chunk) + (thr_id < b_remainder);
  }

  while (istart + 1 <= nthreads) {
    for (b_remainder = 0; b_remainder + 1 <= nrhs; b_remainder++) {
      t = 0.0;
      for (thr_id = 0; thr_id + 1 <= ncols; thr_id++) {
        t += A->data[istart + A->size[0] * thr_id] * x->data[thr_id + x->size[0]
          * b_remainder];
      }

      b->data[istart + b->size[0] * b_remainder] = t;
    }

    istart++;
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

void prodAx(const emxArray_real_T *A, const emxArray_real_T *x, emxArray_real_T *
            b, const emxArray_int32_T *nthreads)
{
  int32_T n;
  int32_T b_n;
  if ((b->size[0] < A->size[0]) || (b->size[1] < x->size[1])) {
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

    prodAx_internal(A, x, b);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    prodAx_internal(A, x, b);
  }
}

void prodAx_initialize(void)
{
}

void prodAx_mpi(const emxArray_real_T *A, const emxArray_real_T *x,
                emxArray_real_T *b, const emxArray_int32_T *nthreads, const
                struct_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i2;
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
      i2 = comm->type->size[k];
      if (i2 != 7 * k + 1) {
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
    i2 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i2, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i2 = 0; i2 < k; i2++) {
      b_comm->data[b_comm->size[0] * i2] = comm->type->data[comm->type->size[0] *
        i2];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    b_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i2 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i2, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i2 = 0; i2 < k; i2++) {
    data->data[i2] = comm->data->data[i2];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  d_prodAx(A, x, b, nthreads, c_comm);
  emxFree_uint8_T(&data);
}

void prodAx_mpip(const emxArray_real_T *A, const emxArray_real_T *x,
                 emxArray_real_T *b, const emxArray_int32_T *nthreads, const
                 struct_T *comm, const emxArray_real_T *pbmsg)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i3;
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
  e_prodAx(A, x, b, nthreads, c_comm, pbmsg);
  emxFree_uint8_T(&data);
}

void prodAx_mpip1(const emxArray_real_T *A, const emxArray_real_T *x,
                  emxArray_real_T *b, const emxArray_int32_T *nthreads, const
                  struct_T *comm, const emxArray_real_T *pbmsg, int32_T pbsz)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i4;
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
  f_prodAx(A, x, b, nthreads, c_comm, pbmsg, pbsz);
  emxFree_uint8_T(&data);
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
