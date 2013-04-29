#include "crs_prodAtAx.h"
#include "spalab_kernel.h"
#include "mpi.h"
#include "omp.h"
#include "m2c.h"

static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
static void accu_partsum(emxArray_real_T *b, int32_T ncols);
static void allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1);
static void b_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2);
static void b_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void b_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b);
static define_emxInit(b_emxInit, real_T, real_T)
static void b_get_local_chunk(int32_T m, int32_T varargin_2, int32_T *istart,
  int32_T *iend);
static void b_msg_error(void);
static void c_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3);
static void c_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, MPI_Comm
  varargin_1);
static void c_msg_error(void);
static void crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void crs_prodAtx_internal(const emxArray_int32_T *row_ptr, const
  emxArray_int32_T *col_ind, const emxArray_real_T *val, const emxArray_real_T
  *x, emxArray_real_T *b, int32_T offset, int32_T nrows, int32_T ncols,
  boolean_T ismt);
static void crs_prodAx(const emxArray_int32_T *A_row_ptr, const emxArray_int32_T
  *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows, const
  emxArray_real_T *x, emxArray_real_T *b);
static void crs_prodAx_internal(int32_T nrows, const emxArray_int32_T *row_ptr,
  const emxArray_int32_T *col_ind, const emxArray_real_T *val, const
  emxArray_real_T *x, emxArray_real_T *b);
static void d_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, emxArray_real_T
  *Ax, const emxArray_int32_T *nthreads, MPI_Comm varargin_1);
static void d_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2);
static void d_msg_error(void);
static void e_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, emxArray_real_T
  *Ax, const emxArray_int32_T *nthreads, MPI_Comm varargin_1, const
  emxArray_real_T *varargin_2);
static void e_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3);
static void e_msg_error(const emxArray_char_T *varargin_3);
static void f_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, emxArray_real_T
  *Ax, const emxArray_int32_T *nthreads, MPI_Comm varargin_1, const
  emxArray_real_T *varargin_2, int32_T varargin_3);
static void f_msg_error(void);
static void get_local_chunk(int32_T m, int32_T *istart, int32_T *iend);
static void msg_error(void);
static void msg_warn(void);
static real_T rt_roundd(real_T u);
static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce(sptr, rptr, count, datatype, op, comm);
}

static void accu_partsum(emxArray_real_T *b, int32_T ncols)
{
  int32_T iend;
  int32_T istart;
  int32_T nthreads;
  int32_T offset;
  int32_T j;
  int32_T i4;
  int32_T k;
  int32_T i;
  if (b->size[0] >= (ncols << 1)) {
    get_local_chunk(ncols, &istart, &iend);
    nthreads = omp_get_num_threads();
    offset = (int32_T)rt_roundd(floor((real_T)b->size[0] / (real_T)ncols));
    if (nthreads <= offset) {
    } else {
      nthreads = offset;
    }

    offset = ncols;
    for (j = 2; j <= nthreads; j++) {
      i4 = b->size[1];
      for (k = 0; k + 1 <= i4; k++) {
        for (i = istart - 1; i + 1 <= iend; i++) {
          b->data[i + b->size[0] * k] += b->data[(offset + i) + b->size[0] * k];
        }
      }

      offset += ncols;
    }
  }
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
    f_msg_error();
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

static void b_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T i0;
  int32_T A_ncols_idx_1;
  emxArray_real_T *Ax;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A_ncols;
  emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = x->size[1];
  i0 = b->size[0] * b->size[1];
  b->size[1] = A_ncols_idx_1;
  emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = A_ncols * x->size[1];
  for (i0 = 0; i0 < A_ncols_idx_1; i0++) {
    b->data[i0] = 0.0;
  }

  b_emxInit_real_T(&Ax, 2);
  i0 = Ax->size[0] * Ax->size[1];
  Ax->size[0] = A_nrows;
  emxEnsureCapacity((emxArray__common *)Ax, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = x->size[1];
  i0 = Ax->size[0] * Ax->size[1];
  Ax->size[1] = A_ncols_idx_1;
  emxEnsureCapacity((emxArray__common *)Ax, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = A_nrows * x->size[1];
  for (i0 = 0; i0 < A_ncols_idx_1; i0++) {
    Ax->data[i0] = 0.0;
  }

  b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
  b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  emxFree_real_T(&Ax);
}

static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T n;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, FALSE);
}

static void b_crs_prodAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  const emxArray_real_T *x, emxArray_real_T *b)
{
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    c_msg_error();
  }

  crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
}

static void b_get_local_chunk(int32_T m, int32_T varargin_2, int32_T *istart,
  int32_T *iend)
{
  int32_T thr_id;
  int32_T chunk;
  if (varargin_2 == 1) {
    *istart = 1;
    *iend = m;
  } else {
    thr_id = omp_get_thread_num();
    if (varargin_2 >= m) {
      *istart = thr_id + 1;
      *iend = thr_id + (thr_id < m);
    } else {
      chunk = (int32_T)rt_roundd(ceil((real_T)m / (real_T)varargin_2));
      *istart = thr_id * chunk + 1;
      thr_id = (*istart + chunk) - 1;
      if (m <= thr_id) {
        *iend = m;
      } else {
        *iend = thr_id;
      }
    }
  }
}

static void b_msg_error(void)
{
  const char * msgid;
  msgid = "prodAAtx:IncorrectBuffer";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer Ax has incorrect size.");
  } else {
    printf("Error %s\nBuffer Ax has incorrect size.", msgid);
  }
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
    f_msg_error();
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

static void c_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  emxArray_real_T *Ax;
  int32_T i5;
  int32_T A_nrows_idx_1;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  b_emxInit_real_T(&Ax, 2);
  i5 = Ax->size[0] * Ax->size[1];
  Ax->size[0] = A_nrows;
  emxEnsureCapacity((emxArray__common *)Ax, i5, (int32_T)sizeof(real_T));
  A_nrows_idx_1 = x->size[1];
  i5 = Ax->size[0] * Ax->size[1];
  Ax->size[1] = A_nrows_idx_1;
  emxEnsureCapacity((emxArray__common *)Ax, i5, (int32_T)sizeof(real_T));
  A_nrows_idx_1 = A_nrows * x->size[1];
  for (i5 = 0; i5 < A_nrows_idx_1; i5++) {
    Ax->data[i5] = 0.0;
  }

  b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
  b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  emxFree_real_T(&Ax);
}

static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, MPI_Comm
  varargin_1)
{
  int32_T n;
  boolean_T b3;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b3 = (n > 1);
  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, b3);
  if (b3) {

#pragma omp barrier

    accu_partsum(b, A_ncols);
  }

  if (b3) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]), MPI_SUM,
            varargin_1);

            M2C_END_REGION(/*omp single*/)

}

static void c_msg_error(void)
{
  const char * msgid;
  msgid = "crs_prodAx:BufferTooSmal";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer space for output b is too small.");
  } else {
    printf("Error %s\nBuffer space for output b is too small.", msgid);
  }
}

static void crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T n;
  boolean_T b1;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b1 = (n > 1);
  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, b1);
  if (b1) {

#pragma omp barrier

    accu_partsum(b, A_ncols);
  }
}

static void crs_prodAtx_internal(const emxArray_int32_T *row_ptr, const
  emxArray_int32_T *col_ind, const emxArray_real_T *val, const emxArray_real_T
  *x, emxArray_real_T *b, int32_T offset, int32_T nrows, int32_T ncols,
  boolean_T ismt)
{
  int32_T n;
  int32_T u1;
  int32_T iend;
  int32_T istart;
  int32_T i3;
  int32_T k;
  if (ismt) {
    n = omp_get_num_threads();
    u1 = (int32_T)rt_roundd(floor((real_T)b->size[0] / (real_T)ncols));
    if (n <= u1) {
    } else {
      n = u1;
    }

    b_get_local_chunk(nrows, n, &istart, &iend);
  } else {
    istart = 1;
    iend = nrows;
  }

  i3 = x->size[1];
  for (k = 0; k + 1 <= i3; k++) {
    n = offset + ncols;
    for (u1 = offset; u1 + 1 <= n; u1++) {
      b->data[u1 + b->size[0] * k] = 0.0;
    }

    for (n = istart - 1; n + 1 <= iend; n++) {
      SPL_daxpy(x->data[n + x->size[0] * k], &val->data[row_ptr->data[n] - 1],
                &b->data[offset + b->size[0] * k], &col_ind->data[row_ptr->
                data[n] - 1], row_ptr->data[n + 1] - row_ptr->data[n]);
    }
  }
}

static void crs_prodAx(const emxArray_int32_T *A_row_ptr, const emxArray_int32_T
  *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows, const
  emxArray_real_T *x, emxArray_real_T *b)
{
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    c_msg_error();
  }

  omp_get_num_threads();
  crs_prodAx_internal(A_nrows, A_row_ptr, A_col_ind, A_val, x, b);
}

static void crs_prodAx_internal(int32_T nrows, const emxArray_int32_T *row_ptr,
  const emxArray_int32_T *col_ind, const emxArray_real_T *val, const
  emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T iend;
  int32_T i;
  int32_T i2;
  int32_T k;
  get_local_chunk(nrows, &i, &iend);
  while (i <= iend) {
    i2 = x->size[1];
    for (k = 0; k + 1 <= i2; k++) {
      b->data[(i + b->size[0] * k) - 1] = SPL_ddot(&val->data[row_ptr->data[i -
        1] - 1], &col_ind->data[row_ptr->data[i - 1] - 1], &x->data[x->size[0] *
        k], row_ptr->data[i] - row_ptr->data[i - 1]);
    }

    i++;
  }
}

static void d_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, emxArray_real_T
  *Ax, const emxArray_int32_T *nthreads, MPI_Comm varargin_1)
{
  int32_T n;
  boolean_T b2;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A_nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    M2C_END_REGION(/*omp master*/)

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

    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    c_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1);

                  M2C_END_REGION(/*omp parallel*/)

  } else if (b2) {
    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    c_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1);
  } else {
    b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
    b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  }
}

static void d_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2)
{
  int32_T n;
  boolean_T b5;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b5 = (n > 1);
  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, b5);
  if (b5) {

#pragma omp barrier

    accu_partsum(b, A_ncols);
  }

  if (b5) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  b_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2);

              M2C_END_REGION(/*omp single*/)

}

static void d_msg_error(void)
{
  const char * msgid;
  msgid = "crs_prodAtx:BufferTooSmal";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer space for output b is too small.");
  } else {
    printf("Error %s\nBuffer space for output b is too small.", msgid);
  }
}

static void e_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, emxArray_real_T
  *Ax, const emxArray_int32_T *nthreads, MPI_Comm varargin_1, const
  emxArray_real_T *varargin_2)
{
  int32_T n;
  boolean_T b4;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A_nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  b4 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b4 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    d_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2);

                  M2C_END_REGION(/*omp parallel*/)

  } else if (b4) {
    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    d_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2);
  } else {
    b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
    b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  }
}

static void e_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b7;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b7 = (n > 1);
  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, b7);
  if (b7) {

#pragma omp barrier

    accu_partsum(b, A_ncols);
  }

  if (b7) {

#pragma omp barrier

  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  c_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2, varargin_3);

              M2C_END_REGION(/*omp single*/)

}

static void e_msg_error(const emxArray_char_T *varargin_3)
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

static void f_crs_prodAtAx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, emxArray_real_T
  *Ax, const emxArray_int32_T *nthreads, MPI_Comm varargin_1, const
  emxArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b6;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A_nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  b6 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b6 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    e_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2, varargin_3);

                  M2C_END_REGION(/*omp parallel*/)

  } else if (b6) {
    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    e_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2, varargin_3);
  } else {
    b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
    b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  }
}

static void f_msg_error(void)
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

static void get_local_chunk(int32_T m, int32_T *istart, int32_T *iend)
{
  int32_T nthreads;
  int32_T thr_id;
  nthreads = omp_get_num_threads();
  if (nthreads == 1) {
    *istart = 1;
    *iend = m;
  } else {
    thr_id = omp_get_thread_num();
    if (nthreads >= m) {
      *istart = thr_id + 1;
      *iend = thr_id + (thr_id < m);
    } else {
      nthreads = (int32_T)rt_roundd(ceil((real_T)m / (real_T)nthreads));
      *istart = thr_id * nthreads + 1;
      nthreads = (*istart + nthreads) - 1;
      if (m <= nthreads) {
        *iend = m;
      } else {
        *iend = nthreads;
      }
    }
  }
}

static void msg_error(void)
{
  const char * msgid;
  msgid = "crs_prodAAtx:IncorrectBuffer";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer b has incorrect size.");
  } else {
    printf("Error %s\nBuffer b has incorrect size.", msgid);
  }
}

static void msg_warn(void)
{
  const char * msgid;
  msgid = "crs_prodAtAx:NestedParallel";
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

void crs_prodAtAx(const struct_T *A, const emxArray_real_T *x, emxArray_real_T
                  *b, emxArray_real_T *Ax, const emxArray_int32_T *nthreads)
{
  int32_T n;
  boolean_T b0;
  if ((b->size[0] < A->ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    msg_error();

    M2C_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A->nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    M2C_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    M2C_END_REGION(/*omp master*/)

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

    crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, Ax);

#pragma omp barrier

    crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, Ax, b);

    M2C_END_REGION(/*omp parallel*/)

  } else if (b0) {
    crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, Ax);

#pragma omp barrier

    crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, Ax, b);
  } else {
    b_crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, Ax);
    b_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, Ax, b);
  }
}

void crs_prodAtAx_initialize(void)
{
}

void crs_prodAtAx_mpi(const struct_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b, emxArray_real_T *Ax, const
                      emxArray_int32_T *nthreads, const b_struct_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i6;
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
      i6 = comm->type->size[k];
      if (i6 != 7 * k + 1) {
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
    i6 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i6, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i6 = 0; i6 < k; i6++) {
      b_comm->data[b_comm->size[0] * i6] = comm->type->data[comm->type->size[0] *
        i6];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    e_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i6 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i6, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i6 = 0; i6 < k; i6++) {
    data->data[i6] = comm->data->data[i6];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  d_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b, Ax,
                 nthreads, c_comm);
  emxFree_uint8_T(&data);
}

void crs_prodAtAx_mpip(const struct_T *A, const emxArray_real_T *x,
  emxArray_real_T *b, emxArray_real_T *Ax, const emxArray_int32_T *nthreads,
  const b_struct_T *comm, const emxArray_real_T *pbmsg)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i7;
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
      i7 = comm->type->size[k];
      if (i7 != 7 * k + 1) {
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
    i7 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i7, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i7 = 0; i7 < k; i7++) {
      b_comm->data[b_comm->size[0] * i7] = comm->type->data[comm->type->size[0] *
        i7];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    e_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i7 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i7, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i7 = 0; i7 < k; i7++) {
    data->data[i7] = comm->data->data[i7];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  e_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b, Ax,
                 nthreads, c_comm, pbmsg);
  emxFree_uint8_T(&data);
}

void crs_prodAtAx_mpip1(const struct_T *A, const emxArray_real_T *x,
  emxArray_real_T *b, emxArray_real_T *Ax, const emxArray_int32_T *nthreads,
  const b_struct_T *comm, const emxArray_real_T *pbmsg, int32_T pbsz)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i8;
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
      i8 = comm->type->size[k];
      if (i8 != 7 * k + 1) {
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
    i8 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i8, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i8 = 0; i8 < k; i8++) {
      b_comm->data[b_comm->size[0] * i8] = comm->type->data[comm->type->size[0] *
        i8];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    e_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i8 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i8, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i8 = 0; i8 < k; i8++) {
    data->data[i8] = comm->data->data[i8];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  f_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b, Ax,
                 nthreads, c_comm, pbmsg, pbsz);
  emxFree_uint8_T(&data);
}

void crs_prodAtAx_ser(const struct_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b)
{
  b_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtAx_ser1(const struct_T *A, const emxArray_real_T *x,
  emxArray_real_T *b)
{
  c_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtAx_terminate(void)
{
}

