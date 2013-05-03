#include "crs_prodAtx.h"
#include "mpi.h"
#include "omp.h"
#include "m2c.h"

static int32_T MMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
static void allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1);
static void b_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2);
static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void b_msg_error(const emxArray_char_T *varargin_3);
static void c_allreduce(emxArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const emxArray_real_T *varargin_2, int32_T varargin_3);
static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b);
static void c_msg_error(void);
static void crs_prodAtx_kernel(const int32_T *row_ptr, const
  int32_T *col_ind, const real_T *val, const real_T *x, int32_T x_m, real_T *b, int32_T b_m, int32_T nrows, int32_T ncols,
  int32_T nrhs, boolean_T ismt);
static void d_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1);
static void e_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2);
static void f_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2, int32_T varargin_3);
static void get_local_chunk(int32_T m, int32_T varargin_2, int32_T *istart,
  int32_T *iend);
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
  s1 = sz_in + varargin_2->size[0];
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

static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T i0;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A_ncols;
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
  crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, A_ncols, A_nrows, A_ncols, x->size[1], FALSE);
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

  MPI_Comm_size(varargin_1, &b_varargin_2);
  if (b_varargin_2 > 1) {

    ptr = (void *)(&b->data[0]);
    MMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
  }
}

static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T i4;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  i4 = b->size[0];
  crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i4, A_nrows, A_ncols, x->size[1], FALSE);
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

static void crs_prodAtx_kernel(const int32_T *row_ptr, const
  int32_T *col_ind, const real_T *val, const real_T *x, int32_T x_m, real_T *b, int32_T b_m, int32_T nrows, int32_T ncols,
  int32_T nrhs, boolean_T ismt)
{
  int32_T offset;
  int32_T xoffset;
  int32_T nthreads;
  int32_T boffset;
  int32_T iend;
  int32_T istart;
  int32_T k;
  int32_T i3;
  int32_T j;
  int32_T i;
  int32_T thr_id;
  real_T u;

  if (ismt) {
    offset = omp_get_num_threads();
    xoffset = (int32_T)rt_roundd(floor((real_T)b_m / (real_T)ncols));
    if (offset <= xoffset) {
      nthreads = offset;
    } else {
      nthreads = xoffset;
    }

    offset = omp_get_thread_num();
    boffset = offset * ncols;
    get_local_chunk(nrows, nthreads, &istart, &iend);
  } else {
    nthreads = 1;
    boffset = 0;
    istart = 1;
    iend = nrows;
  }

  if (istart <= iend) {
    xoffset = -1;
    for (k = 1; k <= nrhs; k++) {
      i3 = boffset + ncols;
      for (j = boffset; j + 1 <= i3; j++) {
        b[j] = 0.0;
      }

      for (i = istart; i <= iend; i++) {
        i3 = row_ptr[i] - 1;
        for (j = row_ptr[i - 1]; j <= i3; j++) {
          thr_id = (boffset + col_ind[j - 1]) - 1;
          b[thr_id] += x[i + xoffset] * val[j - 1];
        }
      }

      xoffset += x_m;
      boffset += b_m;
    }
  }

  if (nthreads > 1) {

#pragma omp barrier

    xoffset = omp_get_num_threads();
    if (xoffset == 1) {
      istart = 0;
      iend = ncols;
    } else {
      thr_id = omp_get_thread_num();
      u = (real_T)ncols / (real_T)xoffset;
      if (u < 0.0) {
        u = ceil(u);
      } else {
        u = floor(u);
      }

      boffset = (int32_T)rt_roundd(u);
      xoffset = ncols - xoffset * boffset;
      if (xoffset <= thr_id) {
        offset = xoffset;
      } else {
        offset = thr_id;
      }

      istart = thr_id * boffset + offset;
      iend = (istart + boffset) + (thr_id < xoffset);
    }

    offset = ncols;
    for (j = 2; j <= nthreads; j++) {
      boffset = 0;
      for (k = 1; k <= nrhs; k++) {
        i3 = boffset + iend;
        for (i = boffset + istart; i + 1 <= i3; i++) {
          b[i] += b[offset + i];
        }

        boffset += b_m;
      }

      offset += ncols;
    }
  }
}

static void d_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1)
{
  int32_T n;
  boolean_T b1;
  int32_T i6;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
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

    n = omp_get_num_threads();
    i6 = b->size[0];
    crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i6, A_nrows, A_ncols, x->size[1], n > 1);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    i6 = b->size[0];
    crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i6, A_nrows, A_ncols, x->size[1], b1);
    if (b1) {

#pragma omp barrier

    }
  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]), MPI_SUM,
            varargin_1);

  M2C_END_REGION(/*omp single*/)

}

static void e_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2)
{
  int32_T n;
  boolean_T b2;
  int32_T i8;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
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

    n = omp_get_num_threads();
    i8 = b->size[0];
    crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i8, A_nrows, A_ncols, x->size[1], n > 1);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    i8 = b->size[0];
    crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i8, A_nrows, A_ncols, x->size[1], b2);
    if (b2) {

#pragma omp barrier

    }
  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  b_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2);

  M2C_END_REGION(/*omp single*/)

}

static void f_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b3;
  int32_T i10;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b3 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b3 && (nthreads->data[0] > 1)) {

#pragma omp master

      M2C_BEGIN_REGION(/*omp master*/)

      msg_warn();

      M2C_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    M2C_BEGIN_REGION(/*omp parallel*/)

    n = omp_get_num_threads();
    i10 = b->size[0];
    crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i10, A_nrows, A_ncols, x->size[1], n > 1);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    i10 = b->size[0];
    crs_prodAtx_kernel(A_row_ptr->data, A_col_ind->data, A_val->data, x->data, x->size[0], b->data, i10, A_nrows, A_ncols, x->size[1], b3);
    if (b3) {

#pragma omp barrier

    }
  }

#pragma omp single

  M2C_BEGIN_REGION(/*omp single*/)

  c_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2, varargin_3);

  M2C_END_REGION(/*omp single*/)

}

static void get_local_chunk(int32_T m, int32_T varargin_2, int32_T *istart,
  int32_T *iend)
{
  int32_T thr_id;
  real_T u;
  int32_T chunk;
  int32_T b_remainder;
  int32_T y;
  if (varargin_2 == 1) {
    *istart = 1;
    *iend = m;
  } else {
    thr_id = omp_get_thread_num();
    u = (real_T)m / (real_T)varargin_2;
    if (u < 0.0) {
      u = ceil(u);
    } else {
      u = floor(u);
    }

    chunk = (int32_T)rt_roundd(u);
    b_remainder = m - varargin_2 * chunk;
    if (b_remainder <= thr_id) {
      y = b_remainder;
    } else {
      y = thr_id;
    }

    *istart = (thr_id * chunk + y) + 1;
    *iend = ((*istart + chunk) + (thr_id < b_remainder)) - 1;
  }
}

static void msg_error(void)
{
  const char * msgid;
  msgid = "crs_prodAtx:BufferTooSmal";
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    mexErrMsgIdAndTxt(msgid, "Buffer space for output b is too small.");
  } else {
    printf("Error %s\nBuffer space for output b is too small.", msgid);
  }
}

static void msg_warn(void)
{
  const char * msgid;
  msgid = "crs_prodAtx:NestedParallel";
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

void crs_prodAtx(const struct_T *A, const emxArray_real_T *x, emxArray_real_T *b,
                 const emxArray_int32_T *nthreads)
{
  int32_T n;
  boolean_T b0;
  int32_T i2;
  if ((b->size[0] < A->ncols) || (b->size[1] < x->size[1])) {
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

    n = omp_get_num_threads();
    i2 = b->size[0];
    crs_prodAtx_kernel(A->row_ptr->data, A->col_ind->data, A->val->data, x->data, x->size[0], b->data, i2, A->nrows, A->ncols, x->size[1], n > 1);

    M2C_END_REGION(/*omp parallel*/)

  } else {
    i2 = b->size[0];
    crs_prodAtx_kernel(A->row_ptr->data, A->col_ind->data, A->val->data, x->data, x->size[0], b->data, i2, A->nrows, A->ncols, x->size[1], b0);
  }
}

void crs_prodAtx_initialize(void)
{
}

void crs_prodAtx_mpi(const struct_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b, const emxArray_int32_T *nthreads, const
                     b_struct_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i5;
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
  d_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b,
                nthreads, c_comm);
  emxFree_uint8_T(&data);
}

void crs_prodAtx_mpip(const struct_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b, const emxArray_int32_T *nthreads,
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
    b_msg_error(b_comm);
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
  e_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b,
                nthreads, c_comm, pbmsg);
  emxFree_uint8_T(&data);
}

void crs_prodAtx_mpip1(const struct_T *A, const emxArray_real_T *x,
  emxArray_real_T *b, const emxArray_int32_T *nthreads, const b_struct_T *comm,
  const emxArray_real_T *pbmsg, int32_T pbsz)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i9;
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
      i9 = comm->type->size[k];
      if (i9 != 7 * k + 1) {
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
    i9 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_comm, i9, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i9 = 0; i9 < k; i9++) {
      b_comm->data[b_comm->size[0] * i9] = comm->type->data[comm->type->size[0] *
        i9];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    b_msg_error(b_comm);
    emxFree_char_T(&b_comm);
  }

  emxInit_uint8_T(&data, 2);
  i9 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  emxEnsureCapacity((emxArray__common *)data, i9, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i9 = 0; i9 < k; i9++) {
    data->data[i9] = comm->data->data[i9];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  f_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b,
                nthreads, c_comm, pbmsg, pbsz);
  emxFree_uint8_T(&data);
}

void crs_prodAtx_ser(const struct_T *A, const emxArray_real_T *x,
                     emxArray_real_T *b)
{
  b_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtx_ser1(const struct_T *A, const emxArray_real_T *x,
                      emxArray_real_T *b)
{
  c_crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtx_terminate(void)
{
}

