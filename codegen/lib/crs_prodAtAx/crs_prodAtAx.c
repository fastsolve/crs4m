#include "crs_prodAtAx.h"
#include "mpi.h"
#include "omp.h"
#include "palc.h"

static void allreduce(plcArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1);
static void b_allreduce(plcArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2);
static void b_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b);
static void b_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b);
static void b_crs_prodAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  const plcArray_real_T *x, plcArray_real_T *b);
static define_plcInit(b_plcInit_real_T, real_T);
static void b_msg_error(void);
static void c_allreduce(plcArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2, int32_T varargin_3);
static void c_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b);
static void c_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, MPI_Comm
  varargin_1);
static void c_msg_error(void);
static void crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b);
static void crs_prodAtx_kernel(const plcArray_int32_T *row_ptr, const
  plcArray_int32_T *col_ind, const plcArray_real_T *val, const plcArray_real_T
  *x, int32_T x_m, plcArray_real_T *b, int32_T b_m, int32_T nrows, int32_T ncols,
  int32_T nrhs, boolean_T ismt);
static void crs_prodAx(const plcArray_int32_T *A_row_ptr, const plcArray_int32_T
  *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows, const
  plcArray_real_T *x, plcArray_real_T *b);
static void crs_prodAx_kernel(const plcArray_int32_T *row_ptr, const
  plcArray_int32_T *col_ind, const plcArray_real_T *val, const plcArray_real_T
  *x, int32_T x_m, plcArray_real_T *b, int32_T b_m, int32_T nrows, int32_T nrhs,
  boolean_T ismt);
static void d_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, plcArray_real_T
  *Ax, const plcArray_int32_T *nthreads, MPI_Comm varargin_1);
static void d_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2);
static void d_msg_error(void);
static void e_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, plcArray_real_T
  *Ax, const plcArray_int32_T *nthreads, MPI_Comm varargin_1, const
  plcArray_real_T *varargin_2);
static void e_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2, int32_T varargin_3);
static void e_msg_error(const plcArray_char_T *varargin_3);
static void f_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, plcArray_real_T
  *Ax, const plcArray_int32_T *nthreads, MPI_Comm varargin_1, const
  plcArray_real_T *varargin_2, int32_T varargin_3);
static void f_msg_error(void);
static void get_local_chunk(int32_T m, int32_T varargin_2, int32_T *istart,
  int32_T *iend);
static void msg_error(void);
static void msg_warn(void);
static int32_T pMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
static real_T rt_roundd(real_T u);
static void allreduce(plcArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
                      varargin_1)
{
  int32_T size;
  void * ptr;
  MPI_Comm_size(varargin_1, &size);
  if (size > 1) {

    ptr = (void *)(&b->data[0]);
    pMPI_Allreduce(MPI_IN_PLACE, ptr, sz_in, MPI_DOUBLE, op, varargin_1);
  }
}

static void b_allreduce(plcArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2)
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

  MPI_Comm_size(varargin_1, &i);
  if (i > 1) {

    ptr = (void *)(&b->data[0]);
    pMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
  }
}

static void b_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b)
{
  int32_T i0;
  int32_T A_ncols_idx_1;
  plcArray_real_T *Ax;
  i0 = b->size[0] * b->size[1];
  b->size[0] = A_ncols;
  plcEnsureCapacity((plcArray__common *)b, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = x->size[1];
  i0 = b->size[0] * b->size[1];
  b->size[1] = A_ncols_idx_1;
  plcEnsureCapacity((plcArray__common *)b, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = A_ncols * x->size[1];
  for (i0 = 0; i0 < A_ncols_idx_1; i0++) {
    b->data[i0] = 0.0;
  }

  b_plcInit_real_T(&Ax, 2);
  i0 = Ax->size[0] * Ax->size[1];
  Ax->size[0] = A_nrows;
  plcEnsureCapacity((plcArray__common *)Ax, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = x->size[1];
  i0 = Ax->size[0] * Ax->size[1];
  Ax->size[1] = A_ncols_idx_1;
  plcEnsureCapacity((plcArray__common *)Ax, i0, (int32_T)sizeof(real_T));
  A_ncols_idx_1 = A_nrows * x->size[1];
  for (i0 = 0; i0 < A_ncols_idx_1; i0++) {
    Ax->data[i0] = 0.0;
  }

  b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
  b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  plcFree_real_T(&Ax);
}

static void b_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b)
{
  int32_T i7;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  i7 = b->size[0];
  crs_prodAtx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, i7, A_nrows, A_ncols, x->size[1], FALSE);
}

static void b_crs_prodAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  const plcArray_real_T *x, plcArray_real_T *b)
{
  int32_T i6;
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    c_msg_error();
  }

  i6 = b->size[0];
  crs_prodAx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, i6, A_nrows, x->size[1], FALSE);
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

static void c_allreduce(plcArray_real_T *b, int32_T sz_in, MPI_Op op, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2, int32_T varargin_3)
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

  MPI_Comm_size(varargin_1, &b_varargin_2);
  if (b_varargin_2 > 1) {

    ptr = (void *)(&b->data[0]);
    pMPI_Allreduce(MPI_IN_PLACE, ptr, s1, MPI_DOUBLE, op, varargin_1);
  }
}

static void c_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b)
{
  plcArray_real_T *Ax;
  int32_T i8;
  int32_T A_nrows_idx_1;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  b_plcInit_real_T(&Ax, 2);
  i8 = Ax->size[0] * Ax->size[1];
  Ax->size[0] = A_nrows;
  plcEnsureCapacity((plcArray__common *)Ax, i8, (int32_T)sizeof(real_T));
  A_nrows_idx_1 = x->size[1];
  i8 = Ax->size[0] * Ax->size[1];
  Ax->size[1] = A_nrows_idx_1;
  plcEnsureCapacity((plcArray__common *)Ax, i8, (int32_T)sizeof(real_T));
  A_nrows_idx_1 = A_nrows * x->size[1];
  for (i8 = 0; i8 < A_nrows_idx_1; i8++) {
    Ax->data[i8] = 0.0;
  }

  b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
  b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  plcFree_real_T(&Ax);
}

static void c_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, MPI_Comm
  varargin_1)
{
  int32_T n;
  boolean_T b2;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b2 = (n > 1);
  n = b->size[0];
  crs_prodAtx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, n, A_nrows, A_ncols, x->size[1], b2);
  if (b2) {

#pragma omp barrier

  }

#pragma omp single

  PLC_BEGIN_REGION(/*omp single*/)

  allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]), MPI_SUM,
            varargin_1);

  PLC_END_REGION(/*omp single*/)

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

static void crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b)
{
  int32_T n;
  int32_T i4;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  i4 = b->size[0];
  crs_prodAtx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, i4, A_nrows, A_ncols, x->size[1], n > 1);
}

static void crs_prodAtx_kernel(const plcArray_int32_T *row_ptr, const
  plcArray_int32_T *col_ind, const plcArray_real_T *val, const plcArray_real_T
  *x, int32_T x_m, plcArray_real_T *b, int32_T b_m, int32_T nrows, int32_T ncols,
  int32_T nrhs, boolean_T ismt)
{
  int32_T offset;
  int32_T xoffset;
  int32_T nthreads;
  int32_T boffset;
  int32_T iend;
  int32_T istart;
  int32_T k;
  int32_T i5;
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
      i5 = boffset + ncols;
      for (j = boffset; j + 1 <= i5; j++) {
        b->data[j] = 0.0;
      }

      for (i = istart; i <= iend; i++) {
        i5 = row_ptr->data[i] - 1;
        for (j = row_ptr->data[i - 1]; j <= i5; j++) {
          thr_id = (boffset + col_ind->data[j - 1]) - 1;
          b->data[thr_id] += x->data[i + xoffset] * val->data[j - 1];
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
        i5 = boffset + iend;
        for (i = boffset + istart; i + 1 <= i5; i++) {
          b->data[i] += b->data[offset + i];
        }

        boffset += b_m;
      }

      offset += ncols;
    }
  }
}

static void crs_prodAx(const plcArray_int32_T *A_row_ptr, const plcArray_int32_T
  *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows, const
  plcArray_real_T *x, plcArray_real_T *b)
{
  int32_T n;
  int32_T i2;
  if ((b->size[0] < A_nrows) || (b->size[1] < x->size[1])) {
    c_msg_error();
  }

  n = omp_get_num_threads();
  i2 = b->size[0];
  crs_prodAx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, i2, A_nrows, x->size[1], n > 1);
}

static void crs_prodAx_kernel(const plcArray_int32_T *row_ptr, const
  plcArray_int32_T *col_ind, const plcArray_real_T *val, const plcArray_real_T
  *x, int32_T x_m, plcArray_real_T *b, int32_T b_m, int32_T nrows, int32_T nrhs,
  boolean_T ismt)
{
  int32_T iend;
  int32_T istart;
  int32_T thr_id;
  real_T t;
  int32_T chunk;
  int32_T b_remainder;
  int32_T i;
  int32_T i3;
  int32_T j;

  if (ismt) {
    iend = omp_get_num_threads();
    if (iend == 1) {
      istart = 0;
      iend = nrows;
    } else {
      thr_id = omp_get_thread_num();
      t = (real_T)nrows / (real_T)iend;
      if (t < 0.0) {
        t = ceil(t);
      } else {
        t = floor(t);
      }

      chunk = (int32_T)rt_roundd(t);
      b_remainder = nrows - iend * chunk;
      if (b_remainder <= thr_id) {
        iend = b_remainder;
      } else {
        iend = thr_id;
      }

      istart = thr_id * chunk + iend;
      iend = (istart + chunk) + (thr_id < b_remainder);
    }
  } else {
    istart = 0;
    iend = nrows;
  }

  b_remainder = -1;
  thr_id = -1;
  for (chunk = 1; chunk <= nrhs; chunk++) {
    for (i = istart + 1; i <= iend; i++) {
      t = 0.0;
      i3 = row_ptr->data[i] - 1;
      for (j = row_ptr->data[i - 1]; j <= i3; j++) {
        t += val->data[j - 1] * x->data[b_remainder + col_ind->data[j - 1]];
      }

      b->data[thr_id + i] = t;
    }

    b_remainder += x_m;
    thr_id += b_m;
  }
}

static void d_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, plcArray_real_T
  *Ax, const plcArray_int32_T *nthreads, MPI_Comm varargin_1)
{
  int32_T n;
  boolean_T b1;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A_nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  b1 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b1 && (nthreads->data[0] > 1)) {

#pragma omp master

      PLC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      PLC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    PLC_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    c_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1);

    PLC_END_REGION(/*omp parallel*/)

  } else if (b1) {
    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    c_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1);
  } else {
    b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
    b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  }
}

static void d_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2)
{
  int32_T n;
  boolean_T b4;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b4 = (n > 1);
  n = b->size[0];
  crs_prodAtx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, n, A_nrows, A_ncols, x->size[1], b4);
  if (b4) {

#pragma omp barrier

  }

#pragma omp single

  PLC_BEGIN_REGION(/*omp single*/)

  b_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2);

  PLC_END_REGION(/*omp single*/)

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

static void e_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, plcArray_real_T
  *Ax, const plcArray_int32_T *nthreads, MPI_Comm varargin_1, const
  plcArray_real_T *varargin_2)
{
  int32_T n;
  boolean_T b3;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A_nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  b3 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b3 && (nthreads->data[0] > 1)) {

#pragma omp master

      PLC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      PLC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    PLC_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    d_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2);

    PLC_END_REGION(/*omp parallel*/)

  } else if (b3) {
    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    d_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2);
  } else {
    b_crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);
    b_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b);
  }
}

static void e_crs_prodAtx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, MPI_Comm
  varargin_1, const plcArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b6;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    d_msg_error();
  }

  n = omp_get_num_threads();
  b6 = (n > 1);
  n = b->size[0];
  crs_prodAtx_kernel(A_row_ptr, A_col_ind, A_val, x, x->size[0], b, n, A_nrows, A_ncols, x->size[1], b6);
  if (b6) {

#pragma omp barrier

  }

#pragma omp single

  PLC_BEGIN_REGION(/*omp single*/)

  c_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2, varargin_3);

  PLC_END_REGION(/*omp single*/)

}

static void e_msg_error(const plcArray_char_T *varargin_3)
{
  const char * msgid;
  plcArray_char_T *b_varargin_3;
  int32_T i1;
  int32_T loop_ub;
  msgid = "MPI_Comm:WrongType";
  plcInit_char_T(&b_varargin_3, 2);
  if (emlrtIsMATLABThread(emlrtRootTLSGlobal)) {
    i1 = b_varargin_3->size[0] * b_varargin_3->size[1];
    b_varargin_3->size[0] = 1;
    b_varargin_3->size[1] = varargin_3->size[1];
    plcEnsureCapacity((plcArray__common *)b_varargin_3, i1, (int32_T)sizeof
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
    plcEnsureCapacity((plcArray__common *)b_varargin_3, i1, (int32_T)sizeof
                      (char_T));
    loop_ub = varargin_3->size[0] * varargin_3->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_varargin_3->data[i1] = varargin_3->data[i1];
    }

    printf("Error %s\nIncorrect data type %s. Expected MPI_Datatype.", msgid,
           &b_varargin_3->data[0]);
  }

  plcFree_char_T(&b_varargin_3);
}

static void f_crs_prodAtAx(const plcArray_int32_T *A_row_ptr, const
  plcArray_int32_T *A_col_ind, const plcArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const plcArray_real_T *x, plcArray_real_T *b, plcArray_real_T
  *Ax, const plcArray_int32_T *nthreads, MPI_Comm varargin_1, const
  plcArray_real_T *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b5;
  if ((b->size[0] < A_ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A_nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  b5 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b5 && (nthreads->data[0] > 1)) {

#pragma omp master

      PLC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      PLC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    PLC_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx(A_row_ptr, A_col_ind, A_val, A_nrows, x, Ax);

#pragma omp barrier

    e_crs_prodAtx(A_row_ptr, A_col_ind, A_val, A_nrows, A_ncols, Ax, b,
                  varargin_1, varargin_2, varargin_3);

    PLC_END_REGION(/*omp parallel*/)

  } else if (b5) {
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

static int32_T pMPI_Allreduce(void * sptr, void * rptr, int32_T count,
  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce(sptr, rptr, count, datatype, op, comm);
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

void crs_prodAtAx(const struct_T *A, const plcArray_real_T *x, plcArray_real_T
                  *b, plcArray_real_T *Ax, const plcArray_int32_T *nthreads)
{
  int32_T n;
  boolean_T b0;
  if ((b->size[0] < A->ncols) || (b->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  if ((Ax->size[0] < A->nrows) || (Ax->size[1] != x->size[1])) {

#pragma omp master

    PLC_BEGIN_REGION(/*omp master*/)

    b_msg_error();

    PLC_END_REGION(/*omp master*/)

  }

  n = omp_get_num_threads();
  b0 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b0 && (nthreads->data[0] > 1)) {

#pragma omp master

      PLC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      PLC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    PLC_BEGIN_REGION(/*omp parallel*/)

    crs_prodAx(A->row_ptr, A->col_ind, A->val, A->nrows, x, Ax);

#pragma omp barrier

    crs_prodAtx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, Ax, b);

    PLC_END_REGION(/*omp parallel*/)

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

void crs_prodAtAx_mpi(const struct_T *A, const plcArray_real_T *x,
                      plcArray_real_T *b, plcArray_real_T *Ax, const
                      plcArray_int32_T *nthreads, const b_struct_T *comm)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i9;
  boolean_T exitg1;
  static const char_T cv0[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  plcArray_char_T *b_comm;
  plcArray_uint8_T *data;
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
    plcInit_char_T(&b_comm, 2);
    i9 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    plcEnsureCapacity((plcArray__common *)b_comm, i9, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i9 = 0; i9 < k; i9++) {
      b_comm->data[b_comm->size[0] * i9] = comm->type->data[comm->type->size[0] *
        i9];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    e_msg_error(b_comm);
    plcFree_char_T(&b_comm);
  }

  plcInit_uint8_T(&data, 2);
  i9 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  plcEnsureCapacity((plcArray__common *)data, i9, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i9 = 0; i9 < k; i9++) {
    data->data[i9] = comm->data->data[i9];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  d_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b, Ax,
                 nthreads, c_comm);
  plcFree_uint8_T(&data);
}

void crs_prodAtAx_mpip(const struct_T *A, const plcArray_real_T *x,
  plcArray_real_T *b, plcArray_real_T *Ax, const plcArray_int32_T *nthreads,
  const b_struct_T *comm, const plcArray_real_T *pbmsg)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i10;
  boolean_T exitg1;
  static const char_T cv1[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  plcArray_char_T *b_comm;
  plcArray_uint8_T *data;
  MPI_Comm c_comm;
  p = FALSE;
  b_p = FALSE;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i10 = comm->type->size[k];
      if (i10 != 7 * k + 1) {
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
    plcInit_char_T(&b_comm, 2);
    i10 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    plcEnsureCapacity((plcArray__common *)b_comm, i10, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i10 = 0; i10 < k; i10++) {
      b_comm->data[b_comm->size[0] * i10] = comm->type->data[comm->type->size[0]
        * i10];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    e_msg_error(b_comm);
    plcFree_char_T(&b_comm);
  }

  plcInit_uint8_T(&data, 2);
  i10 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  plcEnsureCapacity((plcArray__common *)data, i10, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i10 = 0; i10 < k; i10++) {
    data->data[i10] = comm->data->data[i10];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  e_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b, Ax,
                 nthreads, c_comm, pbmsg);
  plcFree_uint8_T(&data);
}

void crs_prodAtAx_mpip1(const struct_T *A, const plcArray_real_T *x,
  plcArray_real_T *b, plcArray_real_T *Ax, const plcArray_int32_T *nthreads,
  const b_struct_T *comm, const plcArray_real_T *pbmsg, int32_T pbsz)
{
  boolean_T p;
  boolean_T b_p;
  int32_T k;
  int32_T exitg2;
  int32_T i11;
  boolean_T exitg1;
  static const char_T cv2[8] = { 'M', 'P', 'I', '_', 'C', 'o', 'm', 'm' };

  plcArray_char_T *b_comm;
  plcArray_uint8_T *data;
  MPI_Comm c_comm;
  p = FALSE;
  b_p = FALSE;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i11 = comm->type->size[k];
      if (i11 != 7 * k + 1) {
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
    plcInit_char_T(&b_comm, 2);
    i11 = b_comm->size[0] * b_comm->size[1];
    b_comm->size[0] = 1;
    b_comm->size[1] = comm->type->size[1] + 1;
    plcEnsureCapacity((plcArray__common *)b_comm, i11, (int32_T)sizeof(char_T));
    k = comm->type->size[1];
    for (i11 = 0; i11 < k; i11++) {
      b_comm->data[b_comm->size[0] * i11] = comm->type->data[comm->type->size[0]
        * i11];
    }

    b_comm->data[b_comm->size[0] * comm->type->size[1]] = '\x00';
    e_msg_error(b_comm);
    plcFree_char_T(&b_comm);
  }

  plcInit_uint8_T(&data, 2);
  i11 = data->size[0] * data->size[1];
  data->size[0] = comm->data->size[0];
  data->size[1] = comm->data->size[1];
  plcEnsureCapacity((plcArray__common *)data, i11, (int32_T)sizeof(uint8_T));
  k = comm->data->size[0] * comm->data->size[1];
  for (i11 = 0; i11 < k; i11++) {
    data->data[i11] = comm->data->data[i11];
  }

  c_comm = *(MPI_Comm*)(&data->data[0]);
  f_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b, Ax,
                 nthreads, c_comm, pbmsg, pbsz);
  plcFree_uint8_T(&data);
}

void crs_prodAtAx_ser(const struct_T *A, const plcArray_real_T *x,
                      plcArray_real_T *b)
{
  b_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtAx_ser1(const struct_T *A, const plcArray_real_T *x,
  plcArray_real_T *b)
{
  c_crs_prodAtAx(A->row_ptr, A->col_ind, A->val, A->nrows, A->ncols, x, b);
}

void crs_prodAtAx_terminate(void)
{
}

