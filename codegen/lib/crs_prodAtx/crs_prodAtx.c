#include "crs_prodAtx.h"
#include "spalab_kernel.h"
#include "mpi.h"
#include "omp.h"
#include "stdio.h"
#define M2C_OFFSET_PTR(a,b)            ((char *)a)+(b)
#define MACC_BEGIN_REGION()            {
#define MACC_END_REGION()              }
#ifdef BUILD_MEX

#include "mex.h"
#define malloc                         mxMalloc
#define calloc                         mxCalloc
#define realloc                        mxRealloc
#define free                           mxFree
#define emlrtIsMATLABThread(s)         1
#define M2C_CHK_OPAQUE_PTR(ptr,parent,offset) \
 if ((parent) && (ptr) != ((char*)mxGetData(parent))+(offset)) \
 mexErrMsgIdAndTxt("opaque_ptr:ParentObjectChanged", \
 "The parent mxArray has changed. Avoid changing a MATLAB variable when dereferenced by an opaque_ptr.");
#else
#define emlrtIsMATLABThread(s)         0
#define mexErrMsgIdAndTxt(a,b)
#define M2C_CHK_OPAQUE_PTR(ptr,parent,offset)
#endif

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
static void accu_partsum(emxArray_real_T *b, int32_T ncols);
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
static void crs_prodAtx_internal(const emxArray_int32_T *row_ptr, const
  emxArray_int32_T *col_ind, const emxArray_real_T *val, const emxArray_real_T
  *x, emxArray_real_T *b, int32_T offset, int32_T nrows, int32_T ncols,
  boolean_T ismt);
static void d_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1);
static void e_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2);
static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize);
static void emxFree_char_T(emxArray_char_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray);
static void emxInit_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions);
static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int32_T numDimensions);
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

static void accu_partsum(emxArray_real_T *b, int32_T ncols)
{
  int32_T nthreads;
  int32_T thr_id;
  int32_T iend;
  int32_T offset;
  int32_T j;
  int32_T i2;
  int32_T k;
  int32_T i;
  if (b->size[0] >= (ncols << 1)) {
    nthreads = omp_get_num_threads();
    if (nthreads == 1) {
      thr_id = 0;
      iend = ncols;
    } else {
      thr_id = omp_get_thread_num();
      if (nthreads >= ncols) {
        iend = thr_id + (thr_id < ncols);
      } else {
        nthreads = (int32_T)rt_roundd(ceil((real_T)ncols / (real_T)nthreads));
        thr_id *= nthreads;
        nthreads += thr_id;
        if (ncols <= nthreads) {
          iend = ncols;
        } else {
          iend = nthreads;
        }
      }
    }

    offset = omp_get_num_threads();
    nthreads = (int32_T)rt_roundd(floor((real_T)b->size[0] / (real_T)ncols));
    if (offset <= nthreads) {
      nthreads = offset;
    }

    offset = ncols;
    for (j = 2; j <= nthreads; j++) {
      i2 = b->size[1];
      for (k = 0; k + 1 <= i2; k++) {
        for (i = thr_id; i + 1 <= iend; i++) {
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

static void b_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T n;
  n = b->size[0] * b->size[1];
  b->size[0] = A_ncols;
  b->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b, n, (int32_T)sizeof(real_T));
  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, FALSE);
}

static void b_msg_error(const emxArray_char_T *varargin_3)
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

static void c_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b)
{
  int32_T n;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_thread_num();
  crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                       A_ncols, FALSE);
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

static void crs_prodAtx_internal(const emxArray_int32_T *row_ptr, const
  emxArray_int32_T *col_ind, const emxArray_real_T *val, const emxArray_real_T
  *x, emxArray_real_T *b, int32_T offset, int32_T nrows, int32_T ncols,
  boolean_T ismt)
{
  int32_T n;
  int32_T u1;
  int32_T iend;
  int32_T istart;
  int32_T i1;
  int32_T k;
  if (ismt) {
    n = omp_get_num_threads();
    u1 = (int32_T)rt_roundd(floor((real_T)b->size[0] / (real_T)ncols));
    if (n <= u1) {
    } else {
      n = u1;
    }

    get_local_chunk(nrows, n, &istart, &iend);
  } else {
    istart = 1;
    iend = nrows;
  }

  i1 = x->size[1];
  for (k = 0; k + 1 <= i1; k++) {
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

static void d_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1)
{
  int32_T n;
  boolean_T b1;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b1 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b1 && (nthreads->data[0] > 1)) {

#pragma omp master

      MACC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      MACC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    MACC_BEGIN_REGION(/*omp parallel*/)

    n = omp_get_thread_num();
    crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                         A_ncols, b1);
    if (b1) {

#pragma omp barrier

      accu_partsum(b, A_ncols);
    }

    MACC_END_REGION(/*omp parallel*/)

  } else {
    n = omp_get_thread_num();
    crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                         A_ncols, b1);
    if (b1) {

#pragma omp barrier

      accu_partsum(b, A_ncols);
    }

    if (b1) {

#pragma omp barrier

    }
  }

#pragma omp single

  MACC_BEGIN_REGION(/*omp single*/)

  allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]), MPI_SUM,
            varargin_1);

            MACC_END_REGION(/*omp single*/)

}

static void e_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2)
{
  int32_T n;
  boolean_T b2;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b2 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b2 && (nthreads->data[0] > 1)) {

#pragma omp master

      MACC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      MACC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    MACC_BEGIN_REGION(/*omp parallel*/)

    n = omp_get_thread_num();
    crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                         A_ncols, b2);
    if (b2) {

#pragma omp barrier

      accu_partsum(b, A_ncols);
    }

    MACC_END_REGION(/*omp parallel*/)

  } else {
    n = omp_get_thread_num();
    crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                         A_ncols, b2);
    if (b2) {

#pragma omp barrier

      accu_partsum(b, A_ncols);
    }

    if (b2) {

#pragma omp barrier

    }
  }

#pragma omp single

  MACC_BEGIN_REGION(/*omp single*/)

  b_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2);

              MACC_END_REGION(/*omp single*/)

}

static void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel,
  int32_T elementSize)
{
  int32_T newNumel;
  int32_T i;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    if (emxArray->allocatedSize==0)
      i = newNumel;
    else {
      i = emxArray->allocatedSize;
      if (i < 16) {
        i = 16;
      }

      while (i < newNumel) {
        i <<= 1;
      }
    }

    newData = calloc((uint32_T)i, (uint32_T)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (uint32_T)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = TRUE;
  }
}


static void emxFree_char_T(emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxFree_uint8_T(emxArray_uint8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint8_T *)NULL) {
    if ((*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_uint8_T *)NULL;
  }
}

static void emxInit_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions)
{
  emxArray_char_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_char_T *)malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions)
{
  emxArray_int32_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int32_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void emxInit_uint8_T(emxArray_uint8_T **pEmxArray, int32_T numDimensions)
{
  emxArray_uint8_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_uint8_T *)malloc(sizeof(emxArray_uint8_T));
  emxArray = *pEmxArray;
  emxArray->data = (uint8_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = TRUE;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void f_crs_prodAtx(const emxArray_int32_T *A_row_ptr, const
  emxArray_int32_T *A_col_ind, const emxArray_real_T *A_val, int32_T A_nrows,
  int32_T A_ncols, const emxArray_real_T *x, emxArray_real_T *b, const
  emxArray_int32_T *nthreads, MPI_Comm varargin_1, const emxArray_real_T
  *varargin_2, int32_T varargin_3)
{
  int32_T n;
  boolean_T b3;
  if ((b->size[0] < A_ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b3 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b3 && (nthreads->data[0] > 1)) {

#pragma omp master

      MACC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      MACC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    MACC_BEGIN_REGION(/*omp parallel*/)

    n = omp_get_thread_num();
    crs_prodAtx_internal(A_row_ptr, A_col_ind, A_val, x, b, n * A_ncols, A_nrows,
                         A_ncols, b3);
    if (b3) {

#pragma omp barrier

      accu_partsum(b, A_ncols);
    }

    MACC_END_REGION(/*omp parallel*/)

  } else {
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
  }

#pragma omp single

  MACC_BEGIN_REGION(/*omp single*/)

  c_allreduce(b, (int32_T)rt_roundd((real_T)A_ncols * (real_T)x->size[1]),
              MPI_SUM, varargin_1, varargin_2, varargin_3);

              MACC_END_REGION(/*omp single*/)

}

static void get_local_chunk(int32_T m, int32_T varargin_2, int32_T *istart,
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
  if ((b->size[0] < A->ncols) || (b->size[1] < x->size[1])) {
    msg_error();
  }

  n = omp_get_num_threads();
  b0 = (n > 1);
  if (!(nthreads->size[0] == 0)) {
    n = omp_get_nested();
    if ((!(n != 0)) && b0 && (nthreads->data[0] > 1)) {

#pragma omp master

      MACC_BEGIN_REGION(/*omp master*/)

      msg_warn();

      MACC_END_REGION(/*omp master*/)

    }

#pragma omp parallel default(shared) num_threads(nthreads->data[0])
    MACC_BEGIN_REGION(/*omp parallel*/)

    n = omp_get_thread_num();
    crs_prodAtx_internal(A->row_ptr, A->col_ind, A->val, x, b, n * A->ncols,
                         A->nrows, A->ncols, b0);
    if (b0) {

#pragma omp barrier

      accu_partsum(b, A->ncols);
    }

    MACC_END_REGION(/*omp parallel*/)

  } else {
    n = omp_get_thread_num();
    crs_prodAtx_internal(A->row_ptr, A->col_ind, A->val, x, b, n * A->ncols,
                         A->nrows, A->ncols, b0);
    if (b0) {

#pragma omp barrier

      accu_partsum(b, A->ncols);
    }
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

emxArray_char_T *emxCreateND_char_T(int32_T numDimensions, int32_T *size)
{
  emxArray_char_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_char_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (char_T *)calloc((uint32_T)numEl, sizeof(char_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_int32_T *emxCreateND_int32_T(int32_T numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_uint8_T *emxCreateND_uint8_T(int32_T numDimensions, int32_T *size)
{
  emxArray_uint8_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_uint8_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (uint8_T *)calloc((uint32_T)numEl, sizeof(uint8_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_char_T *emxCreateWrapperND_char_T(char_T *data, int32_T numDimensions,
  int32_T *size)
{
  emxArray_char_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_char_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_int32_T *emxCreateWrapperND_int32_T(int32_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_int32_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_int32_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions,
  int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_uint8_T *emxCreateWrapperND_uint8_T(uint8_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_uint8_T *emx;
  int32_T numEl;
  int32_T i;
  emxInit_uint8_T(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
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

emxArray_char_T *emxCreate_char_T(int32_T rows, int32_T cols)
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

  emx->data = (char_T *)calloc((uint32_T)numEl, sizeof(char_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_int32_T *emxCreate_int32_T(int32_T rows, int32_T cols)
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

  emx->data = (int32_T *)calloc((uint32_T)numEl, sizeof(int32_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols)
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

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_uint8_T *emxCreate_uint8_T(int32_T rows, int32_T cols)
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

  emx->data = (uint8_T *)calloc((uint32_T)numEl, sizeof(uint8_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_char_T(emxArray_char_T *emxArray)
{
  emxFree_char_T(&emxArray);
}

void emxDestroyArray_int32_T(emxArray_int32_T *emxArray)
{
  emxFree_int32_T(&emxArray);
}

void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}

void emxDestroyArray_uint8_T(emxArray_uint8_T *emxArray)
{
  emxFree_uint8_T(&emxArray);
}
