#include "vecDot.h"
#include "momp.h"
#include "mspack.h"
#include "m2c.h"

static void b_m2c_error(void);
static void b_m2c_warn(void);
static double b_vecDot(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *buf, int varargin_1);
static void b_vecDot_partial(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *prod, int nthreads, int n);
static void c_m2c_error(void);
static void c_vecDot_partial(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *prod, int nthreads, int n);
static void cuBlasGetErrorCode(int errCode, emxArray_char_T *cstr);
static void d_m2c_error(const emxArray_char_T *varargin_3);
static void e_m2c_error(const emxArray_char_T *varargin_3);
static void emxFreeStruct_struct1_T(struct1_T *pStruct);
static void emxInitStruct_struct1_T(struct1_T *pStruct);
static void f_m2c_error(void);
static void g_m2c_error(const emxArray_char_T *varargin_3);
static boolean_T isequal(const emxArray_char_T *varargin_1);
static void m2c_error(void);
static void m2c_warn(void);
static void vecDot_partial(const emxArray_real_T *x, const emxArray_real_T *y,
  double *prod, int nthreads, int n);
static void b_m2c_error(void)
{
  M2C_error("vecDot:BufferTooSmall",
            "Array buf must hold one entry per threads.\n");
}

static void b_m2c_warn(void)
{
  M2C_warn("vecDot:WrongInput", "Wrong mode name. cuda is assumed.\n");
}

static double b_vecDot(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *buf, int varargin_1)
{
  double prod;
  int flag;
  int u1;
  int nthreads;
  emxArray_real_T *b_buf;
  boolean_T guard1 = false;
  int j;
  int i4;
  int m;
  if (x->size[0] != y->size[0]) {
    flag = (M2C_DEBUG);
    if (flag != 0) {
      m2c_error();
    }
  }

  if (!((buf->size[0] == 0) || (buf->size[1] == 0))) {
    if ((buf->size[0] == 0) || (buf->size[1] == 0)) {
      u1 = 0;
    } else {
      flag = buf->size[0];
      u1 = buf->size[1];
      if (flag >= u1) {
        u1 = flag;
      }
    }

    if (u1 < varargin_1) {
      flag = (M2C_DEBUG);
      if (flag != 0) {
        b_m2c_error();
      }
    }
  }

  flag = omp_get_max_threads();
  if (varargin_1 <= flag) {
    nthreads = varargin_1;
  } else {
    nthreads = flag;
  }

  emxInit_real_T(&b_buf, 1);
  guard1 = false;
  if ((buf->size[0] == 0) || (buf->size[1] == 0)) {
    flag = omp_get_num_threads();
    if (flag > 1) {
      b_vecDot_partial(x, y, buf, nthreads, x->size[0]);

#if defined(M2C_OPENMP)
#pragma omp barrier
#endif
#if defined(M2C_OPENMP)
#pragma omp master
#endif

      M2C_BEGIN_REGION(/*omp master*/)

      for (j = 2; j <= nthreads; j++) {
        buf->data[0] += buf->data[j - 1];
      }

      prod = buf->data[0];

      M2C_END_REGION(/*omp master*/)

#if defined(M2C_OPENMP)
#pragma omp barrier
#endif

    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    if (nthreads > 1) {
      flag = omp_get_nested();
      if (!(flag != 0)) {
        flag = omp_get_num_threads();
        if (flag > 1) {
          flag = (M2C_DEBUG);
          if (flag != 0) {

#if defined(M2C_OPENMP)
#pragma omp master
#endif

            M2C_BEGIN_REGION(/*omp master*/)

            m2c_warn();

            M2C_END_REGION(/*omp master*/)

          }
        }
      }

      i4 = b_buf->size[0];
      b_buf->size[0] = nthreads;
      emxEnsureCapacity((emxArray__common *)b_buf, i4, (int)sizeof(double));
      for (i4 = 0; i4 < nthreads; i4++) {
        b_buf->data[i4] = 0.0;
      }

      m = x->size[0];

#if defined(M2C_OPENMP)
#pragma omp parallel default(none) num_threads(nthreads) shared(x, y, b_buf, nthreads, m)
#endif
      M2C_BEGIN_REGION(/*omp parallel*/)

      c_vecDot_partial(x, y, b_buf, nthreads, x->size[0]);

      M2C_END_REGION(/*omp parallel*/)

      for (j = 2; j <= nthreads; j++) {
        b_buf->data[0] += b_buf->data[j - 1];
      }

      prod = b_buf->data[0];
      flag = b_buf->size[0];
      i4 = buf->size[0] * buf->size[1];
      buf->size[0] = flag;
      buf->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)buf, i4, (int)sizeof(double));
      for (i4 = 0; i4 < 1; i4++) {
        for (u1 = 0; u1 < flag; u1++) {
          buf->data[u1] = b_buf->data[u1];
        }
      }
    } else {
      prod = 0.0;
      vecDot_partial(x, y, &prod, 1, x->size[0]);
    }
  }

  emxFree_real_T(&b_buf);
  return prod;
}

static void b_vecDot_partial(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *prod, int nthreads, int n)
{
  int threadId;
  int istart;
  int chunk;
  int iend;
  int b_remainder;
  boolean_T guard1 = false;
  int c_remainder;
  int flag;
  double* ptr;
  int i;
  double* b_ptr;
  if (nthreads > 1) {
    threadId = omp_get_thread_num();
    chunk = M2C_INTDIV(n, nthreads);
    b_remainder = n - nthreads * chunk;
    if (b_remainder <= threadId) {
      c_remainder = b_remainder;
    } else {
      c_remainder = threadId;
    }

    istart = threadId * chunk + c_remainder;
    iend = (istart + chunk) + (threadId < b_remainder);
  } else {
    threadId = 0;
    istart = 0;
    iend = n;
  }

  guard1 = false;
  if (iend >= istart + 1001) {
    flag = (M2C_BLAS);
    if (flag != 0) {
      ptr = (double*)(&x->data[0]);
      if (istart != 0) {
        ptr = M2C_OFFSET_PTR(ptr, istart);
      }

      b_ptr = (double*)(&y->data[0]);
      if (istart != 0) {
        b_ptr = M2C_OFFSET_PTR(b_ptr, istart);
      }

      prod->data[threadId] = cblas_ddot(iend - istart, ptr, 1, b_ptr, 1);
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    prod->data[threadId] = 0.0;
    for (i = istart; i + 1 <= iend; i++) {
      prod->data[threadId] += x->data[i] * y->data[i];
    }
  }
}

static void c_m2c_error(void)
{
  M2C_error("vecDot:WrongInput",
            "When using mex functions, vecDot only supports double.\n");
}

static void c_vecDot_partial(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *prod, int nthreads, int n)
{
  int threadId;
  int istart;
  int chunk;
  int iend;
  int b_remainder;
  boolean_T guard1 = false;
  int c_remainder;
  int flag;
  double* ptr;
  int i;
  double* b_ptr;
  if (nthreads > 1) {
    threadId = omp_get_thread_num();
    chunk = M2C_INTDIV(n, nthreads);
    b_remainder = n - nthreads * chunk;
    if (b_remainder <= threadId) {
      c_remainder = b_remainder;
    } else {
      c_remainder = threadId;
    }

    istart = threadId * chunk + c_remainder;
    iend = (istart + chunk) + (threadId < b_remainder);
  } else {
    threadId = 0;
    istart = 0;
    iend = n;
  }

  guard1 = false;
  if (iend >= istart + 1001) {
    flag = (M2C_BLAS);
    if (flag != 0) {
      ptr = (double*)(&x->data[0]);
      if (istart != 0) {
        ptr = M2C_OFFSET_PTR(ptr, istart);
      }

      b_ptr = (double*)(&y->data[0]);
      if (istart != 0) {
        b_ptr = M2C_OFFSET_PTR(b_ptr, istart);
      }

      prod->data[threadId] = cblas_ddot(iend - istart, ptr, 1, b_ptr, 1);
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    prod->data[threadId] = 0.0;
    for (i = istart; i + 1 <= iend; i++) {
      prod->data[threadId] += x->data[i] * y->data[i];
    }
  }
}

static void cuBlasGetErrorCode(int errCode, emxArray_char_T *cstr)
{
  int varargin_1;
  int varargin_2;
  int varargin_3;
  int varargin_4;
  int varargin_5;
  int varargin_6;
  int varargin_7;
  int varargin_8;
  static const char cv1[14] = { 'U', 'n', 'k', 'n', 'o', 'w', 'n', ' ', 'e', 'r',
    'r', 'o', 'r', '\x00' };

  static const char cv2[22] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'S', 'U', 'C', 'C', 'E', 'S', 'S', '\x00' };

  static const char cv3[30] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'N', 'O', 'T', '_', 'I', 'N', 'I', 'T', 'I', 'A', 'L',
    'I', 'Z', 'E', 'D', '\x00' };

  static const char cv4[27] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'A', 'L', 'L', 'O', 'C', '_', 'F', 'A', 'I', 'L', 'E',
    'D', '\x00' };

  static const char cv5[28] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'I', 'N', 'V', 'A', 'L', 'I', 'D', '_', 'V', 'A', 'L',
    'U', 'E', '\x00' };

  static const char cv6[28] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'A', 'R', 'C', 'H', '_', 'M', 'I', 'S', 'M', 'A', 'T',
    'C', 'H', '\x00' };

  static const char cv7[28] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'M', 'A', 'P', 'P', 'I', 'N', 'G', '_', 'E', 'R', 'R',
    'O', 'R', '\x00' };

  static const char cv8[31] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'E', 'X', 'E', 'C', 'U', 'T', 'I', 'O', 'N', '_', 'F',
    'A', 'I', 'L', 'E', 'D', '\x00' };

  static const char cv9[29] = { 'C', 'U', 'B', 'L', 'A', 'S', '_', 'S', 'T', 'A',
    'T', 'U', 'S', '_', 'I', 'N', 'T', 'E', 'R', 'N', 'A', 'L', '_', 'E', 'R',
    'R', 'O', 'R', '\x00' };

  varargin_1 = (CUBLAS_STATUS_SUCCESS);
  varargin_2 = (CUBLAS_STATUS_NOT_INITIALIZED);
  varargin_3 = (CUBLAS_STATUS_ALLOC_FAILED);
  varargin_4 = (CUBLAS_STATUS_INVALID_VALUE);
  varargin_5 = (CUBLAS_STATUS_ARCH_MISMATCH);
  varargin_6 = (CUBLAS_STATUS_MAPPING_ERROR);
  varargin_7 = (CUBLAS_STATUS_EXECUTION_FAILED);
  varargin_8 = (CUBLAS_STATUS_INTERNAL_ERROR);
  if (varargin_1 == errCode) {
    varargin_1 = 0;
  } else if (varargin_2 == errCode) {
    varargin_1 = 1;
  } else if (varargin_3 == errCode) {
    varargin_1 = 2;
  } else if (varargin_4 == errCode) {
    varargin_1 = 3;
  } else if (varargin_5 == errCode) {
    varargin_1 = 4;
  } else if (varargin_6 == errCode) {
    varargin_1 = 5;
  } else if (varargin_7 == errCode) {
    varargin_1 = 6;
  } else if (varargin_8 == errCode) {
    varargin_1 = 7;
  } else {
    varargin_1 = -1;
  }

  switch (varargin_1) {
   case 0:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 22;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 22; varargin_1++) {
      cstr->data[varargin_1] = cv2[varargin_1];
    }
    break;

   case 1:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 30;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 30; varargin_1++) {
      cstr->data[varargin_1] = cv3[varargin_1];
    }
    break;

   case 2:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 27;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 27; varargin_1++) {
      cstr->data[varargin_1] = cv4[varargin_1];
    }
    break;

   case 3:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 28;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 28; varargin_1++) {
      cstr->data[varargin_1] = cv5[varargin_1];
    }
    break;

   case 4:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 28;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 28; varargin_1++) {
      cstr->data[varargin_1] = cv6[varargin_1];
    }
    break;

   case 5:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 28;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 28; varargin_1++) {
      cstr->data[varargin_1] = cv7[varargin_1];
    }
    break;

   case 6:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 31;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 31; varargin_1++) {
      cstr->data[varargin_1] = cv8[varargin_1];
    }
    break;

   case 7:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 29;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 29; varargin_1++) {
      cstr->data[varargin_1] = cv9[varargin_1];
    }
    break;

   default:
    varargin_1 = cstr->size[0] * cstr->size[1];
    cstr->size[0] = 1;
    cstr->size[1] = 14;
    emxEnsureCapacity((emxArray__common *)cstr, varargin_1, (int)sizeof(char));
    for (varargin_1 = 0; varargin_1 < 14; varargin_1++) {
      cstr->data[varargin_1] = cv1[varargin_1];
    }
    break;
  }
}

static void d_m2c_error(const emxArray_char_T *varargin_3)
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

  M2C_error("m2c_opaque_obj:WrongInput",
            "Incorrect data type %s. Expected cublasHandle_t.\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

static void e_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int i2;
  int loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i2 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i2, (int)sizeof(char));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    b_varargin_3->data[i2] = varargin_3->data[i2];
  }

  M2C_error("CUDA:RuntimeError", "cublasGetPointerMode returned error code %s\n",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

static void emxFreeStruct_struct1_T(struct1_T *pStruct)
{
  emxFree_uint8_T(&pStruct->data);
  emxFree_char_T(&pStruct->type);
}

static void emxInitStruct_struct1_T(struct1_T *pStruct)
{
  emxInit_uint8_T(&pStruct->data, 1);
  emxInit_char_T(&pStruct->type, 2);
}

static void f_m2c_error(void)
{
  M2C_error("vecDot:WrongPointerMode",
            "The given cuBLAS handle has incorrect pointer mode.\n.");
}

static void g_m2c_error(const emxArray_char_T *varargin_3)
{
  emxArray_char_T *b_varargin_3;
  int i3;
  int loop_ub;
  emxInit_char_T(&b_varargin_3, 2);
  i3 = b_varargin_3->size[0] * b_varargin_3->size[1];
  b_varargin_3->size[0] = 1;
  b_varargin_3->size[1] = varargin_3->size[1];
  emxEnsureCapacity((emxArray__common *)b_varargin_3, i3, (int)sizeof(char));
  loop_ub = varargin_3->size[0] * varargin_3->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    b_varargin_3->data[i3] = varargin_3->data[i3];
  }

  M2C_error("cuBLAS:RuntimeError", "cuBLAS returned an error code %s\n.",
            &b_varargin_3->data[0]);
  emxFree_char_T(&b_varargin_3);
}

static boolean_T isequal(const emxArray_char_T *varargin_1)
{
  boolean_T p;
  boolean_T b_p;
  int k;
  int exitg2;
  int i0;
  boolean_T exitg1;
  static const char cv0[4] = { 'c', 'u', 'd', 'a' };

  p = false;
  b_p = false;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i0 = varargin_1->size[k];
      if (i0 != 3 * k + 1) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      b_p = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(varargin_1->size[1] == 0))) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k < 4)) {
      if (!(varargin_1->data[k] == cv0[k])) {
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

  return p;
}

static void m2c_error(void)
{
  M2C_error("vecDot:IncorrectSize",
            "Vectors x and y must have the same type and size.\n");
}

static void m2c_warn(void)
{
  M2C_warn("vecDot:NestedParallel",
           "You are trying to use nested parallel regions, but nested parallelism is not enabled.\n");
}

static void vecDot_partial(const emxArray_real_T *x, const emxArray_real_T *y,
  double *prod, int nthreads, int n)
{
  int istart;
  int threadId;
  int iend;
  int chunk;
  int b_remainder;
  boolean_T guard1 = false;
  int c_remainder;
  int flag;
  double* ptr;
  int i;
  double* b_ptr;
  if (nthreads > 1) {
    threadId = omp_get_thread_num();
    chunk = M2C_INTDIV(n, nthreads);
    b_remainder = n - nthreads * chunk;
    if (b_remainder <= threadId) {
      c_remainder = b_remainder;
    } else {
      c_remainder = threadId;
    }

    istart = threadId * chunk + c_remainder;
    iend = (istart + chunk) + (threadId < b_remainder);
  } else {
    istart = 0;
    iend = n;
  }

  guard1 = false;
  if (iend >= istart + 1001) {
    flag = (M2C_BLAS);
    if (flag != 0) {
      ptr = (double*)(&x->data[0]);
      if (istart != 0) {
        ptr = M2C_OFFSET_PTR(ptr, istart);
      }

      b_ptr = (double*)(&y->data[0]);
      if (istart != 0) {
        b_ptr = M2C_OFFSET_PTR(b_ptr, istart);
      }

      *prod = cblas_ddot(iend - istart, ptr, 1, b_ptr, 1);
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    *prod = 0.0;
    for (i = istart; i + 1 <= iend; i++) {
      *prod += x->data[i] * y->data[i];
    }
  }
}

void emxDestroy_struct1_T(struct1_T emxArray)
{
  emxFreeStruct_struct1_T(&emxArray);
}

void emxInitArray_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxInit_char_T(pEmxArray, numDimensions);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T(pEmxArray, numDimensions);
}

void emxInit_struct1_T(struct1_T *pStruct)
{
  emxInitStruct_struct1_T(pStruct);
}

void vecDot(const emxArray_real_T *x, const emxArray_real_T *y, const
            emxArray_real_T *buf, double *prod, boolean_T *toplevel)
{
  (void)buf;
  *toplevel = true;
  if (x->size[0] != y->size[0]) {
    m2c_error();
  }

  *prod = 0.0;
  vecDot_partial(x, y, prod, 1, x->size[0]);
}

double vecDot_cublas(const struct0_T *u, const struct0_T *v, const struct0_T
                     *buf, const emxArray_char_T *mode, const struct1_T
                     *cublasHdl)
{
  boolean_T p;
  boolean_T b_p;
  int k;
  int exitg4;
  int i5;
  boolean_T exitg3;
  emxArray_char_T *b_cublasHdl;
  static const char cv10[14] = { 'c', 'u', 'b', 'l', 'a', 's', 'H', 'a', 'n',
    'd', 'l', 'e', '_', 't' };

  emxArray_uint8_T *data;
  cublasHandle_t hdl;
  cublasPointerMode_t t_mode;
  int errCode;
  int b_mode;
  emxArray_char_T *r0;
  int exitg2;
  boolean_T exitg1;
  emxArray_char_T *c_cublasHdl;
  unsigned int b_data[2];
  double* output;
  double* b_output;
  double* c_output;
  if (!isequal(mode)) {
    b_m2c_warn();
  }

  if ((u->len != v->len) || (u->type != v->type)) {
    m2c_error();
  }

  if (u->type != 2) {
    c_m2c_error();
  }

  p = false;
  b_p = false;
  k = 0;
  do {
    exitg4 = 0;
    if (k < 2) {
      i5 = cublasHdl->type->size[k];
      if (i5 != 13 * k + 1) {
        exitg4 = 1;
      } else {
        k++;
      }
    } else {
      b_p = true;
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  if (b_p && (!(cublasHdl->type->size[1] == 0))) {
    k = 0;
    exitg3 = false;
    while ((!exitg3) && (k < 14)) {
      if (!(cublasHdl->type->data[k] == cv10[k])) {
        b_p = false;
        exitg3 = true;
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
    emxInit_char_T(&b_cublasHdl, 2);
    i5 = b_cublasHdl->size[0] * b_cublasHdl->size[1];
    b_cublasHdl->size[0] = 1;
    b_cublasHdl->size[1] = cublasHdl->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_cublasHdl, i5, (int)sizeof(char));
    k = cublasHdl->type->size[1];
    for (i5 = 0; i5 < k; i5++) {
      b_cublasHdl->data[b_cublasHdl->size[0] * i5] = cublasHdl->type->
        data[cublasHdl->type->size[0] * i5];
    }

    b_cublasHdl->data[b_cublasHdl->size[0] * cublasHdl->type->size[1]] = '\x00';
    d_m2c_error(b_cublasHdl);
    emxFree_char_T(&b_cublasHdl);
  }

  emxInit_uint8_T(&data, 1);
  i5 = data->size[0];
  data->size[0] = cublasHdl->data->size[0];
  emxEnsureCapacity((emxArray__common *)data, i5, (int)sizeof(unsigned char));
  k = cublasHdl->data->size[0];
  for (i5 = 0; i5 < k; i5++) {
    data->data[i5] = cublasHdl->data->data[i5];
  }

  hdl = *(cublasHandle_t*)(&data->data[0]);
  errCode = cublasGetPointerMode(hdl, &t_mode);
  b_mode = (t_mode);
  emxInit_char_T(&r0, 2);
  if (errCode != 0) {
    k = (M2C_DEBUG);
    if (k != 0) {
      cuBlasGetErrorCode(errCode, r0);
      e_m2c_error(r0);
    }
  }

  if (errCode != 0) {
  } else {
    k = (CUBLAS_POINTER_MODE_DEVICE);
    if (b_mode != k) {
      errCode = -1;
    } else {
      p = false;
      b_p = false;
      k = 0;
      do {
        exitg2 = 0;
        if (k < 2) {
          i5 = cublasHdl->type->size[k];
          if (i5 != 13 * k + 1) {
            exitg2 = 1;
          } else {
            k++;
          }
        } else {
          b_p = true;
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (b_p && (!(cublasHdl->type->size[1] == 0))) {
        k = 0;
        exitg1 = false;
        while ((!exitg1) && (k < 14)) {
          if (!(cublasHdl->type->data[k] == cv10[k])) {
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
        emxInit_char_T(&c_cublasHdl, 2);
        i5 = c_cublasHdl->size[0] * c_cublasHdl->size[1];
        c_cublasHdl->size[0] = 1;
        c_cublasHdl->size[1] = cublasHdl->type->size[1] + 1;
        emxEnsureCapacity((emxArray__common *)c_cublasHdl, i5, (int)sizeof(char));
        k = cublasHdl->type->size[1];
        for (i5 = 0; i5 < k; i5++) {
          c_cublasHdl->data[c_cublasHdl->size[0] * i5] = cublasHdl->type->
            data[cublasHdl->type->size[0] * i5];
        }

        c_cublasHdl->data[c_cublasHdl->size[0] * cublasHdl->type->size[1]] =
          '\x00';
        d_m2c_error(c_cublasHdl);
        emxFree_char_T(&c_cublasHdl);
      }

      i5 = data->size[0];
      data->size[0] = cublasHdl->data->size[0];
      emxEnsureCapacity((emxArray__common *)data, i5, (int)sizeof(unsigned char));
      k = cublasHdl->data->size[0];
      for (i5 = 0; i5 < k; i5++) {
        data->data[i5] = cublasHdl->data->data[i5];
      }

      hdl = *(cublasHandle_t*)(&data->data[0]);
      for (k = 0; k < 2; k++) {
        b_data[k] = u->data[k];
      }

      output = *(double**)(b_data);
      for (k = 0; k < 2; k++) {
        b_data[k] = v->data[k];
      }

      b_output = *(double**)(b_data);
      for (k = 0; k < 2; k++) {
        b_data[k] = buf->data[k];
      }

      c_output = *(double**)(b_data);
      errCode = cublasDdot(hdl, u->len, output, 1, b_output, 1, c_output);
    }
  }

  emxFree_uint8_T(&data);
  if (errCode != 0) {
    if (errCode < 0) {
      f_m2c_error();
    } else {
      cuBlasGetErrorCode(errCode, r0);
      g_m2c_error(r0);
    }
  }

  emxFree_char_T(&r0);
  return 0.0;
}

double vecDot_cublas_sync(const struct0_T *u, const struct0_T *v, const
  emxArray_real_T *buf, const emxArray_char_T *mode, const struct1_T *cublasHdl,
  const emxArray_char_T *sync)
{
  double prod;
  boolean_T p;
  boolean_T b_p;
  int k;
  int exitg4;
  int i6;
  boolean_T exitg3;
  emxArray_char_T *b_cublasHdl;
  static const char cv11[14] = { 'c', 'u', 'b', 'l', 'a', 's', 'H', 'a', 'n',
    'd', 'l', 'e', '_', 't' };

  emxArray_uint8_T *data;
  cublasHandle_t hdl;
  cublasPointerMode_t t_mode;
  int errCode;
  int b_mode;
  emxArray_char_T *r1;
  int exitg2;
  boolean_T exitg1;
  emxArray_char_T *c_cublasHdl;
  unsigned int b_data[2];
  double* output;
  double* b_output;
  (void)buf;
  (void)sync;
  if (!isequal(mode)) {
    b_m2c_warn();
  }

  if ((u->len != v->len) || (u->type != v->type)) {
    m2c_error();
  }

  if (u->type != 2) {
    c_m2c_error();
  }

  prod = 0.0;
  p = false;
  b_p = false;
  k = 0;
  do {
    exitg4 = 0;
    if (k < 2) {
      i6 = cublasHdl->type->size[k];
      if (i6 != 13 * k + 1) {
        exitg4 = 1;
      } else {
        k++;
      }
    } else {
      b_p = true;
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  if (b_p && (!(cublasHdl->type->size[1] == 0))) {
    k = 0;
    exitg3 = false;
    while ((!exitg3) && (k < 14)) {
      if (!(cublasHdl->type->data[k] == cv11[k])) {
        b_p = false;
        exitg3 = true;
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
    emxInit_char_T(&b_cublasHdl, 2);
    i6 = b_cublasHdl->size[0] * b_cublasHdl->size[1];
    b_cublasHdl->size[0] = 1;
    b_cublasHdl->size[1] = cublasHdl->type->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_cublasHdl, i6, (int)sizeof(char));
    k = cublasHdl->type->size[1];
    for (i6 = 0; i6 < k; i6++) {
      b_cublasHdl->data[b_cublasHdl->size[0] * i6] = cublasHdl->type->
        data[cublasHdl->type->size[0] * i6];
    }

    b_cublasHdl->data[b_cublasHdl->size[0] * cublasHdl->type->size[1]] = '\x00';
    d_m2c_error(b_cublasHdl);
    emxFree_char_T(&b_cublasHdl);
  }

  emxInit_uint8_T(&data, 1);
  i6 = data->size[0];
  data->size[0] = cublasHdl->data->size[0];
  emxEnsureCapacity((emxArray__common *)data, i6, (int)sizeof(unsigned char));
  k = cublasHdl->data->size[0];
  for (i6 = 0; i6 < k; i6++) {
    data->data[i6] = cublasHdl->data->data[i6];
  }

  hdl = *(cublasHandle_t*)(&data->data[0]);
  errCode = cublasGetPointerMode(hdl, &t_mode);
  b_mode = (t_mode);
  emxInit_char_T(&r1, 2);
  if (errCode != 0) {
    k = (M2C_DEBUG);
    if (k != 0) {
      cuBlasGetErrorCode(errCode, r1);
      e_m2c_error(r1);
    }
  }

  if (errCode != 0) {
  } else {
    k = (CUBLAS_POINTER_MODE_HOST);
    if (b_mode != k) {
      errCode = -1;
    } else {
      p = false;
      b_p = false;
      k = 0;
      do {
        exitg2 = 0;
        if (k < 2) {
          i6 = cublasHdl->type->size[k];
          if (i6 != 13 * k + 1) {
            exitg2 = 1;
          } else {
            k++;
          }
        } else {
          b_p = true;
          exitg2 = 1;
        }
      } while (exitg2 == 0);

      if (b_p && (!(cublasHdl->type->size[1] == 0))) {
        k = 0;
        exitg1 = false;
        while ((!exitg1) && (k < 14)) {
          if (!(cublasHdl->type->data[k] == cv11[k])) {
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
        emxInit_char_T(&c_cublasHdl, 2);
        i6 = c_cublasHdl->size[0] * c_cublasHdl->size[1];
        c_cublasHdl->size[0] = 1;
        c_cublasHdl->size[1] = cublasHdl->type->size[1] + 1;
        emxEnsureCapacity((emxArray__common *)c_cublasHdl, i6, (int)sizeof(char));
        k = cublasHdl->type->size[1];
        for (i6 = 0; i6 < k; i6++) {
          c_cublasHdl->data[c_cublasHdl->size[0] * i6] = cublasHdl->type->
            data[cublasHdl->type->size[0] * i6];
        }

        c_cublasHdl->data[c_cublasHdl->size[0] * cublasHdl->type->size[1]] =
          '\x00';
        d_m2c_error(c_cublasHdl);
        emxFree_char_T(&c_cublasHdl);
      }

      i6 = data->size[0];
      data->size[0] = cublasHdl->data->size[0];
      emxEnsureCapacity((emxArray__common *)data, i6, (int)sizeof(unsigned char));
      k = cublasHdl->data->size[0];
      for (i6 = 0; i6 < k; i6++) {
        data->data[i6] = cublasHdl->data->data[i6];
      }

      hdl = *(cublasHandle_t*)(&data->data[0]);
      for (k = 0; k < 2; k++) {
        b_data[k] = u->data[k];
      }

      output = *(double**)(b_data);
      for (k = 0; k < 2; k++) {
        b_data[k] = v->data[k];
      }

      b_output = *(double**)(b_data);
      errCode = cublasDdot(hdl, u->len, output, 1, b_output, 1, &prod);
    }
  }

  emxFree_uint8_T(&data);
  if (errCode != 0) {
    if (errCode < 0) {
      f_m2c_error();
    } else {
      cuBlasGetErrorCode(errCode, r1);
      g_m2c_error(r1);
    }
  }

  emxFree_char_T(&r1);
  return prod;
}

void vecDot_initialize(void)
{
}

double vecDot_omp(const emxArray_real_T *u, const emxArray_real_T *v,
                  emxArray_real_T *buf, int nthreads)
{
  return b_vecDot(u, v, buf, nthreads);
}

double vecDot_ser(const emxArray_real_T *u, const emxArray_real_T *v)
{
  double prod;
  int flag;
  if (u->size[0] != v->size[0]) {
    flag = (M2C_DEBUG);
    if (flag != 0) {
      m2c_error();
    }
  }

  prod = 0.0;
  vecDot_partial(u, v, &prod, 1, u->size[0]);
  return prod;
}

void vecDot_terminate(void)
{
}
