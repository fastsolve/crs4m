/*
 * Some kernel macros for sparse matrices
 */

#ifndef _MSP_KERNEL_
#define _MSP_KERNEL_

/* Set M2C_BLAS to false by default. It can be overwritten by m2c options
 * -blas. When enabled, user must specify the library for cblas. */
#ifndef M2C_BLAS  
#define M2C_BLAS    0
#endif

/* Set M2C_SPARSE_BLAS to false by default. It can be overwritten by m2c options
 * -spblas. When enabled, user must specify the library for spblas. */
#ifndef M2C_SPARSE_BLAS  
#define M2C_SPARSE_BLAS    0
#endif

/* Set M2C_CUDA to false by default. It can be overwritten by m2c options
 * -cuda. When enabled, user must specify the library for cublas. */
#ifndef M2C_CUDA  
#define M2C_CUDA    0
#endif

#if M2C_BLAS
#include "cblas.h"
#else

#define cblas_sdsdot(N, alpha, X, incX, Y, incY)     0
#define cblas_ddsdot(N, alpha, X, incX, Y, incY)     0

#define cblas_sdot(N, X, incX, Y, incY)              0
#define cblas_ddot(N, X, incX, Y, incY)              0

#define cblas_cdotu_sub(N, X, incX, Y, incY, dotu)
#define cblas_cdotc_sub(N, X, incX, Y, incY, dotc)

#define cblas_zdotu_sub(N, X, incX, Y, incY, dotu)
#define cblas_zdotc_sub(N, X, incX, Y, incY, dotc)

#endif

#if M2C_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cublas_v2.h"
#include "cusparse.h"

#else

/* CUDA runtime */
typedef int  cudaError_t;
typedef void *cudaStream_t;

inline cudaError_t cudaMalloc(void *a, int sz) { a = NULL; return -1; }
inline cudaError_t cudaFree(void * a) { return -1; }

inline const char *cudaGetErrorString(cudaError_t err) { return NULL; }

inline cudaError_t cudaStreamCreate(cudaStream_t *pStream)
{ pStream = NULL; return -1; }
inline cudaError_t cudaStreamDestroy(cudaStream_t strm) { return -1; }
inline cudaError_t cudaStreamSynchronize(cudaStream_t  strm) {return -1; }

/* cuBLAS */
typedef enum {CUBLAS_STATUS_SUCCESS, CUBLAS_STATUS_NOT_INITIALIZED,
        CUBLAS_STATUS_ALLOC_FAILED, CUBLAS_STATUS_INVALID_VALUE,
        CUBLAS_STATUS_ARCH_MISMATCH, CUBLAS_STATUS_MAPPING_ERROR,
        CUBLAS_STATUS_EXECUTION_FAILED, CUBLAS_STATUS_INTERNAL_ERROR} cublasStatus_t ;
typedef enum {CUBLAS_POINTER_MODE_HOST, CUBLAS_POINTER_MODE_DEVICE} cublasPointerMode_t;
typedef enum {CUBLAS_ATOMICS_NOT_ALLOWED, CUBLAS_ATOMICS_ALLOWED} cublasAtomicsMode_t;

typedef void *       cublasHandle_t;
typedef creal32_T    cuComplex;
typedef creal64_T    cuDoubleComplex;

inline cublasStatus_t cublasCreate(cublasHandle_t *hdl) { hdl=NULL; return -1; }
inline cublasStatus_t cublasDestroy(cublasHandle_t hdl) { return -1; }

inline cublasStatus_t cublasGetPointerMode(cublasHandle_t handle, 
        cublasPointerMode_t *mode)
{ *mode = CUBLAS_POINTER_MODE_HOST; return -1; }
inline cublasStatus_t cublasSetPointerMode(cublasHandle_t handle, 
        cublasPointerMode_t mode) { return -1; }

inline cublasStatus_t cublasGetAtomicsMode(cublasHandle_t handle, 
        cublasAtomicsMode_t *mode)
{ *mode = CUBLAS_ATOMICS_NOT_ALLOWED; return -1; }
inline cublasStatus_t cublasSetAtomicsMode(cublasHandle_t handle, 
        cublasAtomicsMode_t mode) { return -1; }

inline cublasStatus_t cublasSetStream(cublasHandle_t hdl, cudaStream_t  strm)
{return -1; }
inline cublasStatus_t cublasGetStream(cublasHandle_t  strm, cudaStream_t *p)
{ *p = NULL; return -1; }

inline cublasStatus_t cublasSetVector(int n, int elemSize, 
        const void *hostPtr, int incx, void *devicePtr, int incy)
{ return -1; }
inline cublasStatus_t cublasSetVectorAsync(int n, int elemSize, 
        const void *hostPtr, int incx, void *devicePtr, int incy, cudaStream_t stream)
{ return -1; }

inline cublasStatus_t cublasGetVector(int n, int elemSize, 
        const void *hostPtr, int incx, void *devicePtr, int incy) 
{ return -1; }
inline cublasStatus_t cublasGetVectorAsync(int n, int elemSize, 
        const void *hostPtr, int incx, void *devicePtr, int incy, cudaStream_t stream)
{ return -1; }

inline cublasStatus_t cublasGetMatrix(int rows, int cols, int elemSize, 
        const void *A, int lda, void *B, int ldb)
{ return -1; }
inline cublasStatus_t cublasGetMatrixAsync(int rows, int cols, int elemSize, 
        const void *A, int lda, void *B, int ldb, cudaStream_t stream)
{return -1; }
inline cublasStatus_t cublasSetMatrix(int rows, int cols, int elemSize, 
        const void *A, int lda, void *B, int ldb)
{return -1; }
inline cublasStatus_t cublasSetMatrixAsync(int rows, int cols, int elemSize, 
        const void *A, int lda, void *B, int ldb, cudaStream_t stream)
{return -1; }

inline int cublasSdot(cublasHandle_t handle, int n,
        const void *x, int incx, void *y, int incy, float *result)
{ *result = 0; return -1; }
inline int cublasDdot(cublasHandle_t handle, int n,
        const void *x, int incx, void *y, int incy, double *result)
{ *result = 0; return -1; }

#define cublasCdotu(handle, n, x, incx, y, incy, result)   -1
#define cublasCdotc(handle, n, x, incx, y, incy, result)   -1

#define cublasZdotu(handle, n, x, incx, y, incy, result)   -1
#define cublasZdotc(handle, n, x, incx, y, incy, result)   -1

/* cuSparse*/
typedef enum {CUSPARSE_STATUS_SUCCESS, CUSPARSE_STATUS_NOT_INITIALIZED,
        CUSPARSE_STATUS_ALLOC_FAILED, CUSPARSE_STATUS_INVALID_VALUE,
        CUSPARSE_STATUS_ARCH_MISMATCH, CUSPARSE_STATUS_MAPPING_ERROR,
        CUSPARSE_STATUS_EXECUTION_FAILED, CUSPARSE_STATUS_INTERNAL_ERROR, 
        CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED} cusparseStatus_t;
typedef enum {CUSPARSE_POINTER_MODE_HOST, 
        CUSPARSE_POINTER_MODE_DEVICE} cusparsePointerMode_t;

typedef void *       cusparseHandle_t;

inline cusparseStatus_t cusparseCreate(cusparseHandle_t *hdl) { hdl=NULL; return -1; }
inline cusparseStatus_t cusparseDestroy(cusparseHandle_t hdl) { return -1; }

inline cusparseStatus_t cusparseGetPointerMode(cusparseHandle_t handle, 
        cusparsePointerMode_t *mode)
{ *mode = CUSPARSE_POINTER_MODE_HOST; return -1; }
inline cusparseStatus_t cusparseSetPointerMode(cusparseHandle_t handle, 
        cusparsePointerMode_t mode)
{ return -1; }

#endif /* M2C_CUDA for cusparse */

#endif
