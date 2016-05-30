function [errCode, toplevel] = cudaVecCopySubFromHost(n, vec, ...
    istart, inc, cuVec, incCuVec)
% Copies a sub-vector from MATLAB to CUDA
%
% cudaVecCopyFromHost(n, vec, istart, inc, cuVec, incCuVec) copies
% vec(istart:inc,iend) in MATLAB to cuVec(1:incCuVec:) on CUDA (in MATLAB
% notation). This is more efficient than using vec(istart:stride:iend).
%
% SEE ALSO cudaVecCopyFromHost, cudaVecCopyToHost

% Corresponding low-level cuBLAS function:
% cublasStatus_t cublasSetVector(int n, int elemSize, const void *x, ...
%     int incx, void *y, int incy)
% http://docs.nvidia.com/cuda/cublas/#cublassetvector

%#codegen -args {m2c_int, m2c_vec, m2c_int, m2c_int, CudaVec, m2c_int}

coder.cinclude('mspack.h');

toplevel = nargout>1;

if ~isfloat(vec) && m2c_debug
    m2c_error('cudaVecCopyFromHost:TypeMismatch', 'Expected floating-point numbers.');
elseif  (isreal(vec) && cuVec.type ~= CUDA_DOUBLE && cuVec.type ~= CUDA_SINGLE || ...
        ~isreal(vec) && cuVec.type ~= CUDA_DOUBLE_COMPLEX && cuVec.type ~= CUDA_COMPLEX) && ...
        (toplevel || m2c_debug)
    m2c_error('cudaVecCopyFromHost:TypeMismatch', 'Real and complex numbers mismatch.');
elseif n>m2c_intdiv(cuVec.len,incCuVec) && (toplevel || m2c_debug)
    m2c_error('cudaVecCopyFromHost:SizeMismatch', 'Target array is too small.');
end

if cuVec.type == CUDA_DOUBLE || cuVec.type == CUDA_COMPLEX
    sizepe = int32(8);
elseif cuVec.type == CUDA_SINGLE
    sizepe = int32(4);
elseif cuVec.type == CUDA_DOUBLE_COMPLEX
    sizepe = int32(16);
else
    sizepe = int32(0);
    m2c_error('cudaVecCopyFromHost:TypeMismatch', 'Expected real numbers.');
end

errCode = int32(0); %#ok<NASGU>
errCode = coder.ceval('cublasSetVector', n, sizepe, ...
    m2c_opaque_ptr_const(vec, 'char *', (istart-1)*sizepe), inc, ...
    CudaVec(cuVec, 'void *'), incCuVec);

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cublasSetVector returned error code %s\n', cuBlasGetErrorCode(errCode));
end
