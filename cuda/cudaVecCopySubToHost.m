function [vec, errCode, toplevel] = cudaVecCopySubToHost(n, cuVec, ...
    incCuVec, vec, istart, inc)
% Copies a sub-vector from CUDA to MATLAB.
%
% [vec, errCode] = cudaVecCopySubToHost(n, cuVec, incCuVec, vec, istart, inc)
% copies cuVec(1:incCurVec:incCurVec*len) to MATLAB vec(istart:inc:).
% This is more efficient than creating a subvector for vec.
%
% SEE ALSO cudaVecCopyToHost

% Corresponding low-level cuBLAS function:
% cublasStatus_t cublasGetVector(int n, int elemSize,
%                const void *x, int incx, void *y, int incy)
% http://docs.nvidia.com/cuda/cublas/#cublasgetvector

%#codegen -args {m2c_int, CudaVec, m2c_int, m2c_vec, m2c_int, m2c_int}

coder.cinclude('mspack.h');

toplevel = nargout>2;

if toplevel || m2c_debug
    if ~isfloat(vec)
        m2c_error('cudaVecCopyFromHost:TypeMismatch', 'Expected floating-point numbers.');
    elseif  isreal(vec) && cuVec.type ~= CUDA_DOUBLE && cuVec.type ~= CUDA_SINGLE || ...
            ~isreal(vec) && cuVec.type ~= CUDA_DOUBLE_COMPLEX && cuVec.type ~= CUDA_COMPLEX
        m2c_error('cudaVecCopyFromHost:TypeMismatch', 'Real and complex numbers mismatch.');
    elseif n>m2c_intdiv(int32(length(vec)),inc)
        m2c_error('cudaVecCopyFromHost:SizeMismatch', 'Target array is too small.');
    end
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
errCode = coder.ceval('cublasGetVector', n, sizepe, ...
    CudaVec(cuVec, 'void *'), incCuVec, ...
    m2c_opaque_ptr(vec, 'char *', (istart-1)*sizepe), inc);

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cublasGetVector returned error code %s\n', cuBlasGetErrorCode(errCode));
end
