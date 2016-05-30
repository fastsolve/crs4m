function [errCode, toplevel] = cuVecCopySubToGPU(n, vec, ...
    istart, inc, cuVec, incCuVec, varargin)
% Copies a sub-vector from MATLAB to CUDA
%
% cuVecCopyToGPU(n, vec, istart, inc, cuVec, incCuVec) copies
% vec(istart:inc,iend) in MATLAB to cuVec(1:incCuVec:) on CUDA (in MATLAB
% notation). This is more efficient than using vec(istart:stride:iend).
%
% cuVecCopyToGPU(n, vec, istart, inc, cuVec, incCuVec, 'async') copies
% data asynchronously (i.e., returning before finishing data copying).
%
% SEE ALSO cuVecCopyToGPU, cuVecCopyFromGPU

% Corresponding low-level cuBLAS function:
% cublasStatus_t cublasSetVector(int n, int elemSize, const void *x, ...
%     int incx, void *y, int incy)
% http://docs.nvidia.com/cuda/cublas/#cublassetvector

%#codegen -args {m2c_int, m2c_vec, m2c_int, m2c_int, CuVec, m2c_int}

coder.cinclude('mspack.h');

toplevel = nargout>1;

if ~isfloat(vec) && m2c_debug
    m2c_error('cuVecCopyToGPU:TypeMismatch', 'Expected floating-point numbers.');
elseif  (isreal(vec) && cuVec.type ~= CU_DOUBLE && cuVec.type ~= CU_SINGLE || ...
        ~isreal(vec) && cuVec.type ~= CU_DOUBLE_COMPLEX && cuVec.type ~= CU_COMPLEX) && ...
        (toplevel || m2c_debug)
    m2c_error('cuVecCopyToGPU:TypeMismatch', 'Real and complex numbers mismatch.');
elseif n>m2c_intdiv(cuVec.len,incCuVec) && (toplevel || m2c_debug)
    m2c_error('cuVecCopyToGPU:SizeMismatch', 'Target array is too small.');
end

if cuVec.type == CU_DOUBLE || cuVec.type == CU_COMPLEX
    sizepe = int32(8);
elseif cuVec.type == CU_SINGLE
    sizepe = int32(4);
elseif cuVec.type == CU_DOUBLE_COMPLEX
    sizepe = int32(16);
else
    sizepe = int32(0);
    m2c_error('cuVecCopyToGPU:TypeMismatch', 'Expected real numbers.');
end

errCode = int32(0); %#ok<NASGU>
errCode = coder.ceval('cublasSetVector', n, sizepe, ...
    m2c_opaque_ptr_const(vec, 'char *', (istart-1)*sizepe), inc, ...
    CuVec(cuVec, 'void *'), incCuVec);

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cublasSetVector returned error code %s\n', cuBlasGetErrorCode(errCode));
end
