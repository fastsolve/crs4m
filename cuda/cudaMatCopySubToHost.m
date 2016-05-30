function [mat, errCode, toplevel] = cudaMatCopySubToHost ...
    (nrows, ncols, cuMat, mat)
% Copies a leading sub-matrix from CUDA to MATLAB.
%
% [mat, errCode] = cudaMatCopySubToHost(m, n, cuMat, mat)
% copies cuMat(1:nrows, 1:ncols) to MATLAB mat(1:nrows, 1:ncols).
% This is more efficient than creating a submatrix for mat.
%
% SEE ALSO cudaMatCopyToHost

% Corresponding low-level cuBLAS function:
% cublasStatus_t cublasGetMatrix(int rows, int cols, int elemSize,
%           const void *A, int lda, void *B, int ldb)
% http://docs.nvidia.com/cuda/cublas/#cublasgetmatrix

%#codegen -args {m2c_int, m2c_int, CudaMat, m2c_mat}

coder.cinclude('mspack.h');

toplevel = nargout>2;

if (toplevel || m2c_debug) 
    if ~isfloat(mat)
        m2c_error('cudaMatCopyFromHost:TypeMismatch', 'Expected floating-point numbers.');
    elseif  isreal(mat) && cuMat.type ~= CUDA_DOUBLE && cuMat.type ~= CUDA_SINGLE || ...
            ~isreal(mat) && cuMat.type ~= CUDA_DOUBLE_COMPLEX && cuMat.type ~= CUDA_COMPLEX
        m2c_error('cudaMatCopyFromHost:TypeMismatch', 'Real and complex numbers mismatch.');
    elseif nrows>size(mat,1) || ncols>size(mat,1)
        m2c_error('cudaMatCopyFromHost:SizeMismatch', 'Source matrix is too small.');
    elseif nrows>cuMat.dims(1) || ncols>cuMat.dims(2)
        m2c_error('cudaMatCopyFromHost:SizeMismatch', 'Target matrix is too small.');
    end
end

if cuMat.type == CUDA_DOUBLE || cuMat.type == CUDA_COMPLEX
    sizepe = int32(8);
elseif cuMat.type == CUDA_SINGLE
    sizepe = int32(4);
elseif cuMat.type == CUDA_DOUBLE_COMPLEX
    sizepe = int32(16);
else
    sizepe = int32(0);
    m2c_error('cudaMatCopyFromHost:TypeMismatch', 'Expected real numbers.');
end

errCode = int32(0); %#ok<NASGU>
errCode = coder.ceval('cublasGetMatrix', nrows, ncols, sizepe, ...
    CudaMat(cuMat, 'void *'), cuMat.dims(1), ...
    m2c_opaque_ptr(mat, 'void *'), int32(size(mat,1)));

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cublasGetMatrix returned error code %s\n', ...
        cuBlasGetErrorCode(errCode));
end

