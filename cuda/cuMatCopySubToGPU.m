function [errCode, toplevel] = cuMatCopySubToGPU ...
    (nrows, ncols, mat, cuMat)
% Copies a leading sub-matrix from MATLAB to CUDA.
%
% errCode = cuMatCopySubFromGPU(m, n, mat, cuMat)
% copies mat(1:nrows, 1:ncols) in MATLAB to cuMat(1:nrows, 1:ncols) on
% CUDA. This may be more efficient than creating a submatrix for mat.
%
% SEE ALSO cuMatCopyFromGPU

% Corresponding low-level cuBLAS function:
% cublasStatus_t cublasSetMatrix(int rows, int cols, int elemSize,
%           const void *A, int lda, void *B, int ldb)
% http://docs.nvidia.com/cuda/cublas/#cublassetmatrix

%#codegen -args {m2c_int, m2c_int, m2c_mat, CuMat}

coder.cinclude('mspack.h');

toplevel = nargout>1;

if (toplevel || m2c_debug) 
    if ~isfloat(mat)
        m2c_error('cuMatCopyToGPU:TypeMismatch', 'Expected floating-point numbers.');
    elseif  isreal(mat) && cuMat.type ~= CU_DOUBLE && cuMat.type ~= CU_SINGLE || ...
            ~isreal(mat) && cuMat.type ~= CU_DOUBLE_COMPLEX && cuMat.type ~= CU_COMPLEX
        m2c_error('cuMatCopyToGPU:TypeMismatch', 'Real and complex numbers mismatch.');
    elseif nrows>size(mat,1) || ncols>size(mat,1)
        m2c_error('cuMatCopyToGPU:SizeMismatch', 'Target matrix is too small.');
    elseif nrows>cuMat.dims(1) || ncols>cuMat.dims(2)
        m2c_error('cuMatCopyToGPU:SizeMismatch', 'Sourcd matrix is too small.');
    end
end

if cuMat.type == CU_DOUBLE || cuMat.type == CU_COMPLEX
    sizepe = int32(8);
elseif cuMat.type == CU_SINGLE
    sizepe = int32(4);
elseif cuMat.type == CU_DOUBLE_COMPLEX
    sizepe = int32(16);
else
    sizepe = int32(0);
    m2c_error('cuMatCopyToGPU:TypeMismatch', 'Expected real numbers.');
end

errCode = int32(0); %#ok<NASGU>
errCode = coder.ceval('cublasSetMatrix', nrows, ncols, sizepe, ...
    m2c_opaque_ptr(mat, 'void *'), int32(size(mat,1)), ...
    CuMat(cuMat, 'void *'), cuMat.dims(1));

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cublasSetMatrix returned error code %s\n', ...
        cuBlasGetErrorCode(errCode));
end

