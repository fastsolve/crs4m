function [cuMat, errCode] = cudaMatCopyFromHost(mat, cuMat)
% Copies a matrix from MATLAB to CUDA
%
% cudaMatCopyFromHost(mat, cuMat) copies the whole array from mat in MATLAB
% to cuMat on CUDA.
%
% [cuMat, errCode] = cudaMatCopyFromHost(mat) copies the whole array from
% mat in MATLAB to a new matrix cuMat on CUDA. This is a shorthand for
% calling cudaMatCreate followed by calling cudaMatCopyFromHost. The matrix
% cuMat must be destroied by calling cudaMatDestroy after use.
%
% SEE ALSO cudaMatCopyToHost, cudaMatCopySubFromHost, cudaMatDestroy

if size(mat,2)~=1 && (isempty(coder.target) || m2c_debug)
    m2c_warn('cudaMatCopyFromHost:NotMattrix', 'First input should be a column matrix.\n');
end

if isreal(mat)
    if isa(mat, 'double')
        type = CUDA_DOUBLE;
    else
        type = CUDA_SINGLE;
    end
elseif isa(mat, 'double')
    type = CUDA_DOUBLE_COMPLEX;
else
    type = CUDA_COMPLEX;
end

if nargin==1
    cuMat = cudaMatCreate(int32(size(mat,1)), int32(size(mat,2)), type);
end

errCode = cudaMatCopySubFromHost(int32(size(mat,1)), int32(size(mat,2)), mat, cuMat);
