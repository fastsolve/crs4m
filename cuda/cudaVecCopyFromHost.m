function [cuVec, errCode] = cudaVecCopyFromHost(vec, cuVec)
% Copies a vector from MATLAB to CUDA
%
% cudaVecCopyFromHost(vec, cuVec) copies the whole array from vec in MATLAB
% to cuVec on CUDA.
%
% [cuVec, errCode] = cudaVecCopyFromHost(vec) copies the whole array from
% vec in MATLAB to a new vector cuVec on CUDA. This is a shorthand for
% calling cudaVecCreate followed by calling cudaVecCopyFromHost. The vector
% cuVec must be destroied by calling cudaVecDestroy after use.
%
% SEE ALSO cudaVecCopyToHost, cudaVecCopySubFromHost, cudaVecDestroy

if size(vec,2)~=1 && (isempty(coder.target) || m2c_debug)
    m2c_warn('cudaVecCopyFromHost:NotVector', 'First input should be a column vector.\n');
end

if isreal(vec)
    if isa(vec, 'double')
        type = CUDA_DOUBLE;
    else
        type = CUDA_SINGLE;
    end
elseif isa(vec, 'double')
    type = CUDA_DOUBLE_COMPLEX;
else
    type = CUDA_COMPLEX;
end

if nargin==1
    cuVec = cudaVecCreate(int32(length(vec)), type);
end

errCode = cudaVecCopySubFromHost(int32(length(vec)), vec, ...
    int32(1), int32(1), cuVec, int32(1));
