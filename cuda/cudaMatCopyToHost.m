function [mat, errCode, toplevel] = cudaMatCopyToHost(cuMat, mat)
% Copies a matrix from CUDA to MATLAB.
%
% [mat, errCode] = cudaMatCopyToHost(cuMat, mat) copies the whole matrix
% from cuMat on CUDA device to mat in MATLAB. The size to be copied
% is equal to cuMat.dims.
%
% [mat, errCode] = cudaMatCopyToHost(cuMat) allocates mat in MATLAB
% and then copies cuMat on CUDA device to it. In code-generation mode,
% cuMat.type must be a constant at compile time.
%
% SEE ALSO cudaMatCopySubToHost, cudaMatCopySubFromHost

if nargin==1
    if cuMat.type==CUDA_SINGLE
        mat = zeros(cuMat.dims, 'single');
    elseif cuMat.type==CUDA_DOUBLE_COMPLEX
        mat = zeros(cuMat.dims, 'like', 1i);
    elseif cuMat.type==CUDA_COMPLEX
        mat = zeros(cuMat.dims, 'like', single(1i));
    else
        mat = zeros(cuMat.dims);
    end
end

toplevel = nargout>2;
[mat, errCode] = cudaMatCopySubToHost(cuMat.dims(1), cuMat.dims(2), ...
    cuMat, mat);
