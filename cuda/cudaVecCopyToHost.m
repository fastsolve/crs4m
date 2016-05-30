function [vec, errCode, toplevel] = cudaVecCopyToHost(cuVec, vec)
% Copies a vector from CUDA to MATLAB.
%
% [vec, errCode] = cudaVecCopyToHost(cuVec, vec) copies the whole array
% from cuVec on CUDA device to vec in MATLAB. The length to be copied
% is equal to cuVec.len.
%
% [vec, errCode] = cudaVecCopyToHost(cuVec) allocates vec in MATLAB
% and then copies cuVec on CUDA device to it. In code-generation mode,
% cuVec.type must be a constant at compile time.
%
% SEE ALSO cudaVecCopySubToHost, cudaVecCopySubFromHost

if nargin==1
    if cuVec.type==CUDA_SINGLE
        vec = zeros(cuVec.len, 1, 'single');
    elseif cuVec.type==CUDA_DOUBLE_COMPLEX
        vec = zeros(cuVec.len, 1, 'like', 1i);
    elseif cuVec.type==CUDA_COMPLEX
        vec = zeros(cuVec.len, 1, 'like', single(1i));
    else
        vec = zeros(cuVec.len, 1);
    end
end

toplevel = nargout>2;
[vec, errCode] = cudaVecCopySubToHost(cuVec.len, cuVec, int32(1), ...
    vec, int32(1), int32(1));
