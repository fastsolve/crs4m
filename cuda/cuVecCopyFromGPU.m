function [vec, errCode, toplevel] = cuVecCopyFromGPU(cuVec, vec)
% Copies a vector from CUDA to MATLAB.
%
% [vec, errCode] = cuVecCopyFromGPU(cuVec, vec) copies the whole array
% from cuVec on CUDA device to vec in MATLAB. The length to be copied
% is equal to cuVec.len.
%
% [vec, errCode] = cuVecCopyFromGPU(cuVec) allocates vec in MATLAB
% and then copies cuVec on CUDA device to it. In code-generation mode,
% cuVec.type must be a constant at compile time.
%
% SEE ALSO cuVecCopySubFromGPU, cuVecCopySubToGPU

if nargin==1
    if cuVec.type==CU_SINGLE
        vec = zeros(cuVec.len, 1, 'single');
    elseif cuVec.type==CU_DOUBLE_COMPLEX
        vec = zeros(cuVec.len, 1, 'like', 1i);
    elseif cuVec.type==CU_COMPLEX
        vec = zeros(cuVec.len, 1, 'like', single(1i));
    else
        vec = zeros(cuVec.len, 1);
    end
end

toplevel = nargout>2;
[vec, errCode] = cuVecCopySubFromGPU(cuVec.len, cuVec, int32(1), ...
    vec, int32(1), int32(1));
