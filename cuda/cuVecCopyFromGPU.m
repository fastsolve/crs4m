function [vec, errCode, toplevel] = cuVecCopyFromGPU(cuVec, vec, varargin)
% Copies a vector from CUDA to MATLAB.
%
% [vec, errCode] = cuVecCopyFromGPU(cuVec) allocates vec in MATLAB
% and then copies cuVec on CUDA device to it. The length to be copied
% is equal to cuVec.len. In code-generation mode, cuVec.type must be
% a constant at compile time.
%
% vec = cuVecCopyFromGPU(cuVec, vec) copies the whole array
% from cuVec on CUDA device to vec in MATLAB.
%
% vec = cuVecCopyFromGPU(cuVec, vec, strm) copies in asynchronous mode,
% where strm is a CuStreamHandle. It is important for user not to access
% or deallocate vec before calling cuSynchronize(strm) or another
% synchronous operation on the stream.
%
% Note that the output argument vec is required, and it must be the
% same as the second input argument in the last two modes.
%
% SEE ALSO cuVecCopySubFromGPU, cuVecCopySubToGPU

toplevel = nargout>2;
if nargout<1
    % If no output argument, do nothing.
    if  isempty(coder.target) || m2c_debug
        m2c_error('cuVecCopyFromGPU:NoOutput', ...
            'The output argument vec is required and must be the same as the second input.\n');
    end
    return;
end

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

[vec, errCode] = cuVecCopySubFromGPU(cuVec.len, cuVec, int32(1), ...
    vec, int32(1), int32(1), varargin{:});
