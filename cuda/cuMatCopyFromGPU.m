function [mat, errCode, toplevel] = cuMatCopyFromGPU(cuMat, mat, varargin)
% Copies a dense matrix from CUDA to MATLAB.
%
% mat = cuMatCopyFromGPU(cuMat) allocates mat in MATLAB
% and then copies cuMat on CUDA device to it.  The size to be copied
% is given by cuMat.dims. In code-generation mode, cuMat.type must be 
% a constant at compile time.
%
% [mat, errCode] = cuMatCopyFromGPU(cuMat, mat) copies the whole matrix
% from cuMat on CUDA device to mat in MATLAB.
%
% mat = cuMatCopyFromGPU(cuMat, mat, strm) copies in asynchronous mode,
% where strm is a CuStreamHandle. It is important for user not to access
% or deallocate mat before calling cuSynchronize(strm) or another
% synchronous operation on the stream.
%
% Note that the output argument mat is required, and it must be the
% same as the second input argument in the last two modes.
%
% SEE ALSO cuMatCopySubFromGPU, cuMatCopySubToGPU

toplevel = nargout>2;
if nargout<1
    % If no output argument, do nothing.
    if  isempty(coder.target) || m2c_debug
        m2c_error('cuMecCopyFromGPU:NoOutput', ...
            'The output argument vec is required and must be the same as the second input.\n');
    end
    return;
end

if nargin==1
    if cuMat.type==CU_SINGLE
        mat = zeros(cuMat.dims, 'single');
    elseif cuMat.type==CU_DOUBLE_COMPLEX
        mat = zeros(cuMat.dims, 'like', 1i);
    elseif cuMat.type==CU_COMPLEX
        mat = zeros(cuMat.dims, 'like', single(1i));
    else
        mat = zeros(cuMat.dims);
    end
end

[mat, errCode] = cuMatCopySubFromGPU(cuMat.dims(1), cuMat.dims(2), ...
    cuMat, mat, varargin{:});
