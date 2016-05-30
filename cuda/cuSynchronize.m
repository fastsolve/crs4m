function [errCode, toplevel] = cuSynchronize(stm)
%Blocks until stream has completed all operations.
%
%  errCode = cuSynchronize(stm)
%
%  SEE ALSO cuStreamCreate, cuStreamDestroy
%
% CUDA C interface:
%   cudaError_t cudaStreamSynchronize(cuStream_t handle)

%#codegen -args {CuStreamHandle}

coder.cinclude('mspack.h');

errCode = int32(-1);

if ~isempty(coder.target)
    errCode = coder.ceval('cudaStreamSynchronize', CuStreamHandle(stm));
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('CUDA:RuntimeError', 'cudaStreamSynchronize returned error code %s\n', cuGetErrorString(errCode));
    end
end
end
