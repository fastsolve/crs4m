function [errCode, toplevel] = cudaVecDestroy(vec)
%Destroys a vector on a CUDA device.
%
%  errCode = cudaVecDestroy(vec)
%  vec - Handle to vector
%
%  SEE ALSO: cudaVecCreate

%#codegen -args {CudaVec}

coder.cinclude('mspack.h');

errCode = int32(-1);  %#ok<NASGU>
toplevel = nargout>1;

errCode = coder.ceval('cudaFree', CudaVec(vec, 'void *'));

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaFree returned error code %s\n', cudaGetErrorString(errCode));
end
end
