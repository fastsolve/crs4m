function [errCode, toplevel] = cuVecDestroy(vec)
%Destroys a vector on a CUDA device.
%
%  errCode = cuVecDestroy(vec)
%  vec - Handle to vector
%
%  SEE ALSO: cuVecCreate

%#codegen -args {CuVec}

coder.cinclude('mspack.h');

errCode = int32(-1);  %#ok<NASGU>
toplevel = nargout>1;

errCode = coder.ceval('cudaFree', CuVec(vec, 'void *'));

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaFree returned error code %s\n', cuGetErrorString(errCode));
end
end
