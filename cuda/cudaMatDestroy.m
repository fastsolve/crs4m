function [errCode, toplevel] = cudaMatDestroy(mat)
%Destroys a matrix on a CUDA device.
%
%  errCode = cudaMatDestroy(mat)
%  mat - Handle to matrix
%
%  SEE ALSO: cudaMatCreate

%#codegen -args {CudaMat}

coder.cinclude('mspack.h');

errCode = int32(-1);  %#ok<NASGU>
toplevel = nargout>1;

errCode = coder.ceval('cudaFree', CudaMat(mat, 'void *'));

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaFree returned error code %s\n', cudaGetErrorString(errCode));
end
end
