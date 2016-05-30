function [vec, errCode, toplevel] = cuVecCreate(n, type)
%Creates a vector on a CUDA device.
%
%  [vec, errCode] = cuVecCreate(n, [type]) creates a vector of length n
%  on CUDA, where type is MSP_INT*, MSP_UINT*, MSP_DOUBLE, MSP_SINGLE,
%  MSP_COMPLEX, or MSP_DOUBLE_COMPLEX, MSP_INT*, or MSP_UINT*. If not
%  specified, then type is MSP_DOUBLE.
%
%  SEE ALSO: cuVecDestroy, cuVecCopyToGPU, cuVecCopyFromGPU

%#codegen -args {int32(0), int32(0)}
%#codegen cuVecCreate_1arg -args {int32(0)}

coder.cinclude('mspack.h');

toplevel = nargout>2;

if nargin<2; type = MSP_DOUBLE; end

sizepe = mspGetSizePerElement(type);

t_vec = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_vec), n*sizepe);

vec = CuVec(t_vec, type, n, 'wrap');

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n', cuGetErrorString(errCode));
end
end
