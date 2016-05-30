function [vec, errCode, toplevel] = cuVecCreate(n, type)
%Creates a vector on a CUDA device.
%
%  [vec, errCode] = cuVecCreate(n, [type]) creates a vector of length n
%  on CUDA, where type is CU_INT*, CU_UINT*, CU_DOUBLE, CU_SINGLE,
%  CU_COMPLEX, or CU_DOUBLE_COMPLEX. If not specified, then type is CU_DOUBLE.
%
%  SEE ALSO: cuVecDestroy, cuMatCreate, cuVecCopyToGPU, cuVecCopyFromGPU

%#codegen -args {int32(0), int32(0)}
%#codegen cuVecCreate_1arg -args {int32(0)}

coder.cinclude('mspack.h');

toplevel = nargout>2;

if nargin<2; type = CU_DOUBLE; end

sizepe = cuGetSizePerElement(type);

t_vec = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_vec), n*sizepe);
vec = CuVec(t_vec, type, n, 'wrap');

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n', cuGetErrorString(errCode));
end
end
