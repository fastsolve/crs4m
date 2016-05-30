function [mat, errCode, toplevel] = cuMatCreate(m, n, type)
%Creates a matrix on a CUDA device.
%
%  [mat, errCode] = cuMatCreate(m, n, [type]) creates a matrix of
%  size m-by-n on CUDA, where type is CU_DOUBLE, CU_SINGLE,
%  CU_COMPLEX, or CU_DOUBLE_COMPLEX. If not specified, then
%  type is CU_DOUBLE.
%
%  SEE ALSO: cuMatDestroy, cuMatCreate, cuMatCopyToGPU, cuMatCopyFromGPU

%#codegen -args {int32(0), int32(0), int32(0)}
%#codegen cuMatCreate_2args -args {int32(0), int32(0)}

coder.cinclude('mspack.h');

toplevel = nargout>2;

if nargin<3; type = CU_DOUBLE; end

sizepe = cuGetSizePerElement(type);
t_mat = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_mat), m*n*sizepe);

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n',...
        cuGetErrorString(errCode));
end

mat = CuMat(t_mat, type, m, n, 'wrap');

end
