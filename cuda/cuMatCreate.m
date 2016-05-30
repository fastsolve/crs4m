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

errCode = int32(-1);
toplevel = nargout>2;

if nargin<3 || type == CU_DOUBLE
    [mat, errCode] = alloc(CU_DOUBLE, m, n, int32(8));
elseif type == CU_SINGLE
    [mat, errCode] = alloc(CU_SINGLE, m, n, int32(4));
elseif type == CU_COMPLEX
    [mat, errCode] = alloc(CU_COMPLEX, m, n, int32(8));
elseif type == CU_DOUBLE_COMPLEX
    [mat, errCode] = alloc(CU_DOUBLE_COMPLEX, m, n, int32(16));
elseif toplevel
    m2c_error('cuMatCreate:WrongType', 'Unknow data type %d.\n', type)
    [mat, errCode] = alloc('double', m, n, int32(8));
end

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n',...
        cuGetErrorString(errCode));
end
end

function [mat, errCode] = alloc(type, m, n, sizepe)
coder.inline('always');

t_mat = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_mat), m*n*sizepe);
mat = CuMat(t_mat, type, m, n, 'wrap');
end
