function [vec, errCode, toplevel] = cuVecCreate(n, type)
%Creates a vector on a CUDA device.
%
%  [vec, errCode] = cuVecCreate(n, [type]) creates a vector of length n 
%  on CUDA, where type is CU_DOUBLE, CU_SINGLE, CU_COMPLEX, or 
%  CU_DOUBLE_COMPLEX. If not specified, then type is CU_DOUBLE.
%
%  SEE ALSO: cuVecDestroy, cuMatCreate, cuVecCopyToGPU, cuVecCopyFromGPU

%#codegen -args {int32(0), int32(0)}
%#codegen cuVecCreate_1arg -args {int32(0)}

coder.cinclude('mspack.h');

errCode = int32(-1);
toplevel = nargout>2;

if nargin<2 || type == CU_DOUBLE
    [vec, errCode] = alloc(CU_DOUBLE, n, int32(8));
elseif type == CU_SINGLE
    [vec, errCode] = alloc(type, n, int32(4));
elseif type == CU_COMPLEX
    [vec, errCode] = alloc(type, n, int32(8));
elseif type == CU_DOUBLE_COMPLEX
    [vec, errCode] = alloc(type, n, int32(16));
elseif toplevel
    m2c_error('cuVecCreate:WrongType', 'Unknow data type %d.\n', type)
    [vec, errCode] = alloc('double', n, int32(8));
end

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n', cuGetErrorString(errCode));
end
end

function [vec, errCode] = alloc(type, n, sizepe)
coder.inline('always');

t_vec = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_vec), n*sizepe);
vec = CuVec(t_vec, type, n, 'wrap');
end
