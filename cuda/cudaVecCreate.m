function [vec, errCode, toplevel] = cudaVecCreate(n, type)
%Creates a vector on a CUDA device.
%
%  [vec, errCode] = cudaVecCreate(n, [type]) creates a vector of length n 
%  on CUDA, where type is CUDA_DOUBLE, CUDA_SINGLE, CUDA_COMPLEX, or 
%  CUDA_DOUBLE_COMPLEX. If not specified, then type is CUDA_DOUBLE.
%
%  SEE ALSO: cudaVecDestroy, cudaMatCreate, cudaVecCopyFromHost, cudaVecCopyToHost

%#codegen -args {int32(0), int32(0)}
%#codegen cudaVecCreate_1arg -args {int32(0)}

coder.cinclude('mspack.h');

errCode = int32(-1);
toplevel = nargout>2;

if nargin<2 || type == CUDA_DOUBLE
    [vec, errCode] = alloc(CUDA_DOUBLE, n, int32(8));
elseif type == CUDA_SINGLE
    [vec, errCode] = alloc(type, n, int32(4));
elseif type == CUDA_COMPLEX
    [vec, errCode] = alloc(type, n, int32(8));
elseif type == CUDA_DOUBLE_COMPLEX
    [vec, errCode] = alloc(type, n, int32(16));
elseif toplevel
    m2c_error('cudaVecCreate:WrongType', 'Unknow data type %d.\n', type)
    [vec, errCode] = alloc('double', n, int32(8));
end

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n', cudaGetErrorString(errCode));
end
end

function [vec, errCode] = alloc(type, n, sizepe)
coder.inline('always');

t_vec = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_vec), n*sizepe);
vec = CudaVec(t_vec, type, n, 'wrap');
end
