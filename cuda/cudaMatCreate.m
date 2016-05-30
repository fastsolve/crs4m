function [mat, errCode, toplevel] = cudaMatCreate(m, n, type)
%Creates a matrix on a CUDA device.
%
%  [mat, errCode] = cudaMatCreate(m, n, [type]) creates a matrix of
%  size m-by-n on CUDA, where type is CUDA_DOUBLE, CUDA_SINGLE,
%  CUDA_COMPLEX, or CUDA_DOUBLE_COMPLEX. If not specified, then
%  type is CUDA_DOUBLE.
%
%  SEE ALSO: cudaMatDestroy, cudaMatCreate, cudaMatCopyFromHost, cudaMatCopyToHost

%#codegen -args {int32(0), int32(0), int32(0)}
%#codegen cudaMatCreate_2args -args {int32(0), int32(0)}

coder.cinclude('mspack.h');

errCode = int32(-1);
toplevel = nargout>2;

if nargin<3 || type == CUDA_DOUBLE
    [mat, errCode] = alloc(CUDA_DOUBLE, m, n, int32(8));
elseif type == CUDA_SINGLE
    [mat, errCode] = alloc(CUDA_SINGLE, m, n, int32(4));
elseif type == CUDA_COMPLEX
    [mat, errCode] = alloc(CUDA_COMPLEX, m, n, int32(8));
elseif type == CUDA_DOUBLE_COMPLEX
    [mat, errCode] = alloc(CUDA_DOUBLE_COMPLEX, m, n, int32(16));
elseif toplevel
    m2c_error('cudaMatCreate:WrongType', 'Unknow data type %d.\n', type)
    [mat, errCode] = alloc('double', m, n, int32(8));
end

if errCode && (toplevel || m2c_debug)
    m2c_error('CUDA:RuntimeError', 'cudaMalloc returned error code %s\n',...
        cudaGetErrorString(errCode));
end
end

function [mat, errCode] = alloc(type, m, n, sizepe)
coder.inline('always');

t_mat = coder.opaque('void *');
errCode = int32(0);  %#ok<NASGU>
errCode = coder.ceval('cudaMalloc', coder.wref(t_mat), m*n*sizepe);
mat = CudaMat(t_mat, type, m, n, 'wrap');
end
