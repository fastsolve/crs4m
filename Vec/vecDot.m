function [prod, buf, toplevel] = vecDot(x, y, buf, varargin)
%Computes dot product x'*y for column vectors x and y.
%
% Syntax:
%  Serial modes:
%      prod = vecDot(x, y)
%      [prod,buf] = vecDot(x, y, buf)
%
%      Both x and y must be column vectors.
%
%  OpenMP modes:
%      prod = vecDot(x, y, [], n)
%      [prod,buf] = vecDot(x, y, buf, n)
%
%      Both x and y must be column vectors.
%
%  CUDA modes:
%      vecDot(x, y, hdl, 'cuda', cublasHandle)
%      prod = vecDot(x, y, [], 'cuda', cublasHandle, ['sync'])
%
%      Both x and y must be CuVec handles.
%
% Description:
%
% prod = vecDot(x, y) computes x'*y in serial and returns a scalar.
%
% [prod,buf] = vecDot(x, y, buf) computes x'*y and saves the output into prod.
%
% prod = vecDot(x, y, [], n) computes the dot product by starting n threads.
% If n>OMP_NUM_THREADS, it will be reduced to OMP_NUM_THREADS automatically.
%
% [prod,buf] = vecDot(x, y, buf, n) computes x'*y using m threads, assuming
% the function is called from one particular thread. In this mode, the
% function must be called collectively by all the threads in the team,
% and buf is a vector shared by all the threads in the team. The length
% of buf must be >= nthreads, and it is important that buf is passed both
% as input and output with the same name. The final output will be saved
% in prod. There is an implicit barrier at the beginning and the end
% of this function. This mode is not supported when vecDot is a top-level
% function during code generation.
%
% If blas is enabled, it calls BLAS routines automatically for vectors
% with more than 1000 entries in the CPU mode.
%
% vecDot(x, y, prod, 'cuda', cublasHandle) uses CUDA BLAS for acceleration,
% assuming x, y, and prod are all CuVec handles. It uses the
% asynchronous mode of CUDA BLAS. This should be the default mode on CUDA.
% The output values are meaningless.
%
% prod = vecDot(x, y, [], 'cuda', cublasHandle, ['sync']) uses the
% synchronous mode and returns a scalar to the host.
% The last argument is required when calling the MEX version of vecDot.
%
% The compiled code only supports double-precision real numbers. The
% uncompiled code also works with single-precision real numbers and
% with complex numbers for code generation.
%
% The function checks errors when m2c_debug is on or when the function is
% compiled into a MEX function.
%
% See also vecMDot, VecTDot, VecNorm

%#codegen -args {m2c_vec, m2c_vec, m2c_vec}
%#codegen vecDot_ser -args {m2c_vec, m2c_vec}
%#codegen vecDot_omp -args {m2c_vec, m2c_vec, m2c_mat, m2c_int}
%#codegen vecDot_cublas -args {CuVec, CuVec, CuVec, ...
%#codegen         m2c_string, CuBlasHandle}
%#codegen vecDot_cublas_sync -args {CuVec, CuVec, m2c_mat, ...
%#codegen         m2c_string, CuBlasHandle, m2c_string}

if nargin==4;
    coder.inline('never'); % Never inline in OpenMP mode
end
narginchk(2, 6);

toplevel = nargout>2 || ~isempty(varargin) && islogical(varargin{1});
if (isnumeric(x) && (size(x,1)~=size(y,1) || ~isequal(class(x), class(y))) ||...
        isstruct(x) && (x.len~=y.len || x.type ~= y.type)) ...
        && (toplevel || m2c_debug)
    m2c_error('vecDot:IncorrectSize', 'Vectors x and y must have the same type and size.\n');
end

if nargin==4 && ~isempty(buf) && length(buf)<varargin{1} && (toplevel || m2c_debug)
    m2c_error('vecDot:BufferTooSmall', 'Array buf must hold one entry per threads.\n');
end

if nargin>=5
    if ischar(varargin{1}) && ~isequal(varargin{1}, 'cuda') && (toplevel || m2c_debug)
        m2c_warn('vecDot:WrongInput', 'Wrong mode name. cuda is assumed.\n');
    elseif toplevel && x.type ~= MCU_DOUBLE
        m2c_error('vecDot:WrongInput', 'When using mex functions, vecDot only supports double.\n');
    end
    if ~m2c_cuda && (toplevel || m2c_debug)
        m2c_error('vecDot:WrongInput', 'CUDA was not enabled during compilation.\n');
    end
    
    [prod, errCode] = vecDot_cuda_kernel(x, y, buf, varargin{2}, toplevel, varargin{3:end});
    
    if errCode && (toplevel || m2c_debug)
        if errCode<0
            m2c_error('vecDot:WrongPointerMode', 'The given cuBLAS handle has incorrect pointer mode.\n.');
        else
            m2c_error('cuBLAS:RuntimeError', 'cuBLAS returned an error code %s\n.', ...
                cuBlasGetErrorString(errCode));
        end
    end
    return;
end

if nargin>3;
    nthreads = min(int32(varargin{1}), ompGetMaxThreads);
end

if nargin>3 && isempty(buf) && ompGetNumThreads>1
    %% Compute prod_local = x'*y
    buf = vecDot_partial(x, y, buf, nthreads, int32(size(x,1)));
    
    OMP_barrier; OMP_begin_master
    [prod,buf] = accu_partsum(buf, nthreads);
    OMP_end_master; OMP_barrier
elseif nargin>3 && nthreads>1
    %% Declare parallel region
    if ~ompGetNested && ompGetNumThreads>1 && nthreads>1 && (toplevel || m2c_debug)
        OMP_begin_master
        m2c_warn('vecDot:NestedParallel', ...
            'You are trying to use nested parallel regions, but nested parallelism is not enabled.\n');
        OMP_end_master
    end
    
    buf = zeros(nthreads, 1);
    m = int32(size(x,1));
    
    %% Compute prod_local=x'*y
    m2c_rref(buf); OMP_begin_parallel(nthreads);
    buf = vecDot_partial(x, y, buf, nthreads, m);
    OMP_end_parallel;
    
    [prod,buf] = accu_partsum(buf, nthreads);
else
    % Does not appear to be parallel. Use serial.
    prod = 0;
    prod = vecDot_partial(x, y, prod, int32(1), int32(size(x,1)));
end

function prod = vecDot_partial(x, y, prod, nthreads, n)
% Never inline in OpenMP mode to avoid private variables become shared
coder.inline('never');

if nthreads>1
    [istart, iend, threadId] = OMP_local_chunk(n, nthreads);
else
    istart = int32(1); iend = n; threadId = int32(0);
end

LARGE = 1000;

if iend>=istart+LARGE && m2c_blas  % Use BLAS only for large-enough vectors
    % Call BLAS routine.
    if isreal(x)
        if isa(x, 'double')
            func = 'cblas_ddot'; type = 'double';
        elseif isa(x, 'single')
            func = 'cblas_sdot'; type = 'single';
        end
        prod(threadId+1) = coder.ceval(func, iend-istart+1, ...
            m2c_opaque_ptr_const(x, [type '*'], istart-1), int32(1), ...
            m2c_opaque_ptr_const(y, [type '*'], istart-1), int32(1));
    else
        if isa(x, 'double')
            func = 'cblas_zdotc_sub'; type = 'creal64_T';
        elseif isa(x, 'single')
            func = 'cblas_cdotc_sub'; type = 'creal32_T';
        end
        coder.ceval(func, iend-istart+1, ...
            m2c_opaque_ptr_const(x, [type '*'], istart-1), int32(1), ...
            m2c_opaque_ptr_const(y, [type '*'], istart-1), int32(1), ...
            m2c_opaque_ptr(prod, [type '*'], threadId));
    end
else
    prod(threadId+1) = 0.;
    for i=istart:iend
        prod(threadId+1) = prod(threadId+1) + x(i)' * y(i);
    end
end

function [prod, buf] = accu_partsum(buf, n)
coder.inline('always');
for j=2:n
    buf(1) = buf(1) + buf(j);
end
prod = buf(1);

function [output, errCode] = vecDot_cuda_kernel(x, y, hdl, cublasHdl, toplevel, varargin)
coder.inline('always');

if toplevel || x.type==MCU_DOUBLE
    func = 'cublasDdot'; type = 'double'; mzero = 0;
elseif x.type==MCU_SINGLE
    func = 'cublasSdot'; type = 'single'; mzero = single(0);
elseif x.type==MCU_DOUBLE_COMPLEX
    func = 'cublasZdotc'; type = 'cuDoubleComplex'; mzero = complex(0);
elseif  x.type==MCU_COMPLEX
    func = 'cublasCdotc'; type = 'cuComplex'; mzero = complex(single(0));
else
    % Undefined behavior. x.type must be a constant string.
    % This causes a compilation error.
    return;
end

output = mzero;
n = x.len;
errCode = int32(0); %#ok<NASGU>
if toplevel || m2c_debug
    [mode, errCode] = cuBlasGetPointerMode(CuBlasHandle(cublasHdl));
    if errCode; return; end
    
    if ((isempty(hdl) || ~isempty(varargin)) && mode ~= CUBLAS_POINTER_MODE_HOST) || ...
            (~isempty(hdl) && isempty(varargin) && mode ~= CUBLAS_POINTER_MODE_DEVICE)
        errCode = int32(-1);
        return;
    end
end

if isempty(hdl) || ~isempty(varargin)
    % Calls the synchronous version of cuBLAS
    errCode = coder.ceval(func, CuBlasHandle(cublasHdl), n, ...
        CuVec(x, [type '*']), int32(1), ...
        CuVec(y, [type '*']), int32(1), coder.wref(output));
else
    % Calls the asynchronous version of cuBLAS
    errCode = coder.ceval(func, CuBlasHandle(cublasHdl), n, ...
        CuVec(x, [type '*']), int32(1), ...
        CuVec(y, [type '*']), int32(1), ...
        CuVec(hdl, [type '*']));
end

function test %#ok<DEFNU>
%!shared x, y, b0, m
%!test
%! if ~exist(['vecDot.' mexext], 'file')
%!    m=int32(200);
%! else
%!    m=int32(1000000);
%! end
%! tic; x =rand(m,1); y=rand(m,1);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = x(:,1)'*y(:,1);
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; prod = vecDot(x, y);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(abs(b0-prod)/abs(b0)<=1.e-6);
%!
%! fprintf(1, 'Running  with blas mode (if blas is enabled):\n');
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>ompGetMaxThreads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     b2 = zeros(nthreads,1);
%!     tic; [prod, b2] = vecDot(x, y, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-prod)/abs(b0)<=1.e-6);
%! end

%!test
%! fprintf(1, 'Running with cuda in synchronous mode:\n');
%! cuda = cuBlasCreate;
%! cuda_x = cuVecCopyToGPU(x);
%! cuda_y = cuVecCopyToGPU(y);
%! tic; prod = vecDot(cuda_x, cuda_y, [], 'cuda', cuda, 'sync');
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! cuVecDestroy(cuda_x);
%! cuVecDestroy(cuda_y);
%! cuBlasDestroy(cuda);
%! assert(norm(b0-prod)/norm(b0)<=1.e-6);

%!test
%! fprintf(1, 'Running with cuda in asynchronous mode:\n');
%! prod = 0;
%! cuda = cuBlasCreate; cuBlasSetPointerMode(cuda, CUBLAS_POINTER_MODE_DEVICE);
%! cuda_x = cuVecCopyToGPU(x);
%! cuda_y = cuVecCopyToGPU(y);
%! cuda_prod = cuVecCreate(int32(1));
%! tic; vecDot(cuda_x, cuda_y, cuda_prod, 'cuda', cuda);
%! prod = cuVecCopyFromGPU(cuda_prod);
%! fprintf(1, 'Done in %g seconds (including copying data from GPU)\n', toc);
%! cuVecDestroy(cuda_x);
%! cuVecDestroy(cuda_y);
%! cuVecDestroy(cuda_prod);
%! cuBlasDestroy(cuda);
%! assert(norm(b0-prod)/abs(b0)<=1.e-6);
