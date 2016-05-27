function prod = vecDot(x, y, prod, nthreads, varargin)
%Computes dot product x'*y for column vectors x and y.
%
% Syntax:
%      prod = vecDot(x, y)
%      prod = vecDot(x, y, [], n)
%      prod = vecDot(x, y, prod, n)
%      prod = vecDot(x, y, prod, n, 'acc')
%      prod = vecDot(x, y, prod, n, 'blas')
%
%      Both x and y must be column vectors.
%
% Description:
%
% prod = vecDot(x, y) computes x'*y in serial.
%
% prod = vecDot(x, y, [], n) computes the dot product by starting n threads.
% If n>OMP_NUM_THREADS, it will be reduced to OMP_NUM_THREADS automatically.
%
% prod = vecDot(x, y, prod, n) computes x'*y using m threads, assuming
% the function is called from one particular thread. The function must be
% called collectively by all the threads in the team, and prod is a vector
% shared by all the threads in the team. The length of prod must be >=
% nthreads, or memory error may occur.  The final output will be saved
% in prod(1). There is an implicit barrier at the beginning and the end
% of this function. This mode is not supported when vecDot is a top-level
% function during code generation.
%
% prod = vecDot(x, y, prod, n, 'acc') accelerates the computation using
% OpenACC within each thread. It is left to the copmiler to decide whether
% x and y  are already on GPU.
%
% prod = vecDot(x, y, prod, nthreads, 'blas') uses the BLAS routine to
% evaluate within each thread.
%
% Note: 'acc' and 'blas' must be known at compile time. For other
% strings, the solution is undefined.
%
% The compiled code only supports double-precision real numbers. The
% uncompiled code also works with complex numbers.
%
% See also vecMDot, VecTDot, VecNorm

%#codegen -args {m2c_realcol, m2c_realcol, m2c_realcol, m2c_int}
%#codegen vecDot_ser -args {m2c_realcol, m2c_realcol}
%#codegen vecDot_acc -args {m2c_realcol, m2c_realcol, m2c_realcol, m2c_int, m2c_string}

coder.inline('never');

toplevel = nargout>1;

if size(x,1)~=size(y,1) && (toplevel || m2c_debug)
    momp_begin_master
    m2c_error('vecDot:IncorrectSize', 'Dimensions for x and y must match.');
    momp_end_master
end

if nargin>3;
    nthreads = min(nthreads, momp_get_max_threads);
end

if nargin>3 && isempty(prod) && momp_get_num_threads>1
    %% Compute prod_local=x'*y
    prod = vecDot_partial(x, y, prod, true, varargin{:});
    momp_barrier
    
    momp_begin_master
    prod = accu_partsum(prod, nthreads);
    momp_end_master
    
    momp_barrier
elseif nargin>3 && nthreads>1
    %% Declare parallel region
    if ~momp_get_nested && momp_get_num_threads>1 && nthreads>1 && (toplevel || m2c_debug)
        momp_begin_master
        m2c_warn('vecDot:NestedParallel', ...
            'You are trying to use nested parallel regions, but nested parallelism is not enabled.\n');
        momp_end_master
    end
    
    prod = zeros(nthreads, 1);
    
    momp_begin_parallel; momp_clause_num_threads(nthreads);
    % In general, we should always use momp_clause_default('none')
    % and then explicitly define the shared varialble.
    momp_clause_default('none'); momp_clause_shared(x, y, prod);
    
    %% Compute prod_local=x'*y
    prod = vecDot_partial(x, y, prod, true, varargin{:});

    %% End parallel region
    momp_end_parallel;
    
    prod = accu_partsum(prod, nthreads);
else
    % Does not appear to be parallel. Use serial.
    prod = 0;
    prod = vecDot_partial(x, y, prod, false, varargin{:});
end

function prod = vecDot_partial(x, y, prod, is_mt, varargin)
if is_mt
    coder.inline('never');
    [istart, iend] = get_local_chunk(int32(size(x,1)));
    ind = momp_get_thread_num+1;
else
    ind = int32(1);
    istart = 1; iend = int32(size(x,1));
end

% Compute partial dot product with each therad
if ~isempty(varargin) && isequal(varargin{1}, 'blas')
    mode = 'blas';
elseif ~isempty(varargin) && isequal(varargin{1}, 'acc')
    mode = 'acc';
elseif isempty(varargin)
    mode = 'plane';
else
    % Should never reach here.
    mode = UNDEFINED;
end

if isequal(mode, 'blas') && m2c_blas
    % Call BLAS routine. Not yet implemented
    prod(ind) = coder.ceval('cblas_ddot', iend-istart+1, ...
        m2c_opaque_ptr_const(x, 'double *', istart-1), int32(1), ...
        m2c_opaque_ptr_const(y, 'double *', istart-1), int32(1));
else
    prod(ind) = 0.;
    if isequal(mode, 'acc')
        % Add pragma for acc kernel
        macc_begin_kernels;
    end
    for i=istart:iend
        prod(ind) = prod(ind) + x(i)' * y(i);
    end
    if isequal(mode, 'acc')
        % Add pragma for acc kernel
        macc_end_kernels;
    end
end

function prod = accu_partsum(prod, n)
for j=2:n
    prod(1) = prod(1) + prod(j);
end

function test %#ok<DEFNU>
%!shared x, y, b0
%!test
%! if ~exist(['vecDot.' mexext], 'file')
%!    m=200;
%! else
%!    m=1000000;
%! end
%! tic; x =rand(m,1); y=rand(m,1);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = x(:,1)'*y(:,1);
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = vecDot(x, y);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(abs(b0-b1)/abs(b0)<=1.e-6);
%!
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     b2 = zeros(nthreads,1);
%!     tic; b2 = vecDot(x, y, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2(1,1))/norm(b0)<=1.e-6);
%! end

%!test
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     b2 = zeros(nthreads,1);
%!     tic; b2 = vecDot(x, y, b2, nthreads, 'blas');
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2(1,1))/norm(b0)<=1.e-6);
%! end

%!test
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     b2 = zeros(nthreads,1);
%!     tic; b2 = vecDot(x, y, b2, nthreads, 'acc');
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2(1,1))/norm(b0)<=1.e-6);
%! end
