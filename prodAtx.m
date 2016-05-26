function b = prodAtx(A, x, b, nthreads, varargin)
%prodAtx Compute b=A'*x, assuming that size(A,1)==size(x,1).
%
%      b = prodAtx(A, x [, b])
% Computes b=A'*x in serial.
%
%      b = prodAtx(A, x, b, nthreads)
% Computes b=A'*x using nthreads OpenMP threads, where nthreads is an int32.
% In this mode, there are implicit barriers in the function.
%
%      b = prodAtx(A, x, b, [])
% Computes b=A'*x locally within each OpenMP thread, assuming the
% parallel region has already been initialized. The argument b has
% been preallocated. In this mode, there is no implicit barrier.
%
%      b = prodAtx(A, x, b, nthreads, comm)
% also performs communication within the MPI communicator.
% It assumes that A is partitioned rowwise, and hence A'*x is
% a partial sum on each process. An allreduce is performed on b
% within the communicator.
%
% Note that if A is partitioned columnwise, then b must be duplicated
% on all processes, and prodAtx should be called without the communicator.
%
%      b = prodAtx(A, x, b, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in b. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

%#codegen -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), ...
%#codegen coder.typeof(0, [inf,inf]), int32(1)}
%#codegen prodAtx_ser -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen prodAtx_ser1 -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
% %#codegen prodAtx_mpi -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
% %#codegen coder.typeof(0, [inf,inf]), coder.typeof(int32(1), [1,1], [1,0]), MPI_Comm}
% %#codegen prodAtx_mpip -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,1]),
% %#codegen coder.typeof(0, [inf,inf]), coder.typeof(int32(1), [1,1], [1,0]), MPI_Comm, coder.typeof(0, [inf,inf])}
% %#codegen prodAtx_mpip1 -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
% %#codegen coder.typeof(0, [inf,inf]), coder.typeof(int32(1), [1,1], [1,0]), MPI_Comm, coder.typeof(0, [inf,inf]), int32(0)}

coder.inline('never');

if nargin==2;
    b = nullcopy(zeros(size(A,2),size(x,2)));
else
    if size(b,1)<size(A,2) || size(b,2)<size(x,2)
        m2c_error('prodAtx:BufferTooSmal', 'Buffer space for output b is too small.');
    end
end
ismt = nargin>3 && momp_get_num_threads>1;

if nargin>3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~momp_get_nested && ismt && nthreads(1)>1
        momp_begin_master
        m2c_warn('prodAtx:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        momp_end_master
    end
    
    momp_begin_parallel; momp_clause_default('shared');
    momp_clause_num_threads(int32(nthreads(1)));
    
    %% Compute b=A'*x
    b = prodAtx_internal(A, x, b);
    
    momp_end_parallel
else
    %% Compute b=A'*x
    b = prodAtx_internal(A, x, b);
end

if ~isempty(varargin)
    if isempty(nthreads) && ismt; momp_barrier; end
    
    % Perform MPI allreduce
    momp_begin_single
    s = int32(size(A,2)*size(x,2));
    b = allreduce(b, s, MPI_SUM, varargin{:});
    momp_end_single
end

function b = prodAtx_internal(A, x, b)
% Internal function for computing prodAtx on each thread

coder.inline('never');

nrows = int32(size(A,1)); nrhs = int32(size(x,2));

% Decompose the work manually
[jstart, jend] = get_local_chunk(size(A,2));
for j=jstart:jend
    for k=1:nrhs
        t = 0;
        for i=1:nrows
            t = t + A(i,j) * x(i,k);
        end
        b(j,k) = t;
    end
end

function test %#ok<DEFNU>
%!test
%! if ~exist(['prodAtx.' mexext], 'file')
%!    m=200; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; A = rand(m,n); x=rand(m,2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = A'*x;
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = prodAtx(A, x);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(norm(b0-b1)/norm(b0)<=1.e-12);

%! b2 = zeros(n,2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; b2 = prodAtx(A, x, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2)/norm(b0)<=1.e-12);
%! end

% %! fprintf(1, '\tTesting 2 threads with MPI call: ');
% %! tic; b2 = prodAtx(A, x, b2, int32(2), MPI_COMM_SELF);
% %! fprintf(1, 'Done in %g seconds\n', toc);
% %! assert(norm(b0-b2)/norm(b0)<=1.e-12);
% 
% %! b3 = zeros(n+1,1);
% %! fprintf(1, '\tTesting 2 threads with piggy-back: ');
% %! tic; b3 = prodAtx(A, x(:,1), b3, int32(2), MPI_COMM_WORLD, 1, int32(1));
% %! fprintf(1, 'Done in %g seconds\n', toc);
% %! assert(norm(b0(:,1)-b3(1:end-1))/norm(b0(:,1))<=1.e-12);
% %! assert(b3(end)==1);
