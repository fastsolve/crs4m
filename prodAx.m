function b = prodAx( A, x, b, nthreads, varargin)
%prodAx Compute b=A*x, assuming that size(A,2)==size(x,1).
%
%      b = prodAx( A, x [,b])
% Computes b=A*x in serial.
%
%      b = prodAx( A, x, b, nthreads)
% Computes b=A*x using nthreads OpenMP threads, where nthreads is an int32.
% In this mode, there are implicit barriers in the function.
%
%      b = prodAx( A, x, b, [])
% Computes b=A*x locally within each OpenMP thread, assuming the
% parallel region has already been initialized. The argument b has
% been preallocated. In this mode, there is no implicit barrier.
%
%      b = prodAx( A, x, b, nthreads, comm)
% also performs communication within the MPI communicator.
% It assumes that A is partitioned columnwise, and hence A*x is
% a partial sum on each process. An allreduce is performed on b
% within the communicator.
%
% Note that if A is partitioned rowwise, then b must be duplicated
% on all processes, and prodAx should be called without the communicator.
%
%      b = prodAx( A, x, b, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in b. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

%#codegen -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), ...
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0])}
%#codegen prodAx_ser -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen prodAx_ser1 -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen prodAx_mpi -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), MPI_Comm}
%#codegen prodAx_mpip -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), MPI_Comm, coder.typeof(0, [inf,inf])}
%#codegen prodAx_mpip1 -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), MPI_Comm, coder.typeof(0, [inf,inf]), int32(0)}

coder.inline('never');

if nargin==2;
    b = nullcopy(zeros(size(A,1),size(x,2)));
else
    if size(b,1)<size(A,1) || size(b,2)<size(x,2)
        m2c_error('prodAx:BufferTooSmal', 'Buffer space for output b is too small.');
    end
end
ismt = nargin>=4 && MACC_get_num_threads>1;

if nargin>3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~MACC_get_nested && ismt && nthreads(1)>1
        MACC_begin_master
        m2c_warn('prodAx:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        MACC_end_master
    end
    MACC_begin_parallel;
    MACC_clause_default('shared'); 
    MACC_clause_num_threads( int32(nthreads(1)));
    
    %% Compute b=A*x
    b = prodAx_internal( A, x, b);
    
    %% End parallel region
    MACC_end_parallel;
else
    %% Compute b=A*x
    % Decompose the work manually
    b = prodAx_internal( A, x, b);
end

if ~isempty(varargin)
    if isempty(nthreads) && ismt; MACC_barrier; end
    
    % Perform MPI allreduce
    MACC_begin_single
    s = int32(size(A,1)*size(x,2));
    b = allreduce(b, s, MPI_SUM, varargin{:});
    MACC_end_single
end

function b = prodAx_internal( A, x, b)
% Internal function for computing prodAx on each thread

coder.inline('never');

ncols = int32(size(A,2)); nrhs = int32(size(x,2));

% Decompose the work manually
[istart, iend] = get_local_chunk(size(A,1));
for i=istart:iend
    for k=1:nrhs
        t = 0;
        for j=1:ncols
            t = t + A(i,j) * x(j,k);
        end
        b(i,k) = t;
    end
end

function test %#ok<DEFNU>
%!test
%! if ~exist( ['prodAx.' mexext], 'file')
%!    m=200; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; A = rand(m,n); x=rand(size(A,2),2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = A*x;
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = prodAx( A, x);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( norm(b0-b1)/norm(b0)<=1.e-12);

%! b2 = zeros(size(A,1),2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>MACC_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; b2 = prodAx( A, x, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert( norm(b0-b2)/norm(b0)<=1.e-12);
%! end

%! fprintf(1, '\tTesting 2 threads with MPI call: ');
%! tic; b2 = prodAx( A, x, b2, int32(2), MPI_COMM_SELF);
%! fprintf(1, 'Done in %g seconds\n', toc);
%! assert( norm(b0-b2)/norm(b0)<=1.e-12);

%! b3 = zeros(size(A,1)+1,1);
%! fprintf(1, '\tTesting 2 threads with piggy-back: ');
%! tic; b3 = prodAx( A, x(:,1), b3, int32(2), MPI_COMM_SELF, 1, int32(1));
%! fprintf(1, 'Done in %g seconds\n', toc);
%! assert( norm(b0(:,1)-b3(1:end-1))/norm(b0(:,1))<=1.e-12);
%! assert( b3(end)==1);
