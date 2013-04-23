function b = crs_prodAtx( A, x, b, nthreads, varargin)
%crs_prodAtx Compute b=A'*x for a sparse matrix A in CRS format,
% assuming that number of columns of A is equal to size(x,1).
%
%      b = crs_prodAtx( A, x [, b])
% Computes b=A'*x in serial.
%
%      b = crs_prodAtx( A, x, b, nthreads)
% Computes b=A'*x using nthreads OpenMP threads, where nthreads is an int32.
% In this mode, there are implicit barriers in the function.
%
%      b = crs_prodAtx( A, x, b, [])
% Computes b=A'*x locally within each OpenMP thread, assuming the
% parallel region has already been initialized. The argument b has
% been preallocated. In this mode, there is no implicit barrier.
%
%      b = crs_prodAtx( A, x, b, nthreads, comm)
% also performs communication within the MPI communicator.
% It assumes that A is partitioned rowwise, and hence A'*x is
% a partial sum on each process. An allreduce is performed on b
% within the communicator.
%
% Note that if A is partitioned columnwise, then b must be duplicated
% on all processes, and crs_prodAtx should be called without the communicator.
%
%      b = crs_prodAtx( A, x, b, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in b. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

% See http://www.netlib.org/linalg/html_templates/node98.html

%#codegen -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), int32(1)}
%#codegen crs_prodAtx_ser -args {crs_matrix, coder.typeof(0, [inf,inf])}
%#codegen crs_prodAtx_ser1 -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen crs_prodAtx_mpi -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), int32(1), opaque_obj}
%#codegen crs_prodAtx_mpip -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), int32(1), opaque_obj, coder.typeof(0, [inf,1])}
%#codegen crs_prodAtx_mpip1 -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), int32(1), opaque_obj,coder.typeof(0, [inf,1]), int32(0)}

coder.inline('never');

DEBUG = true;

iend = A.nrows;
if nargin==2;
    b = nullcopy( zeros(A.ncols,size(x,2)));
else
    if size(b,1)<A.ncols || size(b,2)<size(x,2)
        msg_error('crs_prodAtx:BufferTooSmal', 'Buffer space for output b is too small.');
    end
end
ismt = nargin>3 && MACC_get_num_threads>1;

if nargin>3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~MACC_get_nested && ismt && nthreads>1
        MACC_begin_master
        msg_warn('prodAx:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        MACC_end_master
    end
    
    if DEBUG; T = MMPI_Wtime; end
    i = int32(0); j = int32(0); k = int32(0); n = int32(0);
    
    [b, i, j, k, n, iend] = MACC_begin_parallel( b, i, j, k, n, iend);
    MACC_clause_default('shared'); MACC_clause_private( i, j, k, n)
    MACC_clause_num_threads( int32(nthreads));

    MACC_begin_for
    for j=1:A.ncols
        for k=1:int32(size(x,2))
            b(j,k) = 0;
        end
    end
    MACC_end_for
    
    %% Compute b=A'*x
    MACC_begin_for
    for i=1:iend
        n = A.row_ptr(i+1) - 1;
        for j = A.row_ptr(i) : n
            for k=1:int32(size(x,2))
                MACC_atomic;
                b(A.col_ind(j), k) = b(A.col_ind(j), k) + A.val(j) * x(i, k);
            end
        end
    end
    MACC_end_for
    MACC_end_parallel
    
    if DEBUG
        T = MMPI_Wtime-T;
        msg_printf( 'csr_prodAx took %g seconds\n', T);
    end
else
    if DEBUG; T = MMPI_Wtime; end
    %% Compute b=A'*x
    [jstart, jend] = get_local_chunk(A.ncols);
    for j=jstart:jend
        for k=1:int32(size(x,2))
            b(j,k) = 0;
        end
    end
    if nargin>3; MACC_barrier; end

    % Decompose the work manually
    [istart, iend] = get_local_chunk(A.nrows);
    for i=istart:iend
        n = A.row_ptr(i+1) - 1;
        for j = A.row_ptr(i) : n
            for k=1:int32(size(x,2))
                if nargin>3; MACC_atomic; end
                b(A.col_ind(j), k) = b(A.col_ind(j), k) + A.val(j) * x(i, k);
            end
        end
    end
    if DEBUG
        T = MMPI_Wtime-T;
        if nargin>3; MACC_begin_master; end
        msg_printf( 'csr_prodAx took %g seconds\n', T);
        if nargin>3; MACC_end_master; end
    end
end

if ~isempty(varargin)
    if isempty(nthreads) && ismt; MACC_barrier; end
    
    % Perform MPI allreduce
    MACC_begin_single
    s = int32(A.ncols*size(x,2));
    b = allreduce(b, s, MPI_SUM, varargin{:});
    MACC_end_single
end

function test %#ok<DEFNU>

%!test
%! sp = sprand(1000,200,0.5); 
%! A = crs_matrix( sp); x=rand(size(sp,1),2);
%! b0 = sp'*x;
%! b1 = crs_prodAtx( A, x);
%! assert( norm(b0-b1)<=1.e-12);
%! b2 = zeros(size(sp,2),2);
%! b2 = crs_prodAtx( A, x, b2, int32(2));
%! assert( norm(b0-b2)<=1.e-12);
%
%! if ~MMPI_Initialized; MMPI_Init; end
%! nprocs = double(MMPI_Comm_size(MPI_COMM_WORLD));
%! b2 = crs_prodAtx( A, x, b2, MACC_get_max_threads, MPI_COMM_WORLD);
%! assert( nprocs>1 || norm(b0-b2)<=1.e-12);

%! b3 = zeros(size(sp,2)+1,1);
%! b3 = crs_prodAtx( A, x(:,1), b3, MACC_get_max_threads, MPI_COMM_WORLD, 1, int32(1));
%! assert( nprocs>1 || norm(b0(:,1)-b3(1:end-1))<=1.e-12);
%! assert( b3(end)==nprocs);
