function b = crs_prodAx(A, x, b, nthreads, varargin)
%crs_prodAx Compute b=A*x for a sparse matrix A in CRS format,
% assuming that number of columns of A is equal to size(x,1).
%
%      b = crs_prodAx( A, x [,b])
% Computes b=A*x in serial.
%
%      b = crs_prodAx( A, x, b, nthreads)
% Computes b=A*x Testing nthreads OpenMP threads, where nthreads is an int32.
% In this mode, there are implicit barriers in the function.
%
%      b = crs_prodAx( A, x, b, [])
% Computes b=A*x locally within each OpenMP thread, assuming the
% parallel region has already been initialized. The argument b has
% been preallocated. In this mode, there is no implicit barrier.
%
%      b = crs_prodAx( A, x, b, nthreads, comm)
% also performs communication within the MPI communicator.
% It assumes that A is partitioned columnwise, and hence A*x is
% a partial sum on each process. An allreduce is performed on b
% within the communicator.
%
% Note that if A is partitioned rowwise, then b must be duplicated
% on all processes, and prodAx should be called without the communicator.
%
%      b = crs_prodAx( A, x, b, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in b. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

% See http://www.netlib.org/linalg/html_templates/node98.html

%#codegen -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), int32(1)}
%#codegen crs_prodAx_ser -args {crs_matrix, coder.typeof(0, [inf,inf])}
%#codegen crs_prodAx_ser1 -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen crs_prodAx_mpi -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), int32(1), opaque_obj}
%#codegen crs_prodAx_mpip -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), int32(1), opaque_obj, coder.typeof(0, [inf,1])}
%#codegen crs_prodAx_mpip1 -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), int32(1), opaque_obj, coder.typeof(0, [inf,1]), int32(0)}

coder.inline('never');

DEBUG = true;

if nargin==2;
    b = nullcopy(zeros(A.nrows,size(x,2)));
else
    if size(b,1)<A.nrows || size(b,2)<size(x,2)
        msg_error('crs_prodAx:BufferTooSmal', 'Buffer space for output b is too small.');
    end
end
ismt = nargin>=4 && MACC_get_num_threads>1;

if nargin>3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~MACC_get_nested && ismt && nthreads>1
        MACC_begin_master
        msg_warn('crs_prodAx:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        MACC_end_master
    end
    if DEBUG; T = MACC_get_wtime; end
    istart = int32(0); iend = int32(0);

    [b, istart, iend] = MACC_begin_parallel( b, istart, iend);
    MACC_clause_default('shared'); MACC_clause_private( istart, iend)
    MACC_clause_num_threads( int32(nthreads));
    
    %% Compute b=A*x
    [istart, iend] = get_local_chunk(A.nrows);
    b = crs_prodAx_internal( A.row_ptr, A.col_ind, A.val, x, istart, iend, b);
    
    %% End parallel region
    MACC_end_parallel;
    if DEBUG
        T = MACC_get_wtime-T;
        msg_printf( 'crs_prodAx took %g seconds\n', T);
    end
else
    if DEBUG; T = MACC_get_wtime; end

    %% Compute b=A*x
    % Decompose the work manually
    [istart, iend] = get_local_chunk(A.nrows);
    b = crs_prodAx_internal( A.row_ptr, A.col_ind, A.val, x, istart, iend, b);
    
    if DEBUG
        T = MACC_get_wtime-T;
        if nargin>3; MACC_begin_master; end
        msg_printf( 'crs_prodAx took %g seconds\n', T);
        if nargin>3; MACC_end_master; end
    end

end

if ~isempty(varargin)
    if isempty(nthreads) && ismt; MACC_barrier; end
    
    % Perform MPI allreduce
    MACC_begin_single
    s = int32(A.nrows*size(x,2));
    b = allreduce(b, s, MPI_SUM, varargin{:});
    MACC_end_single
end

function b = crs_prodAx_internal( row_ptr, col_ind, val, x, istart, iend, b)

if isempty( coder.target)
    for i=istart:iend
        for k=1:int32(size(x,2))
            t = 0.0;
            for j=row_ptr(i):row_ptr(i+1)-1
                t = t + val(j)*x(col_ind(j),k);
            end
            b(i,k) = t;
        end
    end
else
    coder.inline('never');
    coder.cinclude('spalab_kernel.h');
    
    if isa( val, 'double'); func = 'SPL_ddot'; 
    else func = 'SPL_sdot'; end

    for i=istart:iend
        for k=1:int32(size(x,2))
            b(i,k) = coder.ceval( func, coder.rref(val( row_ptr(i))), ...
                coder.rref( col_ind(row_ptr(i))), coder.rref(x(1,k)), ...
                row_ptr(i+1) - row_ptr(i));
        end
    end
end

function test %#ok<DEFNU>
%!test
%! if ~exist( ['crs_prodAx.' mexext], 'file')
%!    m=100; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; sp = sprand(m,n,0.5); x=rand(size(sp,2),2);
%! [is,js,vs] = find(sp); 
%! fprintf(1, 'Generated random matrix in %g seconds\n', toc);
%! tic; b0 = sp*x;
%! fprintf(1, 'Computed reference solution in %g seconds\n', toc);
%! tic; A = crs_matrix(int32(is), int32(js), vs, int32(size(sp,1)), int32(size(sp,2)));
%! fprintf(1, 'Converted into crs_matrix in %g seconds\n', toc);
%! fprintf(1, 'Testing serial: ');
%! b1 = crs_prodAx( A, x);
%! assert( norm(b0-b1)<=1.e-12);

%! b2 = zeros(size(sp,1),2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>MACC_get_max_threads; break; end
%!     fprintf(1, 'Testing %d threads: ', nthreads);
%!     b2 = crs_prodAx( A, x, b2, nthreads);
%!     assert( norm(b0-b2)<=1.e-12);
%! end

%! if ~MMPI_Initialized; MMPI_Init; end
%! nprocs = double(MMPI_Comm_size(MPI_COMM_WORLD));
%! fprintf(1, 'Testing 2 threads with MPI call: ');
%! b2 = crs_prodAx( A, x, b2, int32(2), MPI_COMM_WORLD);
%! assert( nprocs>1 || norm(b0-b2)<=1.e-12);

%! b3 = zeros(size(sp,1)+1,1);
%! fprintf(1, 'Testing 2 threads with piggy-back: ');
%! b3 = crs_prodAx( A, x(:,1), b3, int32(2), MPI_COMM_WORLD, 1, int32(1));
%! assert( nprocs>1 || norm(b0(:,1)-b3(1:end-1))<=1.e-12);
%! assert( b3(end)==nprocs);
