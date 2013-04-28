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
% It will use min( nthreads, size(b,1)/A.ncols, MACC_get_max_threads) threads.
%
%      b = crs_prodAtx( A, x, b, [])
% Computes b=A'*x locally within each OpenMP thread, assuming the parallel
% region has already been initialized. The argument b has been preallocated.
% It will use min( size(b,1)/A.ncols, MACC_get_mnum_threads)
% threads forcomputations. In this mode, there is no implicit barrier.
%
%      b = crs_prodAtx( A, x, b, nthreads, comm)
% also performs communication within the MPI communicator.
% It assumes that A is partitioned rowwise, and hence A'*x is
% a partial sum on each process. An allreduce is performed on b
% within the communicator.
%
% If A is partitioned columnwise, then b must be duplicated on all 
% processes, and crs_prodAtx should be called without the communicator.
%
%      b = crs_prodAtx( A, x, b, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in b. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

% See http://www.netlib.org/linalg/html_templates/node98.html

%#codegen -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), 
%#codegen coder.typeof( int32(1), [1,1], [1,0])}
%#codegen crs_prodAtx_ser -args {crs_matrix, coder.typeof(0, [inf,inf])}
%#codegen crs_prodAtx_ser1 -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen crs_prodAtx_mpi -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), opaque_obj}
%#codegen crs_prodAtx_mpip -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), opaque_obj, coder.typeof(0, [inf,1])}
%#codegen crs_prodAtx_mpip1 -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), opaque_obj,coder.typeof(0, [inf,1]), int32(0)}

coder.inline('never');

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
    if ~MACC_get_nested && ismt && nthreads(1)>1
        MACC_begin_master
        msg_warn('crs_prodAtx:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        MACC_end_master
    end
    
    MACC_begin_parallel; MACC_clause_default('shared'); 
    MACC_clause_num_threads( int32(nthreads(1)));
    
    %% Compute b=A'*x on each thread
    b = crs_prodAtx_internal( A.row_ptr, A.col_ind, A.val, ...
        x, b, MACC_get_thread_num*A.ncols, A.nrows, A.ncols, ismt);
    
    %% Perfom summation of partial sums in b
    if ismt
        MACC_barrier;
        b = accu_partsum( b, A.ncols);
    end
    
    MACC_end_parallel
else
    %% Compute b=A'*x
    b = crs_prodAtx_internal( A.row_ptr, A.col_ind, A.val, ...
        x, b, MACC_get_thread_num*A.ncols, A.nrows, A.ncols, ismt);
    
    if ismt
        %% Perfom summation of partial sums in b
        MACC_barrier;
        b = accu_partsum( b, A.ncols);
    end

    if ~isempty(varargin) && ismt; MACC_barrier; end
end

if ~isempty(varargin)
    % Perform MPI allreduce
    MACC_begin_single
    s = int32(A.ncols*size(x,2));
    b = allreduce(b, s, MPI_SUM, varargin{:});
    MACC_end_single
end

function b = crs_prodAtx_internal( row_ptr, col_ind, val, x, ...
    b, offset, nrows, ncols, ismt)

if isempty( coder.target)
    for k=1:int32(size(x,2))
        for j=offset+1:offset+ncols; b(j,k) = 0; end
        
        for i=1:nrows
            alpha = x(i, k);
            n = row_ptr(i+1) - 1;
            for j = row_ptr(i) : n
                r = offset+col_ind(j);
                b(r, k) = b(r, k) + alpha * val(j);
            end
        end
    end
else
    coder.inline('never');
    coder.cinclude('spalab_kernel.h');

    if ismt
        nthreads = min(MACC_get_num_threads, ...
            int32(floor(double(size(b,1))/double(ncols))));
        [istart, iend] = get_local_chunk(nrows, [], nthreads);
    else
        istart = int32(1); iend = int32( nrows);
    end
    
    if isa( val, 'double'); func = 'SPL_daxpy';
    else func = 'SPL_saxpy'; end
    
    for k=1:int32(size(x,2))
        for j=offset+1:offset+ncols; b(j,k) = 0; end

        for i=istart:iend
            coder.ceval( func, x(i, k), coder.rref(val( row_ptr(i))), ...
                coder.ref( b(offset+1,k)), coder.rref( col_ind(row_ptr(i))), ...
                row_ptr(i+1)-row_ptr(i));
        end
    end
end

function b = accu_partsum( b, ncols)
% This function should be called only in multi-threading mode

if size(b,1)>=ncols*2
    [istart, iend] = get_local_chunk(ncols);
    nthreads = min(MACC_get_num_threads, ...
        int32(floor(double(size(b,1))/double(ncols))));
    
    offset = ncols;
    for j=2:nthreads
        for k=1:int32(size(b,2))
            for i=istart:iend
                b(i,k) = b(i,k) + b(offset+i,k);
            end
        end
        offset = offset + ncols;
    end
end

function test %#ok<DEFNU>
%!test
%! if ~exist( ['crs_prodAtx.' mexext], 'file')
%!    m=100; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; sp = sprand(m,n,0.5); x=rand(size(sp,1),2);
%! [is,js,vs] = find(sp);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = sp'*x;
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! tic; A = crs_matrix(int32(is), int32(js), vs, int32(size(sp,1)), int32(size(sp,2)));
%! fprintf(1, '\tConverted into crs_matrix in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = crs_prodAtx( A, x);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( norm(b0-b1)/norm(b0)<=1.e-10);

%! for nthreads=int32([1 2 4 8])
%!     if nthreads>MACC_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d threads: ', nthreads);
%!     tic; b2 = zeros(size(sp,2)*nthreads,2);
%!     b2 = crs_prodAtx( A, x, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert( norm(b0-b2(1:size(sp,2),:))/norm(b0)<=1.e-10);
%! end

%! nprocs = double(comm_size(MPI_COMM_WORLD));
%! fprintf(1, '\tTesting 2 threads with MPI call: ');
%! b2 = zeros(size(sp,2)*2,2);
%! tic; b2 = crs_prodAtx( A, x, b2, int32(2), MPI_COMM_WORLD);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( norm(b0-b2(1:size(sp,2),:))/norm(b0)<=1.e-10);

%! b3 = zeros(2*size(sp,2)+1,1);
%! fprintf(1, '\tTesting 2 threads with piggy-back: ');
%! tic; b3 = crs_prodAtx( A, x(:,1), b3, int32(2), MPI_COMM_WORLD, 1, int32(1));
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( nprocs>1 || norm(b0(:,1)-b3(1:size(sp,2)))/norm(b0)<=1.e-10);
%! assert( b3(size(sp,2)+1)==nprocs);
