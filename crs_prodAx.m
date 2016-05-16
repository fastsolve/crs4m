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

%#codegen -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof( int32(1), [1,1], [1,0])}
%#codegen crs_prodAx_ser -args {crs_matrix, coder.typeof(0, [inf,inf])}
%#codegen crs_prodAx_ser1 -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen crs_prodAx_mpi -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), MPI_Comm}
%#codegen crs_prodAx_mpip -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), MPI_Comm, coder.typeof(0, [inf,1])}
%#codegen crs_prodAx_mpip1 -args {crs_matrix, coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof( int32(1), [1,1], [1,0]), MPI_Comm, coder.typeof(0, [inf,1]), int32(0)}

coder.inline('never');

if nargin==2;
    b = nullcopy(zeros(A.nrows,size(x,2)));
else
    if size(b,1)<A.nrows || size(b,2)<size(x,2)
        m2c_error('crs_prodAx:BufferTooSmal', 'Buffer space for output b is too small.');
    end
end
ismt = nargin>=4 && pACC_get_num_threads>1;

if nargin>3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~pACC_get_nested && ismt && nthreads(1)>1
        pACC_begin_master
        m2c_warn('crs_prodAx:NestedParallel', ...
            'You are trying to use nested parallel regions. Solution may be incorrect.');
        pACC_end_master
    end
    
    pACC_begin_parallel; pACC_clause_default('shared');
    pACC_clause_num_threads( int32(nthreads(1)));
    
    %% Compute b=A*x
    b = crs_prodAx_kernel( A.row_ptr, A.col_ind, A.val, x, int32(size(x,1)), ...
        b, int32(size(b,1)), A.nrows, int32(size(x, 2)), pACC_get_num_threads>1);
    
    %% End parallel region
    pACC_end_parallel;
else
    %% Compute b=A*x
    b = crs_prodAx_kernel( A.row_ptr, A.col_ind, A.val, x, int32(size(x,1)), ...
        b, int32(size(b,1)), A.nrows, int32(size(x, 2)), ismt);
end

if ~isempty(varargin)
    if isempty(nthreads) && ismt; pACC_barrier; end
    
    % Perform MPI allreduce
    pACC_begin_single
    s = int32(A.nrows*size(x,2));
    b = allreduce(b, s, MPI_SUM, varargin{:});
    pACC_end_single
end

function b = crs_prodAx_kernel( row_ptr, col_ind, val, ...
    x, x_m, b, b_m, nrows, nrhs, ismt)

pACC_kernel_function

coder.inline('never');
if ismt
    [istart, iend] = get_local_chunk(nrows);
else
    istart = int32(1); iend = nrows;
end

xoffset = int32(0); boffset = int32(0);
for k=1:nrhs
    for i=istart:iend
        t = 0.0;
        for j=row_ptr(i):row_ptr(i+1)-1
            t = t + val(j)*x(xoffset+col_ind(j));
        end
        b(boffset+i) = t;
    end
    xoffset = xoffset + x_m; boffset = boffset + b_m;
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
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = sp*x;
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! tic; A = crs_matrix(int32(is), int32(js), vs, int32(size(sp,1)), int32(size(sp,2)));
%! fprintf(1, '\tConverted into crs_matrix in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = crs_prodAx( A, x);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( norm(b0-b1)<=1.e-12);

%! b2 = zeros(size(sp,1),2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>pACC_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; b2 = crs_prodAx( A, x, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert( norm(b0-b2)<=1.e-12);
%! end

%! nprocs = double(comm_size(MPI_COMM_WORLD));
%! fprintf(1, '\tTesting 2 threads with MPI call: ');
%! tic; b2 = crs_prodAx( A, x, b2, int32(2), MPI_COMM_WORLD);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( nprocs>1 || norm(b0-b2)<=1.e-12);

%! b3 = zeros(size(sp,1)+1,1);
%! fprintf(1, '\tTesting 2 threads with piggy-back: ');
%! tic; b3 = crs_prodAx( A, x(:,1), b3, int32(2), MPI_COMM_WORLD, 1, int32(1));
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( nprocs>1 || norm(b0(:,1)-b3(1:end-1))<=1.e-12);
%! assert( b3(end)==nprocs);
