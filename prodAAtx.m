function [b, Atx] = prodAAtx( A, x, b, Atx, nthreads, varargin)
%prodAAtx Compute A*A'*x.
%
%      b = prodAAtx( A, x), b = prodAAtx( A, x, b), or 
%      [b, Atx] = prodAAtx(A, x, b, Atx)
% computes b=A*A'*x in serial. Atx is a buffer for storing the 
% intermediate result A'*x. If passed as input, b and Atx are 
% assumed to have been preallocated.
%
%      [b, Atx] = prodAAtx( A, x, b, Atx, nthreads)
% computes b=A*A'*x using nthreads OpenMP threads, where nthreads is an int32.
%
%      [b, Atx] = prodAAtx( A, x, b, Atx, [])
% computes b=A*A'*x locally within each OpenMP thread, assuming the
% parallel region has already been initialized. In this mode, Atx must
% be a buffer shared by all threads in the team, and it must be both 
% input and output arguments (i.e., passed by reference).
%
%      [b, Atx] = prodAAtx( A, x, b, Atx, nthreads, comm)
% assumes that A is partitioned rowwise, so that Atx will be composed of
% partial sums on each MPI process. It performs a reduction on Atx
% within the MPI communicator before computing A*Atx.
%
%      [b, Atx] = prodAAtx( A, x, b, Atx, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in Atx. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% In multithread or multiprocess mode, this is a collective function, 
% and there is an implicit barrier at the end of this function.
% 
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

%#codegen -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0])}
%#codegen prodAAtx_ser -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen prodAAtx_ser1 -args {coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
%#codegen prodAAtx_mpi -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0]), MPI_Comm}
%#codegen prodAAtx_mpip -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0]), MPI_Comm, coder.typeof(0, [inf,1])}
%#codegen prodAAtx_mpip1 -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0]), MPI_Comm, coder.typeof(0, [inf,1]), int32(0)}

coder.inline('never');

assert( nargin<=8);

if nargin<3
    b = zeros(size(A,1),size(x,2));
elseif size(b,1)<size(A,1) || size(b,2)~=size(x,2)
    MACC_begin_master
    m2c_error('prodAAtx:IncorrectBuffer', 'Buffer b has incorrect size.');
    MACC_end_master
end
if nargin<4
    Atx = zeros(size(A,2),size(x,2));
elseif size(Atx,1)<size(A,2) || size(Atx,2)~=size(x,2)
    MACC_begin_master
    m2c_error('prodAAtx:IncorrectBuffer', 'Buffer Atx has incorrect size.');
    MACC_end_master
end
ismt = nargin>=5 && MACC_get_num_threads>1;

if nargin>=4 && nargout<2
    MACC_begin_master
    m2c_warn('prodAAtx:MissingBuffer', ...
        'prodAAtx is called Atx as input but Atx not an out argument.');
    MACC_end_master
end

%% Declare parallel region
if nargin>=5 && ~isempty( nthreads)
    if ~MACC_get_nested && ismt && nthreads(1)>1
        MACC_begin_master
        m2c_warn('prodAAtx:NestedParallel', ...
            'You are trying to use nested parallel regions, but nested parallelism is not enabled.');
        MACC_end_master
    end
    
    [b, Atx] = MACC_begin_parallel( b, Atx);
    MACC_clause_default( 'shared');
    MACC_clause_num_threads( int32(nthreads(1)));
    
    % Computes Atx=A'*x
    Atx = prodAtx( A, x, Atx, [], varargin{:});
    if isempty( varargin); MACC_barrier; end
    % Computes b=A*Atx
    b = prodAx( A, Atx, b, []);
    
    MACC_end_parallel(Atx);
elseif ismt
    % Computes Atx=A'*x
    Atx = prodAtx( A, x, Atx, [], varargin{:});
    if isempty( varargin); MACC_barrier; end
    
    % Computes b=A*Atx
    b = prodAx( A, Atx, b, []);
else
    % Computes Atx=A'*x
    Atx = prodAtx( A, x, Atx);
    % Computes b=A*Atx
    b = prodAx( A, Atx, b);    
end

function test %#ok<DEFNU>
%!test
%! if ~exist( ['prodAAtx.' mexext], 'file')
%!    m=200; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; A = rand(m,n); x=rand(size(A,1),2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = A*(A'*x);
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = prodAAtx( A, x);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert( norm(b0-b1)/norm(b0)<=1.e-12);

%! b2 = zeros(size(A,1),2);
%! Atx = zeros(size(A,2),2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>MACC_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; [b2, Atx] = prodAAtx( A, x, b2, Atx, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert( norm(b0-b2)/norm(b0)<=1.e-12);
%! end

%! fprintf(1, '\tTesting 2 threads with MPI call: ');
%! tic; [b2, Atx] = prodAAtx( A, x, b2, Atx, int32(2), MPI_COMM_SELF);
%! fprintf(1, 'Done in %g seconds\n', toc);
%! assert( norm(b0-b2)/norm(b0)<=1.e-12);

%! b3 = zeros(size(A,1),1); Atx = zeros(size(A,2)+1,1);
%! fprintf(1, '\tTesting 2 threads with piggy-back: ');
%! tic; [b3, Atx] = prodAAtx( A, x(:,1), b3, Atx, int32(2), MPI_COMM_SELF, 1, int32(1));
%! fprintf(1, 'Done in %g seconds\n', toc);
%! assert( norm(b0(:,1)-b3)/norm(b0(:,1))<=1.e-12);
%! assert( Atx(size(A,2)+1)==1);
