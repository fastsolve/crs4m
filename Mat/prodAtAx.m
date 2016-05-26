function [b, Ax] = prodAtAx(A, x, b, Ax, nthreads, varargin)
%prodAtAx Compute A'*A*x.
%
%      b = prodAtAx(A, x [, b]) or [b, Ax] = prodAtAx(A, x, b, Ax)
% computes b=A'*A*x in serial. Ax is a buffer for storing the 
% intermediate result Ax. If passed as input, b and Ax are assumed 
% to have been preallocated.
%
%      [b, Ax] = prodAtAx(A, x, b, Ax, nthreads)
% computes b=A'*A*x using nthreads OpenMP threads.
%
%      [b, Ax] = prodAtAx(A, x, b, Ax, [])
% computes b=A'*A*x locally within each OpenMP thread, assuming the 
% parallel region has already been initialized. In this mode, Ax must
% be a buffer shared by all threads in the team, and it must be both 
% input and output arguments (i.e., passed by reference).
%
%      [b, Ax] = prodAtAx(A, x, b, Ax, nthreads, comm)
% assumes that A is partitioned rowwise, and x is dupliplicate on all
% processes. It performs a reduction on b within the MPI communicator.
%
%      [b, Ax] = prodAtAx(A, x, b, Ax, nthreads, comm, pgmsg [, pgsz])
% specifies an additiona message that should be summed during
% allreduce and strored in b. This is to piggy-back on allreduce
% to help aggregate small messages to optimize MPI communication.
%
% In multithread or multiprocess mode, this is a collective function, 
% and there is an implicit barrier at the end of this function.
% 
% The compiled version reqiures A and x to be double precision.
% The M file works for both single and double precision.

%#codegen -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), 
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0])}
%#codegen prodAtAx_ser -args {coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf])}
%#codegen prodAtAx_ser1 -args {coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}
% %#codegen prodAtAx_mpi -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), 
% %#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0]), MPI_Comm}
% %#codegen prodAtAx_mpip -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), 
% %#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0]), MPI_Comm, coder.typeof(0, [inf,inf])}
% %#codegen prodAtAx_mpip1 -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), 
% %#codegen coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(0),[1,1],[1,0]), MPI_Comm, coder.typeof(0, [inf,inf]), int32(0)}

coder.inline('never');

assert(nargin<=8);

if nargin<3
    b = zeros(size(A,2),size(x,2)); 
elseif size(b,1)<size(A,2) || size(b,2)~=size(x,2)
    momp_begin_master
    m2c_error('prodAtAx:IncorrectBuffer', 'Buffer b has incorrect size.');
    momp_end_master
end
if nargin<4; 
    Ax = zeros(size(A,1),size(x,2));
elseif size(Ax,1)<size(A,1) || size(Ax,2)~=size(x,2)
    momp_begin_master
    m2c_error('prodAtAx:IncorrectBuffer', 'Buffer Ax has incorrect size.');
    momp_end_master
end
ismt = nargin>=5 && momp_get_num_threads>1;

if nargin>=4 && nargout<2
    momp_begin_master
    m2c_warn('prodAtAx:MissingBuffer', ...
        'prodAtAx is called within a parallel region but Ax is not an in+out argument.');
    momp_end_master
end

%% Declare parallel region
if nargin>=5 && ~isempty(nthreads)
    if ~momp_get_nested && ismt && nthreads(1)>1
        momp_begin_master
        m2c_warn('prodAtAx:NestedParallel', ...
            'You are trying to use nested parallel regions, but nested parallelism is not enabled.');
        momp_end_master
    end

    [b, Ax] = momp_begin_parallel(b, Ax);
    momp_clause_default('shared');
    momp_clause_num_threads(nthreads(1));
    
    % Computes Ax=A*x
    Ax = prodAx(A, x, Ax, []);
    momp_barrier;
    
    % Computes b=A'*Ax
    b = prodAtx(A, Ax, b, [], varargin{:});
    
    momp_end_parallel(Ax);
elseif ismt
    % Computes Ax=A*x
    Ax = prodAx(A, x, Ax, []);
    momp_barrier;
    
    % Computes b=A'*Ax
    b = prodAtx(A, Ax, b, [], varargin{:});
else
    % Computes Ax=A*x
    Ax = prodAx(A, x, Ax);
    % Computes b=A'*Ax
    b = prodAtx(A, Ax, b);
end

function test %#ok<DEFNU>
%!test
%! if ~exist(['prodAtAx.' mexext], 'file')
%!    m=200; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; A = rand(m,n); x=rand(size(A,2),2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = A'*(A*x);
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = prodAtAx(A, x);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(norm(b0-b1)/norm(b0)<=1.e-12);

%! b2 = zeros(size(A,2),2);
%! Ax = zeros(size(A,1),2);
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; [b2, Ax] = prodAtAx(A, x, b2, Ax, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2)/norm(b0)<=1.e-12);
%! end

% %! fprintf(1, '\tTesting 2 threads with MPI call: ');
% %! tic; [b2, Ax] = prodAtAx(A, x, b2, Ax, int32(2), MPI_COMM_SELF);
% %! fprintf(1, 'Done in %g seconds\n', toc);
% %! assert(norm(b0-b2)/norm(b0)<=1.e-12);
% 
% %! b3 = zeros(size(A,2)+1,1);
% %! fprintf(1, '\tTesting 2 threads with piggy-back: ');
% %! x = x(:,1); Ax = zeros(size(A,1),1);
% %! tic; [b3, Ax] = prodAtAx(A, x, b3, Ax, int32(2), MPI_COMM_SELF, 1, int32(1));
% %! fprintf(1, 'Done in %g seconds\n', toc);
% %! assert(norm(b0(:,1)-b3(1:end-1))/norm(b0(:,1))<=1.e-12);
% %! assert(b3(end)==1);
