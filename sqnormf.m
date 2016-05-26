function sqnrm = sqnormf(A, sqnrm, nthreads, varargin)
%SQNORMF Compute the squared Frobenius norm of A, i.e. trace(A'*A).
%       sqnrm = sqnormf(A)
% computes the squared norm in serial.
%
%       sqnrm = sqnormf(A, sqnrm, nthreads)
% uses nthreads OpenMP threads, where nthreads is an int32.
% Argument sqnrm can be an array, but only its first entry will be accessed.
%
%       sqnrm = sqnormf(A, sqnrm, [])
% computes the squared norm on a local chunk per thread. 
% This assumes that a parallel region has already been started.
% This mode is not supported if sqnormf is a top-level mex function. 
% In this mode, argument sqnrm must be a variable shared by all threads in the team. 
%
%       sqnrm = sqnormf(A, sqnrm, nthreads, comm)
% also performs an allreduce on the trace within the MPI communicator.
%
%      sqnrm = sqnorm(A, sqnrm, mc, comm, s)
% specifies that s additional entries in prod should be summed
% during allreduce. This is to piggy-back on allreduce to
% help aggregate small messages to optimize MPI communication.
%
% In multithread or multiprocess mode, this is a collective function.
% In multithread mode, there is an implicit barrier at the beginning and
% at the end of this function.
% 
% The compiled version reqiures A to be double precision.
% The M file works for both single and double precision.

%#codegen -args {coder.typeof(0, [inf,inf]), 0, coder.typeof(int32(0), [1,1], [1,0])}
%#codegen sqnormf_ser -args {coder.typeof(0, [inf,inf])}
% #codegen sqnormf_mpi -args {coder.typeof(0, [inf,inf]), 0, coder.typeof(int32(0), [1,1], [1,0]), MPI_Comm}

coder.inline('never');

iend=int32(numel(A)); 
if nargin<2; sqnrm=0; end
ismt = nargin>2 && momp_get_num_threads>1;

%% Declare parallel region
if nargin>2 && ~isempty(nthreads)
    if ~momp_get_nested && ismt && nthreads(1)>1
        momp_begin_master
        m2c_warn('sqnormf:NestedParallel', ...
            'You are trying to use nested parallel regions, but nested parallelism is not enabled.');
        momp_end_master
    end
    i = int32(0); s_local = cast(0, class(A));
    [i, s_local, iend] = momp_begin_parfor(i, s_local, iend);
    momp_clause_default('none');
    momp_clause_shared(A, iend); momp_clause_private(i); 
    momp_clause_reduction('+', s_local); 
    momp_clause_num_threads(nthreads(1));    
    
    %% Compute sqnrm=trace(A'*A)
    for i=1:iend
        s_local = s_local + A(i)*A(i);
    end
    
    %% End parallel region
    momp_end_parfor;
    sqnrm(1) = s_local;
else
    if ismt
        momp_barrier; momp_begin_single
        sqnrm(1) = 0;
        momp_end_single
        [istart, iend] = get_local_chunk(iend);
    else
        sqnrm(1) = 0; istart = int32(1);
    end

    %% Compute sqnrm=trace(A'*A)
    s_local = 0;
    for i=istart:iend
        s_local = s_local + A(i)*A(i);
    end    

    if nargin>1 && ismt
        %% Perform reduction across all threads
        momp_atomic
        sqnrm(1) = sqnrm(1) + s_local;
        momp_barrier
    else
        sqnrm(1) = s_local;
    end
end

if ~isempty(varargin)
    % Perform MPI allreduce
    momp_begin_single
    sqnrm = allreduce(sqnrm, int32(1),MPI_SUM, varargin{:});
    momp_end_single
end

function test %#ok<DEFNU>
%!test
%! if ~exist(['sqnormf.' mexext], 'file')
%!    m=200; n = 20;
%! else
%!    m=10000; n = 2000;
%! end
%! tic; A = rand(m,n); x=rand(size(A,2),2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = norm(A,'fro').^2;
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = sqnormf(A);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(abs(b0-b1)/abs(b0)<=1.e-12);

%! b2 = nan;
%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     tic; b2 = sqnormf(A, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(abs(b0-b2)/abs(b0)<=1.e-12);
%! end

% %! fprintf(1, '\tTesting 2 threads with MPI call: ');
% %! tic; b2 = sqnormf(A, b2, int32(2), MPI_COMM_SELF);
% %! fprintf(1, 'Done in %g seconds\n', toc);
% %! assert(abs(b0-b2)/abs(b0)<=1.e-12);
% 

