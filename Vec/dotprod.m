function prod = dotprod(u, v, prod, nthreads)
%DOTPROD  computes dot prod of column vectors.
%
%      prod = dotprod(u, v)
% computes the dot product in serial. Both u and v are assumed to be
% column vectors. If they have multiple columns, then compute the dot
% product of each pair of columns and save the answer into columns of prod.
%
%      prod = dotprod(u, v, prod, nthreads)
% computes the dot product using multicores, as controlled by nthreads.
% Note that size(prod,1) should be >= nthreads to facilitate
% concurrency. Otherwise, it will use up to size(prod,1) threads.
%
% If nthreads is empty, then it computes prod=u'*v on a local chunk per
% thread. This assumes that a parallel region has already been started.
% In this mode, prod must be a variable shared by all threads in the team,
% and size(prod,1) should be >= nthreads to facilitate concurrency.
% Otherwise, it will use up to size(prod,1) threads.
%
% In multithread, this is a collective function.
% In multithread mode, there is an implicit barrier at the beginning and
% at the end of this function.
%
% Note: The compiled version requires u and v to be double precision.
%       The M file works for both single and double precision.
%
% See also momp_Control MPI_Comm

%#codegen -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]),
%#codegen coder.typeof(0, [inf,inf]), coder.typeof(int32(1), [1,1], [1,0])}
%#codegen dotprod_ser -args {coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}

coder.inline('never');

if size(u,1)~=size(v,1) || size(u,2)~=size(v,2)
    momp_begin_master
    m2c_error('dotprod:IncorrectSize', 'Dimensions for u and v must match.');
    momp_end_master
end

if nargin<3
    prod=nullcopy(zeros(1,size(u,2)));
elseif size(prod,1)<1 || size(prod,2)~=size(u,2)
    momp_begin_master
    m2c_error('dotprod:IncorrectSize', 'prod and u must have the same number of columns.');
    momp_end_master
end

ismt = nargin>3 && momp_get_num_threads>1;

if nargin>3 && ~isempty(nthreads)
    %% Declare parallel region
    if ~momp_get_nested && momp_get_num_threads>1 && nthreads(1)>1
        momp_begin_master
        m2c_warn('dotprod:NestedParallel', ...
            'You are trying to use nested parallel regions, but nested parallelism is not enabled.');
        momp_end_master
    end
    
    momp_begin_parallel; momp_clause_default('shared');
    momp_clause_num_threads(nthreads(1));
    
    %% Compute prod_local=u'*v
    prod = dotprod_partial(u, v, prod, '');
    momp_barrier
    
    prod = accu_partsum(prod);    
    %% End parallel region
    momp_end_parallel;
elseif ismt
    %% Compute prod_local=u'*v
    prod = dotprod_partial(u, v, prod, '');
    momp_barrier
    
    prod = accu_partsum(prod);
    momp_barrier
else
    prod = dotprod_partial(u, v, prod);
end

function prod = dotprod_partial(u, v, prod, varargin)
% Compute partial dot product with each therad

if nargin>3
    coder.inline('never');
    ind = momp_get_thread_num+1;
    [istart, iend] = get_local_chunk(int32(size(u,1)));
else
    ind = int32(1);
    istart = int32(1); iend = int32(size(u, 1));
end

for k=1:int32(size(u,2))
    prod(ind, k) = 0.;
    for i=istart:iend
        prod(ind, k) = prod(ind, k) + u(i, k) * v(i, k);
    end
end

function prod = accu_partsum(prod)
% Accumulate the partial sums of all columns
coder.inline('always');

[istart, iend] = get_local_chunk(int32(size(prod,2)));
n = min(momp_get_num_threads, int32(size(prod,1)));

for k=istart:iend
    for j=2:n
        prod(1,k) = prod(1,k) + prod(j,k);
    end
end

function test %#ok<DEFNU>
%!test
%! if ~exist(['dotprod.' mexext], 'file')
%!    m=200;
%! else
%!    m=1000000;
%! end
%! tic; u =rand(m,2); v=rand(m,2);
%! fprintf(1, '\n\tGenerated random matrix in %g seconds\n', toc);
%! tic; b0 = [u(:,1)'*v(:,1), u(:,2)'*v(:,2)];
%! fprintf(1, '\tComputed reference solution in %g seconds\n', toc);
%! fprintf(1, '\tTesting serial: ');
%! tic; b1 = dotprod(u, v);
%! fprintf(1, 'Done in %g seconds\n ', toc);
%! assert(norm(b0-b1)/norm(b0)<=1.e-6);

%! for nthreads=int32([1 2 4 8])
%!     if nthreads>momp_get_max_threads; break; end
%!     fprintf(1, '\tTesting %d thread(s): ', nthreads);
%!     b2 = zeros(nthreads,2);
%!     tic; b2 = dotprod(u, v, b2, nthreads);
%!     fprintf(1, 'Done in %g seconds\n ', toc);
%!     assert(norm(b0-b2(1,:))/norm(b0)<=1.e-6);
%! end
