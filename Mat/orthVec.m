function [ z, buf ] = orthVec(Q, z, ncols, buf, nthreads, varargin) %#codegen
%ORTHVEC  Orthogalize a given vector z by leading columns of Q.
%
%   z = orthVec(Q, z), z = orthVec(Q, z, ncols)
% Perform orthgonalization in serial. Argument ncols specifies the number of
% columns of Q to be used, and these columns are assumed to be unit vectors.
%
%   [ z, buf ] = orthVec(Q, z, ncols, buf )
% Use buf as intermediate buffers for compuation. The number of entries in
% buf must be >= ncols.
%
%   [ z, buf ] = orthVec(Q, z, ncols, buf, nthreads)
% Use nthreads for computation. The number of entries in
% buf must be >= ncols.
%
%   [ z, buf ] = orthVec(Q, z, ncols, buf, nthreads, comm )
% Use allreduce across the MPI communicator.
%
%   [ z, buf ] = orthVec(Q, z, ncols, buf, nthreads, comm, pgmsg, [pgsize] )
% Aggregate piggy-back message when performing allreduce.

coder.inline('never');

if nargin<3; ncols = int32(size(Q,2)); end
if nargin<4; 
    buf = nullcopy(zeros(ncols, 1)); 
elseif numel(buf)<ncols
    m2c_error('orthVec:BufferOverflow', 'Buffer space buf is too small.');
end

nt = momp_get_num_threads;

% Perform re-orthoganlization of z using Gram-Schmidt (part 1).
if nargin<5 || ~isempty(nthreads) || isempty(nthreads) && nt<=1;
    for k=1:ncols
        t = 0;
        for j=int32(1):int32(size(Q,1))
            t = t + Q(j,k)*z(j);
        end
        buf(k) = t;
    end
else
    [rstart, rend] = get_local_chunk(ncols);
    for k=rstart:rend
        t = 0;
        for j=int32(1):int32(size(Q,1))
            t = t + Q(j,k)*z(j);
        end
        buf(k) = t;
    end
    
    momp_barrier;
end

if ~isempty(varargin)
    % Perform MPI reduction.
    momp_begin_single
    buf = allreduce(buf, ncols, MPI_SUM, varargin{:});
    momp_end_single
end

% Perform re-orthoganlization of z using Gram-Schmidt (part 2).
% This part is fully parallel
[cstart, cend] = get_local_chunk(int32(size(Q,1)));
for j=cstart:cend
    t = z(j);
    for k=1:ncols
        t = t - buf(k)*Q(j,k);
    end
    z(j) = t;
end
