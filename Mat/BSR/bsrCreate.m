function mat = bsrCreate(mb, nb, nnzb, blkdim, type_in)
%Creates a sparse matrix in block CRS format.
%
%  [mat, errCode] = bsrCreate(mb, nb, nnzb, [blkdim], [type]) creates a
%  sparse block-matrix of size mb-by-nb. The block sizes are blkdim-by-blkdim,
%  stored in column major. If blkdim is not specified, then it is assumed to
%  be 1. The type can be an integer (MCU_DOUBLE etc.) or a character string
%  for MATLAB type 'double'. If not present, then type is MCU_DOUBLE.

if nargin<4; blkdim = int32(1); end
if nargin<5;
    type = MCU_DOUBLE;
elseif ischar(type_in)
    type = mcuType(type_in, true);
else
    type = type_in;
end

rowptr = coder.nullcopy(zeros(mb+1, 1, 'int32'));
colind = coder.nullcopy(zeros(nnzb, 1, 'int32'));
vals = coder.nullcopy(repmat(mcuZero(type),nnzb*blkdim*blkdim, 1));

mat = MatBSR(rowptr, colind, vals, int32(mb), int32(nb), int32(nnzb), int32(blkdim));
end
