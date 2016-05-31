function [mat, errCode] = cuBCRSCreate(mb, nb, nnzb, blkdim, type)
%Creates a sparse matrix in block CRS format on a CUDA device.
%
%  [mat, errCode] = cuBCRSCreate(mb, nb, nnzb, [blkdim], [type]) creates a
%  sparse block-matrix of size mb-by-nb on CUDA. The block sizes are
%  blkdim-by-blkdim. If blkdim is not specified, then it is assumed to be
%  1. The type is MSP_DOUBLE, MSP_SINGLE, MSP_COMPLEX, MSP_DOUBLE_COMPLEX,
%  MSP_INT*, or MSP_UINT*. If not  specified, then type is MSP_DOUBLE.
%
%  SEE ALSO: cuBCRSDestroy, cuBCRSCopyToGPU, cuBCRSCopyFromGPU

%#codegen -args {int32(0), int32(0), int32(0), int32(0), int32(0)}
%#codegen cuBCRSCreate_3args -args {int32(0), int32(0), int32(0)}
%#codegen cuBCRSCreate_4args -args {int32(0), int32(0), int32(0), int32(0)}

coder.cinclude('mspack.h');

if nargin<4; blkdim = int32(1); end
if nargin<5; type = MSP_DOUBLE; end

[rowptr, errCode] = cuVecCreate(mb+1, MSP_INT32);
if ~errCode;
    [colind, errCode] = cuVecCreate(nnzb, MSP_INT32);
else
    colind = CuVec(0);
end
if ~errCode;
    [vals, errCode] = cuVecCreate(nnzb*blkdim*blkdim, type);
else
    vals = CuVec(0);
end

mat = CuBCRS(rowptr, colind, vals, type, mb, nb, nnzb, blkdim);

end
