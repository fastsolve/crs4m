function [mat, errCode] = cuCRSCreate(m, n, nnz, varargin)
%Creates a sparse matrix in CRS format on a CUDA device.
%
%  [mat, errCode] = cuCRSCreate(m, n, nnz, [type]) creates a sparse
%  matrix of size m-by-n on CUDA, where type is MSP_DOUBLE, MSP_SINGLE,
%  MSP_COMPLEX, MSP_DOUBLE_COMPLEX, MSP_INT*, or MSP_UINT*. If not
%  specified, then type is MSP_DOUBLE.
%
%  SEE ALSO: cuCRSDestroy, cuCRSCopyToGPU, cuCRSCopyFromGPU

%#codegen -args {int32(0), int32(0), int32(0), int32(0)}
%#codegen cuCRSCreate_3args -args {int32(0), int32(0), int32(0)}

coder.cinclude('mspack.h');

if nargin<4; type = MSP_DOUBLE; end

[rowptr, errCode] = cuVecCreate(m+1, MSP_INT32);
if ~errCode;
    [colind, errCode] = cuVecCreate(nnz, MSP_INT32);
else
    colind = CuVec(0);
end
if ~errCode;
    [vals, errCode] = cuVecCreate(nnz, type);
else
    vals = CuVec(0);
end

mat = CuCRS(rowptr, colind, vals, type, m, n, nnz);

end
