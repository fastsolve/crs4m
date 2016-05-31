function [mat, errCode] = cuBCRSDestroy(mat)
%Destroys a sparse matrix in block CRS format on a CUDA device.
%
%  [mat, errCode] = cuBCRSDestroy(mat) destroies the arrays pointed to by mat.
%
%  SEE ALSO: cuBCRSCreate, cuBCRSCopyToGPU, cuBCRSCopyFromGPU

%#codegen -args {CuBCRS}

coder.cinclude('mspack.h');

[mat.vals, errCode] = cuVecDestroy(mat.vals);
if ~errCode; [mat.colind, errCode] = cuVecDestroy(mat.colind); end
if ~errCode; [mat.rowptr, errCode] = cuVecDestroy(mat.rowptr); end
end
