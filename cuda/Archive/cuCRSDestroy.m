function [mat, errCode] = cuCRSDestroy(mat)
%Destroys a sparse matrix in CRS format on a CUDA device.
%
%  [mat, errCode] = cuCRSDestroy(mat) destroies the arrays pointed to by mat.
%
%  SEE ALSO: cuCRSCreate, cuCRSCopyToGPU, cuCRSCopyFromGPU

coder.inline('always');
[mat, errCode] = cuBCRSDestroy(mat);