function [errCode, toplevel] = cudaVecCopyFromHost(vec, cuVec)
% Copies a vector from MATLAB to CUDA
%
% cudaVecCopyFromHost(vec, cuVec) copies the whole array from vec in MATLAB
% to cuVec on CUDA.
%
% SEE ALSO cudaVecCopySubFromHost, cudaVecCopySubToHost

%#codegen -args {m2c_vec, CudaVec}

[errCode, toplevel] = cudaVecCopySubFromHost(int32(length(vec)), vec, ...
    int32(1), int32(1), cuVec, int32(1));
