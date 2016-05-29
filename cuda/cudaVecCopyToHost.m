function [vec, errCode, toplevel] = cudaVecCopyToHost(cuVec, vec)
% Copies a vector from CUDA to MATLAB.
%
% [vec, errCode] = cudaVecCopyToHost(cuVec, vec) copies the whole array
% from cuVec on CUDA device to vec in MATLAB. The length to be copied
% is equal to cuVec.len.
%
% SEE ALSO cudaVecCopySubToHost, cudaVecCopySubFromHost

%#codegen -args {CudaVec, m2c_vec}

[vec, errCode, toplevel] = cudaVecCopySubToHost(cuVec.len, cuVec, int32(1), ...
    vec, int32(1), int32(1));
