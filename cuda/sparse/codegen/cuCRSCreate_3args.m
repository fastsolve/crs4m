function [mat, errCode] = cuCRSCreate_3args(m, n, nnz)
[mat, errCode] = cuCRSCreate(m, n, nnz);
