function [mat, errCode] = cuBCRSCreate_4args(mb, nb, nnzb, blkdim)
[mat, errCode] = cuBCRSCreate(mb, nb, nnzb, blkdim);
