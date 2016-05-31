function [mat, errCode] = cuBCRSCreate_3args(mb, nb, nnzb)
[mat, errCode] = cuBCRSCreate(mb, nb, nnzb);
