function [b, Atx] = prodAAtx_mpip(A, x, b, Atx, nthreads, comm, pbmsg) %#codegen

[b, Atx] = prodAAtx(A, x, b, Atx, nthreads, MPI_Comm(comm), pbmsg);
