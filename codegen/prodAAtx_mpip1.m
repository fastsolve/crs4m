function [b, Atx] = prodAAtx_mpip1(A, x, b, Atx, nthreads, comm, pbmsg, pbsz) %#codegen

[b, Atx] = prodAAtx(A, x, b, Atx, nthreads, MPI_Comm(comm), pbmsg, pbsz);
