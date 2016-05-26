function [b, Ax] = prodAtAx_mpip(A, x, b, Ax, nthreads, comm, pbmsg) %#codegen

[b, Ax] = prodAtAx(A, x, b, Ax, nthreads, MPI_Comm(comm), pbmsg);
