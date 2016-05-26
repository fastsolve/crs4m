function [b, Ax] = prodAtAx_mpi(A, x, b, Ax, nthreads, comm) %#codegen

[b, Ax] = prodAtAx(A, x, b, Ax, nthreads, MPI_Comm(comm));
