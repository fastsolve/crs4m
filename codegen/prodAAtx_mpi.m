function [b, Atx] = prodAAtx_mpi(A, x, b, Atx, nthreads, comm) %#codegen

[b, Atx] = prodAAtx(A, x, b, Atx, nthreads, MPI_Comm(comm));
