function s = sqnormf_mpi(A, s, nthreads, comm) %#codegen

s = sqnormf(A, s, nthreads, MPI_Comm(comm));
