function b = prodAx_mpi(A, x, b, nthreads, comm) %#codegen

b = prodAx(A, x, b, nthreads, MPI_Comm(comm));
