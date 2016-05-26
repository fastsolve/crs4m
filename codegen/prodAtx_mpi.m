function b = prodAtx_mpi(A, x, b, nthreads, comm) %#codegen

b = prodAtx(A, x, b, nthreads, MPI_Comm(comm));

end
