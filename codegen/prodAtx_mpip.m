function b = prodAtx_mpip( A, x, b, nthreads, comm, pbmsg) %#codegen

b = prodAtx( A, x, b, nthreads, MPI_Comm(comm), pbmsg);
