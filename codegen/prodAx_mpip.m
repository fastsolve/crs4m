function b = prodAx_mpip( A, x, b, nthreads, comm, pbmsg) %#codegen

b = prodAx( A, x, b, nthreads, MPI_Comm(comm), pbmsg);
