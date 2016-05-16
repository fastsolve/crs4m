function b = prodAx_mpip1( A, x, b, nthreads, comm, pbmsg, pbsz) %#codegen

b = prodAx( A, x, b, nthreads, MPI_Comm(comm), pbmsg, pbsz);
