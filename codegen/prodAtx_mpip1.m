function b = prodAtx_mpip1( A, x, b, nthreads, comm, pbmsg, pbsz) %#codegen

b = prodAtx( A, x, b, nthreads, MPI_Comm(comm), pbmsg, pbsz);
