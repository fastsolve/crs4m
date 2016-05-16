function [b, Ax] = prodAtAx_mpip1( A, x, b, Ax, nthreads, comm, pbmsg, pbsz) %#codegen

[b, Ax] = prodAtAx( A, x, b, Ax, nthreads, MPI_Comm(comm), pbmsg, pbsz);
