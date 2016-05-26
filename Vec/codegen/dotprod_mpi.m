function prod = dotprod_mpi(u, v, prod, n, comm) %#codegen

prod = dotprod(u, v, prod, n, MPI_Comm(comm));
