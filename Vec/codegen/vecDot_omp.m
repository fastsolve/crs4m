function [prod,buf] = vecDot_omp(u, v, buf, nthreads) %#codegen
[prod,buf] = vecDot(u, v, buf, nthreads);
