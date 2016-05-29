function prod = vecDot_cublas(u, v, prod, nthreads, mode, cublasHdl) %#codegen
if ~isequal(mode, 'cublas')
    m2c_warn('vecDot:WrongMode', 'Incorrect mode. Assuming ''cublas''.\n');
end

prod = vecDot(u, v, prod, nthreads, 'cublas', cublasHdl);
