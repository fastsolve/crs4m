function prod = vecDot_acc(u, v, prod, nthreads, mode) %#codegen
if isequal(mode, 'blas')
    prod = vecDot(u, v, prod, nthreads, 'blas');
elseif isequal(mode, 'acc')
    prod = vecDot(u, v, prod, nthreads, 'acc');
else
    m2c_warn('vecDot:WrongMode', ...
        'Unsupported supported fast mode.\n');
    prod = vecDot(u, v, prod, nthreads);
end
