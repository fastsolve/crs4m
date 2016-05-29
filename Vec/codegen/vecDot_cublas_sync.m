function [prod, buf] = vecDot_cublas_sync(u, v, buf, mode, cublasHdl, sync) %#codegen

if ~isequal(mode, 'cuda') 
    m2c_warn('vecDot:WrongInput', 'Wrong mode name. cuda is assumed.\n');
end

[prod, buf] = vecDot(u, v, buf, true, cublasHdl, sync);
