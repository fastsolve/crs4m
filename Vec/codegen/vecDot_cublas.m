function [prod,buf] = vecDot_cublas(u, v, buf, mode, cublasHdl) %#codegen

if ~isequal(mode, 'cuda')
    m2c_warn('vecDot:WrongInput', 'Wrong mode name. cuda is assumed.\n');
end

[prod,buf] = vecDot(u, v, buf, true, cublasHdl);
