function [errCode, toplevel] = cuBlasDestroy(hdl)
%Destroys the handle to the CUBLAS library context.
%
%  errCode = cuBlasDestroy(hdl)
%
%  SEE ALSO: cuBlasCreate
%
% cuBLAS C interface:
%   cublasStatus_t cublasDestroy(cublasHandle_t handle)

%#codegen -args {CublasHandle}

coder.cinclude('mspack.h');

errCode = int32(-1);

if ~isempty(coder.target)
    errCode = coder.ceval('cublasDestroy', CublasHandle(hdl));
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('CUDA:RuntimeError', 'cublasDestroy returned error code %s\n', cuBlasGetErrorCode(errCode));
    end
end
end
