function [errCode, toplevel] = cuBlasSetPointerMode(hdl, mode)
%Sets the pointer mode used by the cuBLAS library.
%
%  errCode = cuBlasGetPointerMode(hdl, mode) sets the pointer mode to mode,
%  where mode is CUBLAS_POINTER_MODE_HOST or CUBLAS_POINTER_MODE_DEVICE.
%
%  SEE ALSO: cuBlasGetPointerMode
%
% cuBLAS C interface:
%   cublasStatus_t cuBlasSetPointerMode(cublasHandle_t handle, cublasPointerMode_t mode))

%#codegen -args {CublasHandle, m2c_int}

coder.cinclude('mspack.h');

errCode = int32(-1);
if ~isempty(coder.target)
    errCode = coder.ceval('cublasSetPointerMode', CublasHandle(hdl), mode);
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('cuBlas:RuntimeError', 'cublasSetPointerMode returned error code %d\n', errCode)
    end
end
end
