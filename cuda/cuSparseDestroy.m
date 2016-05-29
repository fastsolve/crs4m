function [errCode, toplevel] = cuSparseDestroy(hdl)
%Destroys the handle to the CUSPARSE library context.
%
%  errCode = cuSparseDestroy(hdl)
%
%  SEE ALSO: cuSparseCreate
%
% cuSPARSE C interface:
%   cusparseStatus_t cusparseDestroy(cusparseHandle_t handle)

%#codegen -args {CusparseHandle}

coder.cinclude('mspack.h');

errCode = int32(-1);

if ~isempty(coder.target)
    errCode = coder.ceval('cusparseDestroy', CusparseHandle(hdl));
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('CUDA:RuntimeError', 'cusparseDestroy returned error code %s\n', cuSparseGetErrorCode(errCode));
    end
end
end
