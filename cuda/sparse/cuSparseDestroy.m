function [errCode, toplevel] = cuSparseDestroy(hdl)
%Destroys the cuSPARSE context specified by the handle.
%
%  errCode = cuSparseDestroy(hdl)
%
%  SEE ALSO: cuSparseCreate
%
% cuSPARSE C interface:
%   cusparseStatus_t cusparseDestroy(cusparseHandle_t handle)

%#codegen -args {CuSparseHandle}

coder.cinclude('mspack.h');

errCode = int32(-1);

if ~isempty(coder.target)
    errCode = coder.ceval('cusparseDestroy', CuSparseHandle(hdl));
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('CUDA:RuntimeError', 'cusparseDestroy returned error code %s\n', cuSparseGetErrorString(errCode));
    end
end
end
