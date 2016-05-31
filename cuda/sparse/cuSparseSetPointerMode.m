function [errCode, toplevel] = cuSparseSetPointerMode(hdl, mode)
%Sets the pointer mode used by the cuSPARSE library.
%
%  errCode = cuSparseSetPointerMode(hdl, mode) sets the pointer mode to mode,
%  where mode is CUSPARSE_POINTER_MODE_HOST or CUSPARSE_POINTER_MODE_DEVICE.
%
%  SEE ALSO: cuSparseGetPointerMode
%
% cuSPARSE C interface:
%   cusparseStatus_t cuSparseSetPointerMode(cusparseHandle_t handle, cusparsePointerMode_t mode))

%#codegen -args {CuSparseHandle, m2c_int}

coder.cinclude('mspack.h');

errCode = int32(-1);
if ~isempty(coder.target)
    errCode = coder.ceval('cusparseSetPointerMode', CuSparseHandle(hdl), mode);
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('cuSparse:RuntimeError', 'cusparseSetPointerMode returned error code %s\n', ...
            cuSparseGetErrorString(errCode));
    end
end
end
