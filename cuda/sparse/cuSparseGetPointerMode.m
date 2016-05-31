function [mode, errCode, toplevel] = cuSparseGetPointerMode(hdl)
%Obtains the pointer mode used by the cuSPARSE library.
%
%  [mode, errCode] = cuSparseGetPointerMode(hdl)
%  where mode is either CUSPARSE_POINTER_MODE_HOST or CUSPARSE_POINTER_MODE_DEVICE
%
%  SEE ALSO: cuSparseSetPointerMode
%
% cuSPARSE C interface:
%   cusparseStatus_t cuSparseGetPointerMode(cusparseHandle_t handle, cusparsePointerMode_t *mode))

%#codegen -args {CuSparseHandle}

coder.cinclude('mspack.h');

errCode = int32(-1);
if ~isempty(coder.target)
    t_mode = coder.opaque('cusparsePointerMode_t');
    errCode = coder.ceval('cusparseGetPointerMode', ...
        CuSparseHandle(hdl), coder.wref(t_mode));
    mode = int32(0); %#ok<NASGU>
    mode = coder.ceval(' ', t_mode);
    toplevel = nargout>2;
    
    if errCode && (toplevel || m2c_debug)
        m2c_error('CUDA:RuntimeError', 'cusparseGetPointerMode returned error code %s\n', cuSparseGetErrorString(errCode));
    end
end
end
