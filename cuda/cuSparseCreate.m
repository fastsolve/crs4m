function [hdl, errCode, toplevel] = cuSparseCreate
%Initializes CUSPARSE and creates a handle to the CUSPARSE library context.
%
%  [hdl, errCode] = cuSparseCreate
%
%  SEE ALSO: cuSparseDestroy
%
% cuSPARSE C interface:
%   cusparseStatus_t cusparseCreate(cusparseHandle_t *handle)

%#codegen -args {}

coder.cinclude('mspack.h');
errCode = int32(-1);

if ~isempty(coder.target)
    t_hdl = coder.opaque('cusparseHandle_t');
    errCode = coder.ceval('cusparseCreate', coder.wref(t_hdl));
    
    toplevel = nargout>2;
    hdl = CuSparseHandle(t_hdl, toplevel);
    
    if errCode && (toplevel || m2c_debug)
        m2c_error('CUDA:RuntimeError', 'cusparseCreate returned error code %s\n', cuSparseGetErrorString(errCode));
    end
end
end
