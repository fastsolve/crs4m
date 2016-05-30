function [errCode, toplevel] = cuSparseSetStream(hdl, strm)
%Sets the cudaStream used by the cuSPARSE library.
%
%  errCode = cuSparseSetStream(hdl, strm) sets the stream to strm.
%
%  SEE ALSO: cuSparseGetStream
%
% cuSPARSE C interface:
%   cusparseStatus_t cuSparseSetStream(cusparseHandle_t handle, cudaStream_t strm))

%#codegen -args {CuSparseHandle, CuStreamHandle}

coder.cinclude('mspack.h');

errCode = int32(-1);
if ~isempty(coder.target)
    errCode = coder.ceval('cusparseSetStream', CuSparseHandle(hdl), CuStreamHandle(strm));
    
    toplevel = nargout>1;
    if errCode && (toplevel || m2c_debug)
        m2c_error('cuSparse:RuntimeError', 'cusparseSetStream returned error code %s\n', ...
            cuSparseGetErrorString(errCode))
    end
end
end
