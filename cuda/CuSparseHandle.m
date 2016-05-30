function hdl = CuSparseHandle(varargin) %#codegen
%Map an opaque object into a cusparseHandle_t object
%
%  CuSparseHandle() simply returns a definition of the m2c_opaque_type,
%  suitable in the argument specification for codegen.
%
%  CuSparseHandle(obj) or CuSparseHandle(false) converts obj into a cusparseHandle_t object.
%
%  CuSparseHandle(obj, 'wrap') or CuSparseHandle(obj, true) wraps the obj into an opaque 
%  object. This should be used if the opaque object needs to be returned to MATLAB.

coder.inline('always');

hdl = m2c_opaque_obj('cusparseHandle_t', varargin{:});
