function hdl = CusparseHandle(varargin) %#codegen
%Map an opaque object into a cusparseHandle_t object
%
%  CusparseHandle() simply returns a definition of the m2c_opaque_type,
%  suitable in the argument specification for codegen.
%
%  CusparseHandle(obj) or CusparseHandle(false) converts obj into a cusparseHandle_t object.
%
%  CusparseHandle(obj, 'wrap') or CusparseHandle(obj, true) wraps the obj into an opaque 
%  object. This should be used if the opaque object needs to be returned to MATLAB.

coder.inline('always');

hdl = m2c_opaque_obj('cusparseHandle_t', varargin{:});
