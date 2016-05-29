function vec = CublasHandle(varargin) %#codegen
%Map an opaque object into a cublasHandle_t object
%
%  CublasHandle() simply returns a definition of the m2c_opaque_type,
%  suitable in the argument specification for codegen.
%
%  CublasHandle(obj) or CublasHandle(false) converts obj into a cublasHandle_t object.
%
%  CublasHandle(obj, 'wrap') or CublasHandle(obj, true) wraps the obj into an opaque 
%  object. This should be used if the opaque object needs to be returned to MATLAB.

coder.inline('always');

vec = m2c_opaque_obj('cublasHandle_t', varargin{:});
