function output = CudaVec(varargin) %#codegen
%Createa an opaque vector object on CUDA.
%
%  CudaVec() simply returns a definition of a structure,
%  suitable in the argument specification for codegen.
%
%  CudaVec(obj, type) converts obj into a pointer of the given type, where
%  type is a constant string containing a C pointer (such as 'double *').
%
%  CudaVec(obj, type, 'offset', i) returns a
%          into a CudaVec object.

%  CudaVec(obj, type, n, 'wrap') wraps the given opaque pointer
%          into a CudaVec object.

coder.inline('always');

narginchk(0, 4);

if nargin==0
    output = coder.typeof(struct('data', coder.typeof(uint32(0),[2,1]), ...
        'type', int32(0), 'len', int32(0)));
    return;
elseif nargin==4 && ischar(varargin{4})
    output = struct('data', zeros(2,1,'uint32'), ...
        'type', int32(0), 'len', int32(0));
    
    obj = varargin{1};
    output.type = int32(varargin{2});
    output.len = int32(varargin{3});
    
    if ~isempty(obj)
        ptr = coder.opaque('uint32_T *', 'NULL');
        ptr = coder.ceval('(uint32_T *)', coder.rref(obj));
        for i=int32(1):2
            output.data(i) = coder.ceval('*', ptr);
            ptr = m2c_offset_ptr(ptr, int32(1));
        end
    end
elseif nargin==2
    if isstruct(varargin{1})
        output = castdata(varargin{2}, varargin{1}.data);
    else
        output = varargin{1};
    end
elseif nargin==4 && ischar(varargin{3})
    % Compute pointer with offset
    if isstruct(varargin{1})
        output = m2c_offset_ptr(castdata(varargin{2}, varargin{1}.data), varargin{4});
    else
        output = m2c_offset_ptr(castdata(varargin{2}, varargin{1}), varargin{4});
    end
else
    % Undefined.
end
