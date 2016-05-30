function output = CudaMat(varargin) %#codegen
%Createa an opaque dense matrix object on CUDA.
%
%  CudaMat() simply returns a definition of a structure, suitable in the
%  argument specification for codegen.
%
%  CudaMat(obj, type, m, n, 'wrap') wraps the given opaque pointer
%          into an m-by-n CudaMat object.
%
%  CudaMat(obj, type) converts obj into a CUDA pointer of the given type,
%  where type is a constant string, such as 'double *'.
%
%  CudaMat(obj, type, 'offset', k) returns a CUDA pointer but offset
%  by k elements.
%
% See also CudaVec,  CudaCRS, CudaBCRS, cudaMatCreate

coder.inline('always');

narginchk(0, 5);

if nargin==0
    output = coder.typeof(struct('data', coder.typeof(uint64(0)), ...
        'type', int32(0), 'dims', m2c_introw(2)));
    return;
elseif nargin==5 && ischar(varargin{5})
    output = struct('data', uint64(0), ...
        'type', int32(varargin{2}), ...
        'dims', [int32(varargin{3}), int32(varargin{4})]);
    
    obj = varargin{1};
    if ~isempty(obj)
        output.data = coder.ceval('*(uint64_T *)', coder.rref(obj));
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
