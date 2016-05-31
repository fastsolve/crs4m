function output = CRS(varargin) %#codegen
%Createa a sparse matrix object in CRS format.
%
%  CRS() simply returns a definition of a structure, suitable in the
%  argument specification for codegen. (Note: Any compiled code
%  with CRS as an input varaible in its top-level would support
%  only double precision.)
%
%  CRS(rowptr, colind, vals, [type], [m], [n], [nnz]) wraps the given
%  data into sparse m-by-n matrix in CRS format. If type is not given,
%  it inherits from vals. If m is not given, it is set to lenth(rowptr)-1.
%  If n is not given, it is assumed to be the same as m. If nnz is not given,
%  it is set to length(colind).
%
%  CRS(obj, 'rowptr') or CRS(obj, 'colind') returns an opaque C pointer
%  to the the field rowptr or colind, respectively. This is necessary when
%  passing a field of obj to a C finction in readonly mode, because MATLAB
%  coder does not allow specifying coder.rref on a field of a structure.
%
%  CRS(obj, 'vals', type) returns an opaque C pointer to the field vals
%  of a particular pointer type specified by type, such as 'double *'.
%
% See also CRS, CuVec, CudaBCRS, CRSCreate

coder.inline('always');

narginchk(0, 5);

if nargin==0
    output = coder.typeof(struct(...
        'rowptr', m2c_intvec, 'colind', m2c_intvec, ...
        'vals', m2c_vec, 'type', int32(0), ...
        'dims', m2c_dims(2), 'nnz', int32(0)));
    return;
elseif nargin==2 && ischar(varargin{2}) && ~isequal(varargin{2}, 'vals')
    output = castptr('int *', varargin{1}.(varargin{2}));
elseif nargin==3 && ischar(varargin{2}) && isequal(varargin{2}, 'vals')
    output = castptr(varargin{3}, varargin{1}.(varargin{2}));
elseif nargin>=3
    output = struct(...
        'rowptr', int32(varargin{1}), 'colind', int32(varargin{2}), ...
        'vals', varargin{3}, 'type', int32(0), ...
        'dims', [int32(0), int32(0)], ...
        'nnz', int32(varargin{7}));
    coder.varsize('output.rowptr', 'output.colind', 'output.vals', [inf,1]);
    
    if nargin>=4;
        output.type = int32(varargin{4});
    else
        output.type = cuType(class(varargin{4}), isreal(varargin{4}));
    end
    
    if nargin>=5;
        output.dims(1) = int32(varargin{5});
    else
        output.dims(1) = int32(length(output.rowptr))-1;
    end
    
    if nargin>=6;
        output.dims(2) = int32(varargin{6});
    else
        output.dims(2) = output.dims(1);
    end
    
    if nargin>=7;
        output.nnz = int32(varargin{7});
    else
        output.nnz = int32(length(output.colind));
    end
else
    % Undefined.
end

function ptr = castptr(type, data) %#codegen
%Casts the given data field into an opaque C pointer.
%    ptr = m2c_castptr(type, data)
coder.inline('always')

% Note: This might not work, since MATLAB Coder may make copy of the data
ptr = coder.opaque(type);
ptr = coder.ceval(['(' type ')'], coder.rref(data));
