function output = MatBSR(varargin) %#codegen
%Createa a sparse matrix object in block compressed-row sparse (BSR) format.
%
%  MatBSR() returns a definition of a structure, suitable in the
%  argument specification for codegen. (Note: Any compiled code
%  with BSR as an input variable in its top-level would support
%  only double precision.)
%
%  MatBSR(rowptr, colind, vals, [mb], [nb], [nnzb], [blkdim]) wraps
%  the given data into sparse mb-by-nb block-matrix in BSR format. Each
%  square block is of size blkdim-by-blkdim, stored in column major.
%  If mb is not given, it is set to length(rowptr)-1.
%  If nb is not given, it is assumed to be the same as mb.
%  If nnzb is not given, it is set to length(colind).
%  If blkdim is not given, then it is assumed to be 1.
%
%  MatBSR(obj, 'rowptr') or MatBSR(obj, 'colind') returns an opaque C pointer
%  to the field rowptr or colind, respectively. This is useful when
%  passing a field of obj to a C function in read-only mode, because MATLAB
%  coder does not allow specifying coder.rref on a field of a struct.
%
%  MatBSR(obj, 'vals', ctype) returns an opaque C pointer to the field vals
%  of a particular pointer type specified by type, such as 'double *'.
%
% See also Mat, Vec, CuBSR

coder.inline('always');

narginchk(0, 7);

if nargin==0
    output = coder.typeof(struct(...
        'rowptr', m2c_intvec, 'colind', m2c_intvec, ...
        'vals', m2c_vec, 'type', int32(0), ...
        'dimsb', m2c_dims(2), 'nnzb', int32(0), ...
        'blkdim', int32(0)));
    coder.cstructname(output, 'MSP_BSR');
elseif nargin==2 && ischar(varargin{2}) && ~isequal(varargin{2}, 'vals')
    output = castptr('int *', varargin{1}.(varargin{2}));
elseif nargin==3 && ischar(varargin{2}) && isequal(varargin{2}, 'vals')
    output = castptr(varargin{3}, varargin{1}.(varargin{2}));
elseif nargin>=3
    output = struct(...
        'rowptr', int32(varargin{1}), 'colind', int32(varargin{2}), ...
        'vals', varargin{3}, 'type', int32(0), ...
        'dimsb', [int32(0), int32(0)], 'nnzb', int32(0), 'blkdim', int32(1));
    coder.cstructname(output, 'MSP_BSR');
    coder.varsize('output.rowptr', 'output.colind', 'output.vals', [inf,1]);
        
    if nargin>=4
        output.dimsb(1) = int32(varargin{4});
    else
        output.dimsb(1) = int32(length(output.rowptr))-1;
    end
    
    if nargin>=5
        output.dimsb(2) = int32(varargin{5});
    else
        output.dimsb(2) = output.dimsb(1);
    end
    
    if nargin>=6
        output.nnzb = int32(varargin{6});
    else
        output.nnzb = int32(length(output.colind));
    end
    
    if nargin>=7; 
        output.blkdim = int32(varargin{7});
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
