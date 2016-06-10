function output = CRS(varargin) %#codegen
%Createa a sparse matrix object in CRS format.
%
%  CRS() returns a definition of a structure, suitable in the
%  argument specification for codegen. (Note: Any compiled code
%  with CRS as an input variable in its top-level would support
%  only double precision.)
%
%  CRS(rowptr, colind, vals, [type], [m], [n], [nnz]) wraps the given
%  data into sparse m-by-n matrix in CRS format. If type is not given,
%  it inherits from vals. If m is not given, it is set to length(rowptr)-1.
%  If n is not given, it is assumed to be the same as m. If nnz is not given,
%  it is set to length(colind).
%
%  CRS(obj, 'rowptr') or CRS(obj, 'colind') returns an opaque C pointer
%  to the field rowptr or colind, respectively. This is necessary when
%  passing a field of obj to a C function in read-only mode, because MATLAB
%  coder does not allow specifying coder.rref on a field of a structure.
%
%  CRS(obj, 'vals', type) returns an opaque C pointer to the field vals
%  of a particular pointer type specified by type, such as 'double *'.
%
% See also BCRS, CuCRS

coder.inline('always');
narginchk(0, 7);

[output, ~] = BCRS(varargin{:});
