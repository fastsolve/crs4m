function output = CuCRS(varargin) %#codegen
%Createa an opaque sparse matrix object in CRS format on CUDA.
%
%  CuCRS() simply returns a definition of a structure, suitable in the
%  argument specification for codegen.
%
%  CuCRS(rowptr, colind, vals, [type], [m], [n], [nnz]) wraps the given
%  opaque CUDA pointer into an m-by-n CuCRS object. If type is not given,
%  it inherits from vals. If m is not given, it is set to lenth(rowptr)-1.
%  If n is not given, it is assumed to be the same as m. If nnz is not given,
%  it is set to colind.len.
%
%  CuCRS(obj, 'rowptr') or CuCRS(obj, 'colind') converts the rowptr or
%  colind field into a CUDA pointer, respectively.
%
%  CuCRS(obj, 'vals', type) converts the vals field into a CUDA pointer
%  of a particular pointer type specified by type, such as 'double *'.
%
% See also CuBCRS, cuCRSCreate

coder.inline('always');

narginchk(0, 7);

[output, ~] = CuBCRS(varargin{:});

