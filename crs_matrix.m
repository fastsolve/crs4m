function varargout = crs_matrix( varargin) %#codegen
%crs_matrix  Create a sparse matrix in CRS-format from a matrix sp
%
%    A = crs_matrix( is, js, vs [, nrows, ncols]);
%    [row_ptr, col_ind, val, nrows, ncols] = crs_matrix( ...);
% In the first case, A is a struct with fields row_ptr, col_ind, 
% val, nrows, and ncols.
%
%    A = crs_matrix( sp, [, nrows, ncols]);
%    [row_ptr, col_ind, val, nrows, ncols] = crs_matrix( ...);
% This mode is incompatible with MATLAB Coder. It is provided
% for convenience within MATALB.
%
% See also crs_2sparse, crs_create
%
% Note: This function does not use multithreading. If there is no input,
%   this function creates a type declaration to be used with codegen.

if nargin==0
    varargout{1} = coder.typeof( struct( 'row_ptr', coder.typeof( int32(0), [inf,1]), ...
        'col_ind', coder.typeof( int32(0), [inf,1]), ...
        'val', coder.typeof( 0, [inf,1]), 'nrows', int32(0), 'ncols', int32(0)));
    return;
elseif issparse(varargin{1})
    [is,js,vs] = find(varargin{1});
    
    A = crs_create( int32(is), int32(js), vs, varargin{2:end});
else
    is = int32(varargin{1});
    js = int32(varargin{2});
    vs = varargin{3};
    
    A = crs_create( is, js, vs, varargin{4:end});
end

if nargout<=1
    varargout{1} = A;
else
    varargout{1} = A.row_ptr;
    varargout{2} = A.col_ind;
    if nargout>2; varargout{3} = A.val; end
    if nargout>3; varargout{4} = A.nrows; end
    if nargout>4; varargout{5} = A.ncols; end
end
