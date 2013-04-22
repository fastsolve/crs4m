function A = crs_createSparse(row_ptr, col_ind, val, varargin)
% Create a MATLAB sparse object from a CRS format.

% This function is not incompatible with MATLAB Coder. It is 
% provided for convenience of testing.

row_ind = crs_obtainRowInd(row_ptr, col_ind);

A = sparse(double(row_ind), double(col_ind), val, varargin{:});
