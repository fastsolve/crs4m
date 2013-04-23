function [row_ptr, col_ind, val] = crs_transpose(row_ptr, col_ind, val)
% Transpose a give matrix in compressed-row format {row_ptr, col_ind, val}.
% It is equivalent to constructing a compressed-column format of the matrix
% by treating row_ptr as col_ptr and treating col_ind as row_ind in the
% output.

%#codegen -args {coder.typeof(int32(0),[inf,1]),coder.typeof(int32(0),[inf,1]),
%#codegen coder.typeof(0,[inf,1])}

At_rowind = col_ind;
At_colind = crs_rowind(row_ptr, col_ind);

% Exchange row and col_ind.
[row_ptr, col_ind, val] = crs_create( At_rowind, At_colind, val);
