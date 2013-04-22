function [row_ptr, col_ind, val] = crs_createFromSparse(A,m)
% Create a CRS format {row_ptr, col_ind, val} from a MATLAB sparse matrix A.
%
% See also crs_createFromAIJ

% This function is not incompatible with MATLAB Coder. It is 
% provided for convenience of testing.

[row,col,v] = find(A);

[row_ptr, col_ind, val] = crs_createFromAIJ(int32(row),int32(col),v,int32(size(A,1)));

if nargin>1 && m+1>length(row_ptr)
    row_ptr = [row_ptr; repmat( row_ptr(end),m+1-length(row_ptr),1)];
end
