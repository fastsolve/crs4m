function A = bsrCreateFromSparse(sp)
%Creates a sparse BSR matrix from MATLAB's built-in sparse matrix.
%
%    A = bsrCreateFromSparse(sp)
%
% Note: This function cannot be compiled.

[rows, cols, vs] = find(sp);

A = bsrCreateFromAIJ(rows, cols, vs, size(sp, 1), size(sp,2));
