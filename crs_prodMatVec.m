function y = crs_prodMatVec(row_ptr, col_ind, val, x, y)
% Compute matrix-vector multiplication y=A*x for a sparse matrix A in
% compressed row format (row_ptr, col_ind, val).

% See http://www.netlib.org/linalg/html_templates/node98.html

%#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
%#codegen coder.typeof(0, [inf,1]), coder.typeof(0, [inf,1]), coder.typeof(0, [inf,1])}

if nargin<5;
    y=zeros(length(row_ptr)-1,1);
else
    for i = 1:int32(length(y)); y(i) =  0; end
end

for i = 1:int32(length(row_ptr))-1
    for j = row_ptr(i) : row_ptr(i+1) - 1
        y(i) = y(i) + val(j) * x(col_ind(j));
    end
end

function test  %#ok<DEFNU>
%!test
%! for k=1:100
%!     A = sprand(20,10,0.5); x = rand(10,1);
%!     [row_ptr, col_ind, val] = crs_createFromSparse(A);
%!     y = crs_prodMatVec( row_ptr, col_ind, val, x);
%!     assert( norm(A*x-y)<=1.e-12);
%! end
