function y = crs_prodSymmatVec(row_ptr, col_ind, val, x, y)
% Compute matrix-vector multiplication y=A*x for a symmetric matrix A,
% whose upper triangular part is stored in compressed row format.
% y = crs_matSymVecProd(row_ptr, col_ind, val, x, y)

% See http://www.netlib.org/linalg/html_templates/node98.html

% %#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
% %#codegen coder.typeof(0, [inf,1]), coder.typeof(0, [inf,1]), coder.typeof(0, [inf,1])}

if nargin<5; y=nullcopy(zeros(size(x))); end

for i = 1:int32(length(y))
    y(i)  = 0;
    for j = row_ptr(i) : row_ptr(i+1) - 1
        y(i) = y(i) + val(j) * x(col_ind(j));
    end
end

for j = 1 : int32(length(x))
    if row_ptr(j)<row_ptr(j+1) && col_ind(row_ptr(j))==j
        % Proper CRS format with nonzero diagonal
        for i = row_ptr(j)+1 : row_ptr(j+1)-1
            y(col_ind(i)) = y(col_ind(i)) + val(i) * x(j);
        end
    else
        for i = row_ptr(j) : row_ptr(j+1)-1
            if col_ind(i)>j
                y(col_ind(i)) = y(col_ind(i)) + val(i) * x(j);
            end
        end
    end
end

function test  %#ok<DEFNU>
%!test
%! for k=1:100
%!     A = sprand(20,20,0.5); A=A+A'; x = rand(20,1);
%!     [row_ptr, col_ind, val] = crs_createFromSparse(triu(A),20);
%!     y = crs_prodSymmatVec( row_ptr, col_ind, val, x);
%!     assert( norm(A*x-y)<=1.e-12);
%! end
