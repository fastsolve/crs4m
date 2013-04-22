function y = crs_symMatVecProd(row_ptr, col_ind, val, x, y)
% Compute matrix-vector multiplication y=A*x for a symmetric sparse matrix A,
% where only the upper triangular part of A is stored in compressed row format
% {row_ptr, col_ind, val}.

% See http://www.netlib.org/linalg/html_templates/node98.html

%#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
%#codegen  coder.typeof(0, [inf,1]), coder.typeof(0, [inf,1]), coder.typeof(0, [inf,1])}

% Compute upper-triangular part.
for i = 1:length(y)
    y(i)  = 0;
    for j = row_ptr(i) : row_ptr(i+1) - 1
        y(i) = y(i) + val(j) * x(col_ind(j));
    end
end

% Compute lower-triangular part.
for j = 1 : int32(length(x))
    if row_ptr(j)~=row_ptr(j+1) && col_ind(row_ptr(j))==j
        % put diagonal entry first for better efficiency.
        for i = row_ptr(j)+1 : row_ptr(j+1)-1
            y(col_ind(i)) = y(col_ind(i)) + val(i) * x(j);
        end
    else
        for i = row_ptr(j) : row_ptr(j+1)-1
            if col_ind(i)~=j
                y(col_ind(i)) = y(col_ind(i)) + val(i) * x(j);
            end
        end
    end
end
