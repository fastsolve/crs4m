function [row_ptr, col_ind, val] = crs_triu(row_ptr, col_ind, val)
% Obtain upper-triangular part of matrix in compressed-row format.

% %#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
% %#codegen coder.typeof(0, [inf,1]}

assert(nargin==nargout && nargin>=2);

if nargin==2
    col_ind = crs_sort(row_ptr, col_ind);
else
    [col_ind, val] = crs_sort(row_ptr, col_ind, val);
end

offset = int32(0);
start = int32(1);

for i=1:int32(length(row_ptr))-1
    for j=start : row_ptr(i+1)-1
        if col_ind(j)<i
            offset = offset+1;
        elseif offset
            col_ind( j-offset) = col_ind( j);
            if nargin>2; val(j-offset) = val(j); end
        end
    end
    
    start = row_ptr(i+1);
    row_ptr(i+1) = row_ptr(i+1) - offset;
end

if offset
    newlen = int32(length(col_ind))-offset;
    col_ind = sub_colvec( col_ind, 1, newlen);
    if nargin>2; val = sub_colvec( val, 1, newlen); end
end

function test %#ok<DEFNU>
%!test
%! A = sprand(1000,1000,0.05);
%! B = A'*A;
%! C=triu(B);
%! [row_ptr, col_ind, val] = crs_createFromSparse(B);
%! [row_ptr1,col_ind1,val1]= crs_triu(row_ptr, col_ind, val);
%! assert(isequal(C, crs_createSparse(row_ptr1,col_ind1,val1)));
