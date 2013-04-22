function row_ind = crs_obtainRowInd(row_ptr, col_ind) %#codegen
% Obtain the row indices from the compressed row format.
%
% See also crs_createFromAIJ

%#codegen -args {coder.typeof(int32(0),[inf,1]), coder.typeof(int32(0),[inf,1])}

row_ind = nullcopy( zeros(size(col_ind),'int32'));
for i=1:int32(length(row_ptr))-1
    for j = row_ptr(i) : row_ptr(i+1) - 1
        row_ind(j) = i;
    end
end
