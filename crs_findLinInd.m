function k = crs_findLinInd( row_ptr, col_ind, r, c) %#codegen
% Find linear index of an entry in the compressed-row format of a sparse matrix.

for i= row_ptr(r) : row_ptr(r+1)-1 
    if col_ind(i)==c
        k=i; return;
    end
end

k=int32(0);
error('Could not find (%d, %d) in matrix. The sparse matrix was not initialzied correctly.', r, c);
