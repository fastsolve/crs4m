function [colind, vals] = sortCols(rowptr, colind, vals)
%sortCols  Sort column indices within each row for CRS format.
%  [colind, vals] = sortCols(rowptr, colind, vals)

for i=1:int32(length(rowptr))-1
    ascend = true;
    for j=rowptr(i)+1 : rowptr(i+1)-1
        if colind(j)<colind(j-1)
            ascend = false;
            break;
        end
    end
    
    if ~ascend
        % sort in place
        [colind, vals] = heapsort_tag(colind, vals, ...
            rowptr(i+1)-rowptr(i), rowptr(i)-1);
    end
end
