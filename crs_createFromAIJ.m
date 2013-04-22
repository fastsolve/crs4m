function [row_ptr, col_ind, val] = crs_createFromAIJ(row,col,v,nrow)
% Create a CRS format {row_ptr, col_ind, val} from AIJ format.
% Duplicate entries in the input will be preserved. To remove duplicate
% entries, call crs_uniqueColInd afterwards.
%
% See also crs_obtainRowInd, crs_uniqueColInd

%#codegen -args {coder.typeof(int32(0),[inf,1]), coder.typeof(int32(0),[inf,1])
%#codegen coder.typeof(0,[inf,1])}

if nargin<4;
    nrow = max(row);
end

%% Construct row_ptr
row_ptr = zeros(nrow+1,1,'int32');

for i=1:int32(length(row))
    row_ptr(row(i)+1) = row_ptr(row(i)+1) + 1;
end

row_ptr(1) = 1;
for i=1:nrow
    row_ptr(i+1) = row_ptr(i) + row_ptr(i+1);
end

%% Construct col_ind and val
% Check whether row indices are in ascending order
ascend = true;
for i=2:length(row)
    if row(i)<row(i-1)
        ascend = false;
        break;
    end
end

if ascend
    % if row is already in ascending order, simply return col as col_ind
    % and v as val.
    col_ind = int32(col);
    val = v;
else
    % Construct col_ind and val
    col_ind = nullcopy(zeros(length(col),1,'int32'));
    val = zeros(length(col),1);
   
    for i=1:length(row)
        j = row_ptr(row(i));
        val(j) = v(i);
        col_ind(j) = col(i);
        row_ptr(row(i)) = row_ptr(row(i)) + 1;
    end
    
    % Recover row_ptr
    for i=length(row_ptr):-1:2
        row_ptr(i) = row_ptr(i-1);
    end
    row_ptr(1)=1;
end

end
