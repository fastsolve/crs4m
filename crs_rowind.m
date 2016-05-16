function row_ind = crs_rowind(row_ptr, col_ind, varargin) %#codegen
% Obtain the row indices from the compressed row format.
%
%     row_ind = crs_rowind(row_ptr, col_ind, [usethreads])
% If the third argument is present, it will utilize parfor
% if a parallel section has already been created.
%
% See also crs_create

%#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0),[inf,1])}

row_ind = nullcopy( zeros(size(col_ind),class(col_ind)));

nrows = int32(length(row_ptr))-1;
if nargin>2
    i=int32(0); j=int32(0); [i,j] = m2c_touch(i,j);
    MACC_begin_for; MACC_clause_private(i, j)
end
for i=1:nrows
    for j = row_ptr(i) : row_ptr(i+1) - 1
        row_ind(j) = i;
    end
end
if nargin>2; MACC_end_for; end
