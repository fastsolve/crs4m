function varargout = crs_matrix( varargin) %#codegen
%crs_matrix  Create a sparse matrix in CRS-format from a matrix sp
%
%    [row_ptr, col_ind, val, nrows, ncols] = crs_matrix( row, col, v [, nrows , ncols]);
%    A = crs_matrix( row, col, v [, nrows, ncols]);
% In the second case, A is a struct with fields row_ptr, col_ind, and val.
%
%    [row_ptr, col_ind, val, nrows, ncols] = crs_matrix( sp [, nrows [, ncols]]);
%    A = crs_matrix( sp, [, nrows [, ncols]]);
% This mode is incompatible with MATLAB Coder. It is provided
% for convenience within MATALB.
%
% See also crs_2sparse
%
% Note: This function does not use multithreading. If there is no input,
%   this function creates a type declaration to be used with codegen.

if nargin==0
    varargout{1} = coder.typeof( struct( 'row_ptr', coder.typeof( int32(0), [inf,1]), ...
        'col_ind', coder.typeof( int32(0), [inf,1]), ...
        'val', coder.typeof( 0, [inf,1]), 'nrows', int32(0), 'ncols', int32(0)));
    return;
elseif issparse(varargin{1})
    [row,col,v] = find(varargin{1});
    
    if nargin<2; nrows = int32(size(varargin{1},1));
    else nrows = int32(varargin{2}); end
    if nargin<3; ncols = int32(size(varargin{1},2));
    else ncols = int32(varargin{3}); end
else
    row = int32(varargin{1});
    col = int32(varargin{2});
    v = varargin{3};
    
    if nargin<4; nrows = max(row);
    else nrows = int32(varargin{4}); end
    if nargin<5; ncols = max(col);
    else ncols = int32(varargin{4}); end
end

%% Construct row_ptr
row_ptr = zeros(nrows+1,1,'int32');

for i=1:int32(length(row))
    row_ptr(row(i)+1) = row_ptr(row(i)+1) + 1;
end

row_ptr(1) = 1;
for i=1:nrows
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
    val = nullcopy(zeros(length(col),1, class(v)));
    
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

if nargin==2 && m+1>length(row_ptr)
    row_ptr = [row_ptr; repmat( row_ptr(end),m+1-length(row_ptr),1)];
end

if nargout<=1
    varargout{1} = struct( 'row_ptr', row_ptr, 'col_ind', col_ind, ...
        'val', val, 'nrows', nrows, 'ncols', ncols);
else
    varargout{1} = row_ptr;
    varargout{2} = col_ind;
    varargout{3} = val;
end
