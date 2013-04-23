function A = crs_create( is, js, vs, varargin)
%crs_create  Create a sparse matrix in CRS-format from ijv format
%
%    [A.row_ptr, A.col_ind, A.val, ni, nj] = crs_create( is, js, vs [, ni , nj]);
%    A = crs_create( is, js, vs [, ni, nj]);
% In the second case, A is a struct with fields A.row_ptr, A.col_ind, A.val, ...
%    nrows, and ncols.
%
% See also crs_matrix
%
% Note: This function does not use multithreading. If there is no input,
%   this function creates a type declaration to be used with codegen.

%#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
%#codegen coder.typeof(0, [inf,1])}
%#codegen crs_create1 -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
%#codegen coder.typeof(0, [inf,1]), int32(0), int32(0)}

if nargin<4; nrows = max(is);
else nrows = int32(varargin{1}); end
if nargin<5; ncols = max(js);
else ncols = int32(varargin{2}); end

A = struct( 'row_ptr', nullcopy(zeros(nrows+1,1,'int32')), ...
    'col_ind', nullcopy(zeros(size(js),'int32')), ...
    'val', nullcopy(zeros(size(js), class(vs))), ...
    'nrows', nrows, 'ncols', ncols);

%% Construct A.row_ptr
for i=1:int32(length(is))
    A.row_ptr(is(i)+1) = A.row_ptr(is(i)+1) + 1;
end

A.row_ptr(1) = 1;
for i=1:nrows
    A.row_ptr(i+1) = A.row_ptr(i) + A.row_ptr(i+1);
end

%% Construct A.col_ind and A.val
% Check whether row indices are in ascending order
ascend = true;
for i=2:length(is)
    if is(i)<is(i-1)
        ascend = false;
        break;
    end
end

if ascend
    % if IS is already in ascending order, simply return js as A.col_ind
    % and vs as A.val.
    A.col_ind = int32(js);
    A.val = vs;
else
    % Construct A.col_ind and A.val
    A.col_ind = nullcopy(zeros(length(js),1,'int32'));
    A.val = nullcopy(zeros(length(js),1, class(vs)));
    
    for i=1:length(is)
        j = A.row_ptr(is(i));
        A.val(j) = vs(i);
        A.col_ind(j) = js(i);
        A.row_ptr(is(i)) = A.row_ptr(is(i)) + 1;
    end
    
    % Recover A.row_ptr
    for i=length(A.row_ptr):-1:2
        A.row_ptr(i) = A.row_ptr(i-1);
    end
    A.row_ptr(1)=1;
end
