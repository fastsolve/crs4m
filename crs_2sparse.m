function sp = crs_2sparse(varargin)
%crs_2sparse Create a MATLAB sparse object from CRS format.
%      sp = crs_2sparse(row_ptr, col_ind, val [,m [,n]])
%      sp = crs_2sparse(A [,m [,n]])
%
% See also crs_create

% Note: This function is incompatible with MATLAB Coder.
% It is provided for convenience of testing.

assert( nargin>=1 && nargin<=5);

if isstruct( varargin{1})
    row_ptr = varargin{1}.row_ptr;
    col_ind = varargin{1}.col_ind;
    val = varargin{1}.val;
    s = 2;
else
    row_ptr = varargin{1};
    col_ind = varargin{2};
    val = varargin{3};
    s = 4;
end

row_ind = crs_obtainRowInd(row_ptr, col_ind);

sp = sparse(double(row_ind), double(col_ind), double(val), varargin{s:end});

%!test
%! for k=1:100
%!     sp = sprand(20,10,0.5); x = rand(10,1);
%!     A = crs_create(sp);
%!     sp1 = crs_2sparse(A, 20, 10);
%!     
%!     [row_ptr, col_ind, val] = crs_create(sp);
%!     sp1 = crs_2sparse(row_ptr, col_ind, val, 20, 10);
%!     
%!     assert( isequal( sp, sp1));
%! end
