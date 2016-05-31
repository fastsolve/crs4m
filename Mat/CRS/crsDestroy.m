function mat = crsDestroy(mat)
%Destroys a sparse matrix in CRS format.
%
%  mat = crsDestroy(mat) destroies the arrays pointed to by mat.
%
%  SEE ALSO: crsCreate

%#codegen -args {CRS}

coder.cinclude('mspack.h');

mat.vals(1:end) = [];
mat.colind(1:end) = [];
mat.rowptr(1:end) = [];
end
