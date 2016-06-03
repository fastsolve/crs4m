function mat = bsrCreateFromAIJ(rows, cols, vs, nrows, ncols) %#codegen
%Creates a sparse BSR matrix from AIJ (IJV) format.
%
%    mat = bsrCreateFromSparse(is, js, vs, nrows, ncols)
%
% See also bsrCreate, bsrCreateFromSparse

%#codegen -args {m2c_intvec, m2c_intvec, m2c_vec, int32(0), int32(0)}

t = length(vs)/length(rows);
if t>1
    blkdim = int32(sqrt(t));
else
    blkdim = int32(1);
end

nnz = int32(length(cols));
mat = bsrCreate(int32(nrows), int32(ncols), nnz, blkdim, class(vs));

%% Construct mat.rowptr
for i=1:int32(length(rows))
    mat.rowptr(rows(i)+1) = mat.rowptr(rows(i)+1) + 1;
end

mat.rowptr(1) = 1;
for i=1:nrows
    mat.rowptr(i+1) = mat.rowptr(i) + mat.rowptr(i+1);
end

%% Construct mat.colind and mat.vals
% Check whether row indices are in ascending order
ascend = true;
for i=2:nnz
    if rows(i)<rows(i-1)
        ascend = false;
        break;
    end
end

indx = coder.nullcopy(zeros(nnz, 1, 'int32'));
if ascend
    % if rows are already in sorted, simply copy cols to mat.colind
    mat.colind(:) = int32(cols);
    for i=1:nnz
        indx(i) = i;
    end
else
    for i=1:nnz
        j = mat.rowptr(rows(i));
        indx(j) = i;
        mat.colind(j) = cols(i);
        mat.rowptr(rows(i)) = mat.rowptr(rows(i)) + 1;
    end
    
    % Recover mat.rowptr
    for i=length(mat.rowptr):-1:2
        mat.rowptr(i) = mat.rowptr(i-1);
    end
    mat.rowptr(1)=1;
end

%% Sort colind
[mat.colind, indx] = sortCols(mat.rowptr, mat.colind, indx);

%% Copy values from vs to mat.vals
if blkdim==1
    % This is simply CRS format
    for i=1:nnz
        mat.vals(i) = vs(indx(i));
    end
else
    blklen = blkdim*blkdim;
    for i=1:nnz
        istart_l = blklen*(i-1)+1;
        istart_r = blklen*(indx(i)-1)+1;
        for j=1:blklen
            mat.vals(istart_l + j) = vs(istart_r + j);
        end
    end
end
