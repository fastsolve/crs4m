function mat = crsCreate(m, n, nnz, type)
%Creates a sparse matrix in CRS format.
%
%  [mat, errCode] = crsCreate(m, n, nnz, [type]) creates a sparse
%  matrix of size m-by-n, where type is MSP_DOUBLE, MSP_SINGLE,
%  MSP_COMPLEX, MSP_DOUBLE_COMPLEX, MSP_INT*, or MSP_UINT*. If
%  not specified, then type is MSP_DOUBLE.
%
%  SEE ALSO: crsDestroy

%#codegen -args {int32(0), int32(0), int32(0), int32(0)}
%#codegen crsCreate_3args -args {int32(0), int32(0), int32(0)}

if nargin<4; type = MSP_DOUBLE; end

mat = CuCRS(cpder.nullcopy(zeros(m+1, 1, 'int32')), ...
    coder.nullcopy(zeros(nnz, 1, 'int32')), ...
    coder.nullcopy(zeros(nnz, 1, 'like', mspZero(type))), ...
    int32(type), int32(m), int32(n), int32(nnz));
end
