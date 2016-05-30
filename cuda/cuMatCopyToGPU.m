function [cuMat, errCode] = cuMatCopyToGPU(mat, cuMat)
% Copies a matrix from MATLAB to CUDA
%
% cuMatCopyToGPU(mat, cuMat) copies the whole array from mat in MATLAB
% to cuMat on CUDA.
%
% [cuMat, errCode] = cuMatCopyToGPU(mat) copies the whole array from
% mat in MATLAB to a new matrix cuMat on CUDA. This is a shorthand for
% calling cuMatCreate followed by calling cuMatCopyToGPU. The matrix
% cuMat must be destroied by calling cuMatDestroy after use.
%
% SEE ALSO cuMatCopyFromGPU, cuMatCopySubToGPU, cuMatDestroy

if size(mat,2)~=1 && (isempty(coder.target) || m2c_debug)
    m2c_warn('cuMatCopyToGPU:NotMattrix', 'First input should be a column matrix.\n');
end

if isreal(mat)
    if isa(mat, 'double')
        type = CU_DOUBLE;
    else
        type = CU_SINGLE;
    end
elseif isa(mat, 'double')
    type = CU_DOUBLE_COMPLEX;
else
    type = CU_COMPLEX;
end

if nargin==1
    cuMat = cuMatCreate(int32(size(mat,1)), int32(size(mat,2)), type);
end

errCode = cuMatCopySubToGPU(int32(size(mat,1)), int32(size(mat,2)), mat, cuMat);
