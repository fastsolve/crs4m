function [cuVec, errCode] = cuVecCopyToGPU(vec, cuVec)
% Copies a vector from MATLAB to CUDA
%
% cuVecCopyToGPU(vec, cuVec) copies the whole array from vec in MATLAB
% to cuVec on CUDA.
%
% [cuVec, errCode] = cuVecCopyToGPU(vec) copies the whole array from
% vec in MATLAB to a new vector cuVec on CUDA. This is a shorthand for
% calling cuVecCreate followed by calling cuVecCopyToGPU. The vector
% cuVec must be destroied by calling cuVecDestroy after use.
%
% SEE ALSO cuVecCopyFromGPU, cuVecCopySubToGPU, cuVecDestroy

if size(vec,2)~=1 && (isempty(coder.target) || m2c_debug)
    m2c_warn('cuVecCopyToGPU:NotVector', 'First input should be a column vector.\n');
end

if isreal(vec)
    if isa(vec, 'double')
        type = CU_DOUBLE;
    else
        type = CU_SINGLE;
    end
elseif isa(vec, 'double')
    type = CU_DOUBLE_COMPLEX;
else
    type = CU_COMPLEX;
end

if nargin==1
    cuVec = cuVecCreate(int32(size(vec,1)), type);
end

errCode = cuVecCopySubToGPU(int32(length(vec)), vec, ...
    int32(1), int32(1), cuVec, int32(1));

function test %#ok<DEFNU>
%!test
%! m=int32(1000000);
%! x =rand(m,1);
%! cuda_x = cuVecCopyToGPU(x);
%! y = cuVecCopyFromGPU(cuda_x);
%! cuVecDestroy(cuda_x);
%! assert(isequal(x, y));
