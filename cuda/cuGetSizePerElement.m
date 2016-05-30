function sizepe = cuGetSizePerElement(type)
% Obtains the number of bypets for each data type.

coder.inline('always');

if type == CU_DOUBLE || type == CU_COMPLEX || type == CU_INT64
    sizepe = int32(8);
elseif type == CU_SINGLE || type == CU_INT32 || type == CU_INT32
    sizepe = int32(4);
elseif type == CU_INT16 || type == CU_INT16
    sizepe = int32(2);
elseif type == CU_INT8 || type == CU_INT8
    sizepe = int32(1);
elseif type == CU_DOUBLE_COMPLEX
    sizepe = int32(16);
else
    m2c_error('cuGetSizePerElement:WrongType', 'Unknow data type %d.\n', type);
    sizepe = int32(0);
end
