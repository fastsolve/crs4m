function sizepe = mspGetSizePerElement(type)
% Obtains the number of bypets for each data type.

coder.inline('always');

if type == MSP_DOUBLE || type == MSP_COMPLEX || type == MSP_INT64
    sizepe = int32(8);
elseif type == MSP_SINGLE || type == MSP_INT32 || type == MSP_INT32
    sizepe = int32(4);
elseif type == MSP_INT16 || type == MSP_INT16
    sizepe = int32(2);
elseif type == MSP_INT8 || type == MSP_INT8
    sizepe = int32(1);
elseif type == MSP_DOUBLE_COMPLEX
    sizepe = int32(16);
else
    m2c_error('mspGetSizePerElement:WrongType', 'Unknow data type %d.\n', type);
    sizepe = int32(0);
end
