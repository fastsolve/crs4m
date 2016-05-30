function val = mspZero(type)
% Obtains zero of a particular type

coder.inline('always');

if type == MSP_DOUBLE
    val = 0;
elseif type == MSP_SINGLE
    val = single(0);
elseif type == MSP_INT32
    val = int32(0);
elseif type == MSP_UINT32
    val = uint32(0);
elseif type == MSP_INT8
    val = int8(0);
elseif type == MSP_UINT8
    val = uint8(0);
elseif type == MSP_INT16
    val = int16(0);
elseif type == MSP_UINT16
    val = uint16(0);
elseif type == MSP_INT64
    val = int64(0);
elseif type == MSP_UINT64
    val = uint64(0);
elseif type == MSP_COMPLEX
    val = complex(single(0));
elseif type == MSP_DOUBLE_COMPLEX
    val = complex(0);
else
    m2c_error('mspZero:WrongType', 'Unknow data type %d.\n', type)
end
