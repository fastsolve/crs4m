function type = mspType(mType, isreal)
% Convert from a mType string to MSPACK type.

coder.inline('always');

switch mType
    case 'double'
        if isreal
            type = MSP_DOUBLE;
        else
            type = MSP_DOUBLE_COMPLEX;
        end
    case 'single'
        if isreal
            type = MSP_SINGLE;
        else
            type = MSP_COMPLEX;
        end
    case 'int8'
        type = MSP_INT8;
    case 'uint8'
        type = MSP_UINT8;
    case 'int16'
        type = MSP_INT16;
    case 'uint16'
        type = MSP_UINT16;
    case 'int32'
        type = MSP_INT32;
    case 'uint32'
        type = MSP_UINT32;
    case 'int64'
        type = MSP_INT64;
    case 'uint64'
        type = MSP_UINT64;
    otherwise
        m2c_error('mspZero:WrongType', ['Unknow supported data type ' mType '.\n']);
end
end
