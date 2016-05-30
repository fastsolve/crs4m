function type = cuType(mType, isreal)
% Convert from a mType string to cu type.

coder.inline('always');

switch mType
    case 'double'
        if isreal
            type = CU_DOUBLE;
        else
            type = CU_DOUBLE_COMPLEX;
        end
    case 'single'
        if isreal
            type = CU_SINGLE;
        else
            type = CU_COMPLEX;
        end
    case 'int8'
        type = CU_INT8;
    case 'uint8'
        type = CU_UINT8;
    case 'int16'
        type = CU_INT16;
    case 'uint16'
        type = CU_UINT16;
    case 'int32'
        type = CU_INT32;
    case 'uint32'
        type = CU_UINT32;
    case 'int64'
        type = CU_INT64;
    case 'uint64'
        type = CU_UINT64;
    otherwise
        m2c_error('cuZero:WrongType', ['Unknow supported data type ' mType '.\n']);
end
end
