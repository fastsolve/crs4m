function cstr = cuSparseGetErrorString(errCode) %#codegen

switch errCode
    case CUSPARSE_STATUS_SUCCESS
        cstr = ['CUSPARSE_STATUS_SUCCESS' char(0)];
    case CUSPARSE_STATUS_NOT_INITIALIZED
        cstr = ['CUSPARSE_STATUS_NOT_INITIALIZED' char(0)];
    case CUSPARSE_STATUS_ALLOC_FAILED
        cstr = ['CUSPARSE_STATUS_ALLOC_FAILED' char(0)];
    case CUSPARSE_STATUS_INVALID_VALUE
        cstr = ['CUSPARSE_STATUS_INVALID_VALUE' char(0)];
    case CUSPARSE_STATUS_ARCH_MISMATCH
        cstr = ['CUSPARSE_STATUS_ARCH_MISMATCH' char(0)];
    case CUSPARSE_STATUS_MAPPING_ERROR
        cstr = ['CUSPARSE_STATUS_MAPPING_ERROR' char(0)];
    case CUSPARSE_STATUS_EXECUTION_FAILED
        cstr = ['CUSPARSE_STATUS_EXECUTION_FAILED' char(0)];
    case CUSPARSE_STATUS_INTERNAL_ERROR
        cstr = ['CUSPARSE_STATUS_INTERNAL_ERROR' char(0)];
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED
        cstr = ['CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED' char(0)];
    otherwise
        cstr = ['Unknown error' char(0)];
end
