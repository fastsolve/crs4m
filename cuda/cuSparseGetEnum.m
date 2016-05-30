function val = cuSparseGetEnum(str) %#codegen
% Obtains the value of an enumerate number in cuSPARSE
%
%  val = cuSparseGetEnum(str)
%
% SEE ALSO: cuSparseGetErrorString

%#codegen -args {m2c_string}

coder.cinclude('mspack.h');

val = int32(0);
switch str
    case 'CUSPARSE_POINTER_MODE_HOST'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_POINTER_MODE_HOST'));
    case 'CUSPARSE_POINTER_MODE_DEVICE'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_POINTER_MODE_DEVICE'));
    case 'CUSPARSE_STATUS_SUCCESS'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_SUCCESS'));
    case 'CUSPARSE_STATUS_NOT_INITIALIZED'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_NOT_INITIALIZED'));
    case 'CUSPARSE_STATUS_ALLOC_FAILED'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_ALLOC_FAILED'));
    case 'CUSPARSE_STATUS_INVALID_VALUE'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_INVALID_VALUE'));
    case 'CUSPARSE_STATUS_ARCH_MISMATCH'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_ARCH_MISMATCH'));
    case 'CUSPARSE_STATUS_MAPPING_ERROR'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_MAPPING_ERROR'));
    case 'CUSPARSE_STATUS_EXECUTION_FAILED'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_EXECUTION_FAILED'));
    case 'CUSPARSE_STATUS_INTERNAL_ERROR'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_INTERNAL_ERROR'));
    case 'CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED'));
    otherwise
        val = int32(-1);
end
