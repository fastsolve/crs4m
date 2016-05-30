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
    case 'CUSPARSE_DIAG_TYPE_NON_UNIT'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_DIAG_TYPE_NON_UNIT'));
    case 'CUSPARSE_DIAG_TYPE_UNIT'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_DIAG_TYPE_UNIT'));
    case 'CUSPARSE_FILL_MODE_LOWER'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_FILL_MODE_LOWER'));
    case 'CUSPARSE_FILL_MODE_UPPER'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_FILL_MODE_UPPER'));
    case 'CUSPARSE_INDEX_BASE_ZERO'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_INDEX_BASE_ZERO'));
    case 'CUSPARSE_INDEX_BASE_ONE'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_INDEX_BASE_ONE'));
    case 'CUSPARSE_MATRIX_TYPE_GENERAL'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_MATRIX_TYPE_GENERAL'));
    case 'CUSPARSE_MATRIX_TYPE_SYMMETRIC'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_MATRIX_TYPE_SYMMETRIC'));
    case 'CUSPARSE_MATRIX_TYPE_HERMITIAN'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_MATRIX_TYPE_HERMITIAN'));
    case 'CUSPARSE_MATRIX_TYPE_TRIANGULAR'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_MATRIX_TYPE_TRIANGULAR'));
    case 'CUSPARSE_DIRECTION_ROW'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_DIRECTION_ROW'));
    case 'CUSPARSE_DIRECTION_COLUMN'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_DIRECTION_COLUMN'));
    case 'CUSPARSE_OPERATION_NON_TRANSPOSE'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_OPERATION_NON_TRANSPOSE'));
    case 'CUSPARSE_OPERATION_TRANSPOSE'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_OPERATION_TRANSPOSE'));
    case 'CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE'
        val = coder.ceval(' ', coder.opaque('int',  'CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE'));
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
