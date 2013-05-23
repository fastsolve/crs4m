/*
 * crs_triu_mex.c
 *
 * Auxiliary code generation for function crs_triu
 *
 * C source code generated by m2c.
*
 */

#include "mex.h"

#define BUILD_MEX
/* Include the C file generated by codegen in lib mode */
#include "crs_triu.h"
#include "m2c.c"
#include "crs_triu.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void crs_triu_api(const mxArray ** prhs, const mxArray **plhs) {

    struct_T             A;
    struct_T             U;
    mxArray              *_sub_mx1;

    /* Marshall in function inputs */

    if ( !mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_triu:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if ( mxGetNumberOfFields( prhs[0])!=5)
        mexWarnMsgIdAndTxt("crs_triu:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if ( mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_triu:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField( prhs[0], 0, "row_ptr");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(m2cArray__common));
    alias_mxArray_to_m2cArray(_sub_mx1, (m2cArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "col_ind");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(m2cArray__common));
    alias_mxArray_to_m2cArray(_sub_mx1, (m2cArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "val");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu:WrongInputStruct",
            "Input argument A does not have the field val.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_triu:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(m2cArray__common));
    alias_mxArray_to_m2cArray(_sub_mx1, (m2cArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "nrows");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_triu:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField( prhs[0], 0, "ncols");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_triu:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);

    /* Preallocate output variables */
    *(void **)&U.row_ptr = mxCalloc(1, sizeof(m2cArray__common));    init_m2cArray( (m2cArray__common*)U.row_ptr, 1);
    *(void **)&U.col_ind = mxCalloc(1, sizeof(m2cArray__common));    init_m2cArray( (m2cArray__common*)U.col_ind, 1);
    *(void **)&U.val = mxCalloc(1, sizeof(m2cArray__common));    init_m2cArray( (m2cArray__common*)U.val, 1);

    /* Invoke the target function */
    crs_triu_initialize();
    crs_triu(&A, &U);
    crs_triu_terminate();
    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    m2cSize _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_m2cArray_to_mxArray((m2cArray__common*)U.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_m2cArray_to_mxArray((m2cArray__common*)U.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_m2cArray_to_mxArray((m2cArray__common*)U.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&U.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&U.ncols, mxINT32_CLASS));    /* Free temporary variables */
    free_m2cArray( (m2cArray__common*)A.val); mxFree( A.val);
    free_m2cArray( (m2cArray__common*)A.col_ind); mxFree( A.col_ind);
    free_m2cArray( (m2cArray__common*)A.row_ptr); mxFree( A.row_ptr);

    free_m2cArray( (m2cArray__common*)U.val); mxFree( U.val);
    free_m2cArray( (m2cArray__common*)U.col_ind); mxFree( U.col_ind);
    free_m2cArray( (m2cArray__common*)U.row_ptr); mxFree( U.row_ptr);
}

void crs_triu1_api(const mxArray ** prhs, const mxArray **plhs) {

    struct_T             A;
    struct_T             U;
    mxArray              *_sub_mx1;

    int32_T              k;

    /* Marshall in function inputs */

    if ( !mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if ( mxGetNumberOfFields( prhs[0])!=5)
        mexWarnMsgIdAndTxt("crs_triu1:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if ( mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_triu1:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField( prhs[0], 0, "row_ptr");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(m2cArray__common));
    alias_mxArray_to_m2cArray(_sub_mx1, (m2cArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "col_ind");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(m2cArray__common));
    alias_mxArray_to_m2cArray(_sub_mx1, (m2cArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "val");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputStruct",
            "Input argument A does not have the field val.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(m2cArray__common));
    alias_mxArray_to_m2cArray(_sub_mx1, (m2cArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "nrows");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_triu1:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField( prhs[0], 0, "ncols");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if ( mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_triu1:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);
    if ( mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_triu1:WrongInputType",
            "Input argument k has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("crs_triu1:WrongSizeOfInputArg",
            "Argument k should be a scalar.");
    k = *(int32_T*)mxGetData(prhs[1]);

    /* Preallocate output variables */
    *(void **)&U.row_ptr = mxCalloc(1, sizeof(m2cArray__common));    init_m2cArray( (m2cArray__common*)U.row_ptr, 1);
    *(void **)&U.col_ind = mxCalloc(1, sizeof(m2cArray__common));    init_m2cArray( (m2cArray__common*)U.col_ind, 1);
    *(void **)&U.val = mxCalloc(1, sizeof(m2cArray__common));    init_m2cArray( (m2cArray__common*)U.val, 1);

    /* Invoke the target function */
    crs_triu_initialize();
    crs_triu1(&A, k, &U);
    crs_triu_terminate();
    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    m2cSize _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_m2cArray_to_mxArray((m2cArray__common*)U.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_m2cArray_to_mxArray((m2cArray__common*)U.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_m2cArray_to_mxArray((m2cArray__common*)U.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&U.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&U.ncols, mxINT32_CLASS));    /* Free temporary variables */
    free_m2cArray( (m2cArray__common*)A.val); mxFree( A.val);
    free_m2cArray( (m2cArray__common*)A.col_ind); mxFree( A.col_ind);
    free_m2cArray( (m2cArray__common*)A.row_ptr); mxFree( A.row_ptr);

    free_m2cArray( (m2cArray__common*)U.val); mxFree( U.val);
    free_m2cArray( (m2cArray__common*)U.col_ind); mxFree( U.col_ind);
    free_m2cArray( (m2cArray__common*)U.row_ptr); mxFree( U.row_ptr);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 1) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_triu:TooManyOutputArguments","Too many output arguments for entry-point 'crs_triu'.");
        /* Call the API function. */
        crs_triu_api(prhs, (const mxArray**)plhs);
    }
    else if (nrhs == 2) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_triu1:TooManyOutputArguments","Too many output arguments for entry-point 'crs_triu1'.");
        /* Call the API function. */
        crs_triu1_api(prhs, (const mxArray**)plhs);
    }
    else
        mexErrMsgIdAndTxt("crs_triu:WrongNumberOfInputs","Incorrect number of input variables for entry-point 'crs_triu'.");
}