/*
 * crs_tril_mex.c
 *
 * Auxiliary code generation for function crs_tril
 *
 * C source code generated by lib2mex.
*
 */

#include "mex.h"

#define BUILD_MEX
/* Include the C file generated by codegen in lib mode */
#include "crs_tril.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"


void crs_tril_api(const mxArray ** prhs, const mxArray **plhs) {


    struct_T             A;
    struct_T             L;
    mxArray              *_sub_mx1;


    /* Marshall in function inputs */

    if ( !mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_tril:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if ( mxGetNumberOfFields( prhs[0])!=5)
        mexWarnMsgIdAndTxt("crs_tril:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if ( mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_tril:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField( prhs[0], 0, "row_ptr");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "col_ind");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "val");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril:WrongInputStruct",
            "Input argument A does not have the field val.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_tril:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "nrows");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_tril:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField( prhs[0], 0, "ncols");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_tril:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);

    /* Preallocate output variables */
    *(void **)&L.row_ptr = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)L.row_ptr, 1);
    *(void **)&L.col_ind = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)L.col_ind, 1);
    *(void **)&L.val = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)L.val, 1);

    /* Invoke the target function */
    crs_tril_initialize();
    crs_tril(&A, &L);
    crs_tril_terminate();

    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    int32_T _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_emxArray_to_mxArray((emxArray__common*)L.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_emxArray_to_mxArray((emxArray__common*)L.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_emxArray_to_mxArray((emxArray__common*)L.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&L.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&L.ncols, mxINT32_CLASS));

    /* Free temporary variables */
    free_emxArray( (emxArray__common*)A.val); mxFree( A.val);
    free_emxArray( (emxArray__common*)A.col_ind); mxFree( A.col_ind);
    free_emxArray( (emxArray__common*)A.row_ptr); mxFree( A.row_ptr);

    free_emxArray( (emxArray__common*)L.val); mxFree( L.val);
    free_emxArray( (emxArray__common*)L.col_ind); mxFree( L.col_ind);
    free_emxArray( (emxArray__common*)L.row_ptr); mxFree( L.row_ptr);

}


void crs_tril1_api(const mxArray ** prhs, const mxArray **plhs) {


    struct_T             A;
    struct_T             L;
    mxArray              *_sub_mx1;

    int32_T              k;

    /* Marshall in function inputs */

    if ( !mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if ( mxGetNumberOfFields( prhs[0])!=5)
        mexWarnMsgIdAndTxt("crs_tril1:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if ( mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_tril1:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField( prhs[0], 0, "row_ptr");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "col_ind");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "val");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputStruct",
            "Input argument A does not have the field val.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "nrows");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_tril1:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField( prhs[0], 0, "ncols");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_tril1:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);
    if ( mxGetData(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_tril1:WrongInputType",
            "Input argument k has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("crs_tril1:WrongSizeOfInputArg",
            "Argument k should be a scalar.");
    k = *(int32_T*)mxGetData(prhs[1]);

    /* Preallocate output variables */
    *(void **)&L.row_ptr = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)L.row_ptr, 1);
    *(void **)&L.col_ind = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)L.col_ind, 1);
    *(void **)&L.val = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)L.val, 1);

    /* Invoke the target function */
    crs_tril_initialize();
    crs_tril1(&A, k, &L);
    crs_tril_terminate();

    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    int32_T _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_emxArray_to_mxArray((emxArray__common*)L.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_emxArray_to_mxArray((emxArray__common*)L.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_emxArray_to_mxArray((emxArray__common*)L.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&L.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&L.ncols, mxINT32_CLASS));

    /* Free temporary variables */
    free_emxArray( (emxArray__common*)A.val); mxFree( A.val);
    free_emxArray( (emxArray__common*)A.col_ind); mxFree( A.col_ind);
    free_emxArray( (emxArray__common*)A.row_ptr); mxFree( A.row_ptr);

    free_emxArray( (emxArray__common*)L.val); mxFree( L.val);
    free_emxArray( (emxArray__common*)L.col_ind); mxFree( L.col_ind);
    free_emxArray( (emxArray__common*)L.row_ptr); mxFree( L.row_ptr);

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if (nrhs == 1) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_tril:TooManyOutputArguments","Too many output arguments for entry-point 'crs_tril'.");

        /* Call the API function. */
        crs_tril_api(prhs, (const mxArray**)plhs);
    }

    else if (nrhs == 2) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_tril1:TooManyOutputArguments","Too many output arguments for entry-point 'crs_tril1'.");

        /* Call the API function. */
        crs_tril1_api(prhs, (const mxArray**)plhs);
    }

    else
        mexErrMsgIdAndTxt("crs_tril:WrongNumberOfInputs","Incorrect number of input variables for entry-point 'crs_tril'.");

}

