/*
 * crs_prodPtAP_mex.c
 *
 * Auxiliary code generation for function crs_prodPtAP
 *
 * C source code generated by lib2mex.
*
 */

#include "mex.h"

#define BUILD_MEX
/* Include the C file generated by codegen in lib mode */
#include "crs_prodPtAP.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"


void crs_prodPtAP_api(const mxArray ** prhs, const mxArray **plhs) {


    struct_T             A;
    struct_T             P;
    struct_T             B;
    mxArray              *_sub_mx1;


    /* Marshall in function inputs */

    if ( !mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if ( mxGetNumberOfFields( prhs[0])!=5)
        mexWarnMsgIdAndTxt("crs_prodPtAP:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if ( mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField( prhs[0], 0, "row_ptr");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "col_ind");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "val");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument A does not have the field val.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField( prhs[0], 0, "nrows");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField( prhs[0], 0, "ncols");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);

    if ( !mxIsStruct(prhs[1]))
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument P has incorrect data type. struct is expected.");
    if ( mxGetNumberOfFields( prhs[1])!=5)
        mexWarnMsgIdAndTxt("crs_prodPtAP:InputStructWrongFields",
            "Input argument P has incorrect number of fields.");
    if ( mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongSizeOfInputArg",
            "Argument P must contain 1 items.");

    _sub_mx1 = mxGetField( prhs[1], 0, "row_ptr");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument P does not have the field row_ptr.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument P.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&P.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)P.row_ptr, "P.row_ptr", 1);
    _sub_mx1 = mxGetField( prhs[1], 0, "col_ind");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument P does not have the field col_ind.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument P.col_ind has incorrect data type. int32 is expected.");
    *(void**)&P.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)P.col_ind, "P.col_ind", 1);
    _sub_mx1 = mxGetField( prhs[1], 0, "val");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument P does not have the field val.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument P.val has incorrect data type. double is expected.");
    *(void**)&P.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)P.val, "P.val", 1);
    _sub_mx1 = mxGetField( prhs[1], 0, "nrows");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument P does not have the field nrows.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument P.nrows has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongSizeOfInputArg",
            "Argument P.nrows should be a scalar.");
    P.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField( prhs[1], 0, "ncols");
    if ( _sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputStruct",
            "Input argument P does not have the field ncols.");
    if ( mxGetData(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongInputType",
            "Input argument P.ncols has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongSizeOfInputArg",
            "Argument P.ncols should be a scalar.");
    P.ncols = *(int32_T*)mxGetData(_sub_mx1);

    /* Preallocate output variables */
    *(void **)&B.row_ptr = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)B.row_ptr, 1);
    *(void **)&B.col_ind = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)B.col_ind, 1);
    *(void **)&B.val = mxCalloc(1, sizeof(emxArray__common));    init_emxArray( (emxArray__common*)B.val, 1);

    /* Invoke the target function */
    crs_prodPtAP_initialize();
    crs_prodPtAP(&A, &P, &B);
    crs_prodPtAP_terminate();

    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    int32_T _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_emxArray_to_mxArray((emxArray__common*)B.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_emxArray_to_mxArray((emxArray__common*)B.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_emxArray_to_mxArray((emxArray__common*)B.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&B.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&B.ncols, mxINT32_CLASS));

    /* Free temporary variables */
    free_emxArray( (emxArray__common*)A.val); mxFree( A.val);
    free_emxArray( (emxArray__common*)A.col_ind); mxFree( A.col_ind);
    free_emxArray( (emxArray__common*)A.row_ptr); mxFree( A.row_ptr);

    free_emxArray( (emxArray__common*)P.val); mxFree( P.val);
    free_emxArray( (emxArray__common*)P.col_ind); mxFree( P.col_ind);
    free_emxArray( (emxArray__common*)P.row_ptr); mxFree( P.row_ptr);

    free_emxArray( (emxArray__common*)B.val); mxFree( B.val);
    free_emxArray( (emxArray__common*)B.col_ind); mxFree( B.col_ind);
    free_emxArray( (emxArray__common*)B.row_ptr); mxFree( B.row_ptr);

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if (nrhs == 2) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_prodPtAP:TooManyOutputArguments","Too many output arguments for entry-point 'crs_prodPtAP'.");

        /* Call the API function. */
        crs_prodPtAP_api(prhs, (const mxArray**)plhs);
    }

    else
        mexErrMsgIdAndTxt("crs_prodPtAP:WrongNumberOfInputs","Incorrect number of input variables for entry-point 'crs_prodPtAP'.");

}

