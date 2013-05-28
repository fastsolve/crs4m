/*
 * crs_create_mex.c
 *
 * Auxiliary code generation for function crs_create
 *
 * C source code generated by palc.
*
 */

#include "mex.h"

#define BUILD_MEX
/* Include the C file generated by codegen in lib mode */
#include "crs_create.h"
#include "palc.c"
#include "crs_create.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void crs_create_api(const mxArray ** prhs, const mxArray **plhs) {

    plcArray_int32_T     rows;
    plcArray_int32_T     cols;
    plcArray_real_T      vs;

    struct_T             A;
    mxArray              *_sub_mx1;

    /* Marshall in function inputs */
    if ( mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create:WrongInputType",
            "Input argument rows has incorrect data type. int32 is expected.");
    alias_mxArray_to_plcArray(prhs[0], (plcArray__common *)&rows, "rows", 1);
    if ( mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create:WrongInputType",
            "Input argument cols has incorrect data type. int32 is expected.");
    alias_mxArray_to_plcArray(prhs[1], (plcArray__common *)&cols, "cols", 1);
    if ( mxGetNumberOfElements(prhs[2]) && mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_create:WrongInputType",
            "Input argument vs has incorrect data type. double is expected.");
    alias_mxArray_to_plcArray(prhs[2], (plcArray__common *)&vs, "vs", 1);

    /* Preallocate output variables */
    *(void **)&A.row_ptr = mxCalloc(1, sizeof(plcArray__common));    init_plcArray( (plcArray__common*)A.row_ptr, 1);
    *(void **)&A.col_ind = mxCalloc(1, sizeof(plcArray__common));    init_plcArray( (plcArray__common*)A.col_ind, 1);
    *(void **)&A.val = mxCalloc(1, sizeof(plcArray__common));    init_plcArray( (plcArray__common*)A.val, 1);

    /* Invoke the target function */
    crs_create_initialize();
    crs_create(&rows, &cols, &vs, &A);
    crs_create_terminate();
    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    plcSize _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_plcArray_to_mxArray((plcArray__common*)A.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_plcArray_to_mxArray((plcArray__common*)A.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_plcArray_to_mxArray((plcArray__common*)A.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&A.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&A.ncols, mxINT32_CLASS));    /* Free temporary variables */
    free_plcArray( (plcArray__common*)&rows);
    free_plcArray( (plcArray__common*)&cols);
    free_plcArray( (plcArray__common*)&vs);
    free_plcArray( (plcArray__common*)A.val); mxFree( A.val);
    free_plcArray( (plcArray__common*)A.col_ind); mxFree( A.col_ind);
    free_plcArray( (plcArray__common*)A.row_ptr); mxFree( A.row_ptr);
}

void crs_create0_api(const mxArray ** prhs, const mxArray **plhs) {

    int32_T              ni;
    int32_T              nj;
    b_struct_T           A;

    /* Marshall in function inputs */
    if ( mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create0:WrongInputType",
            "Input argument ni has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_create0:WrongSizeOfInputArg",
            "Argument ni should be a scalar.");
    ni = *(int32_T*)mxGetData(prhs[0]);
    if ( mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create0:WrongInputType",
            "Input argument nj has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("crs_create0:WrongSizeOfInputArg",
            "Argument nj should be a scalar.");
    nj = *(int32_T*)mxGetData(prhs[1]);

    /* Preallocate output variables */

    /* Invoke the target function */
    crs_create_initialize();
    A = crs_create0(ni, nj);
    crs_create_terminate();
    /* Marshall out function outputs */
    {const char *_fields[] = { "nrows", "ncols",  ""};
    plcSize _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 2, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, copy_scalar_to_mxArray(&A.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, copy_scalar_to_mxArray(&A.ncols, mxINT32_CLASS));}

void crs_create1_api(const mxArray ** prhs, const mxArray **plhs) {

    plcArray_int32_T     is;
    plcArray_int32_T     js;
    plcArray_real_T      vs;

    struct_T             A;
    mxArray              *_sub_mx1;

    int32_T              ni;
    int32_T              nj;

    /* Marshall in function inputs */
    if ( mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create1:WrongInputType",
            "Input argument is has incorrect data type. int32 is expected.");
    alias_mxArray_to_plcArray(prhs[0], (plcArray__common *)&is, "is", 1);
    if ( mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create1:WrongInputType",
            "Input argument js has incorrect data type. int32 is expected.");
    alias_mxArray_to_plcArray(prhs[1], (plcArray__common *)&js, "js", 1);
    if ( mxGetNumberOfElements(prhs[2]) && mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_create1:WrongInputType",
            "Input argument vs has incorrect data type. double is expected.");
    alias_mxArray_to_plcArray(prhs[2], (plcArray__common *)&vs, "vs", 1);
    if ( mxGetNumberOfElements(prhs[3]) && mxGetClassID(prhs[3]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create1:WrongInputType",
            "Input argument ni has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgIdAndTxt("crs_create1:WrongSizeOfInputArg",
            "Argument ni should be a scalar.");
    ni = *(int32_T*)mxGetData(prhs[3]);
    if ( mxGetNumberOfElements(prhs[4]) && mxGetClassID(prhs[4]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_create1:WrongInputType",
            "Input argument nj has incorrect data type. int32 is expected.");
    if ( mxGetNumberOfElements(prhs[4]) != 1)
        mexErrMsgIdAndTxt("crs_create1:WrongSizeOfInputArg",
            "Argument nj should be a scalar.");
    nj = *(int32_T*)mxGetData(prhs[4]);

    /* Preallocate output variables */
    *(void **)&A.row_ptr = mxCalloc(1, sizeof(plcArray__common));    init_plcArray( (plcArray__common*)A.row_ptr, 1);
    *(void **)&A.col_ind = mxCalloc(1, sizeof(plcArray__common));    init_plcArray( (plcArray__common*)A.col_ind, 1);
    *(void **)&A.val = mxCalloc(1, sizeof(plcArray__common));    init_plcArray( (plcArray__common*)A.val, 1);

    /* Invoke the target function */
    crs_create_initialize();
    crs_create1(&is, &js, &vs, ni, nj, &A);
    crs_create_terminate();
    /* Marshall out function outputs */
    {const char *_fields[] = { "row_ptr", "col_ind", "val", "nrows", "ncols",  ""};
    plcSize _one=1;
    plhs[0] = create_struct_mxArray( 1, &_one, 5, _fields);}
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 0, move_plcArray_to_mxArray((plcArray__common*)A.row_ptr, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 1, move_plcArray_to_mxArray((plcArray__common*)A.col_ind, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 2, move_plcArray_to_mxArray((plcArray__common*)A.val, mxDOUBLE_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 3, copy_scalar_to_mxArray(&A.nrows, mxINT32_CLASS));
    mxSetFieldByNumber( (mxArray*)(plhs[0]), 0, 4, copy_scalar_to_mxArray(&A.ncols, mxINT32_CLASS));    /* Free temporary variables */
    free_plcArray( (plcArray__common*)&is);
    free_plcArray( (plcArray__common*)&js);
    free_plcArray( (plcArray__common*)&vs);
    free_plcArray( (plcArray__common*)A.val); mxFree( A.val);
    free_plcArray( (plcArray__common*)A.col_ind); mxFree( A.col_ind);
    free_plcArray( (plcArray__common*)A.row_ptr); mxFree( A.row_ptr);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 3) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_create:TooManyOutputArguments","Too many output arguments for entry-point 'crs_create'.");
        /* Call the API function. */
        crs_create_api(prhs, (const mxArray**)plhs);
    }
    else if (nrhs == 2) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_create0:TooManyOutputArguments","Too many output arguments for entry-point 'crs_create0'.");
        /* Call the API function. */
        crs_create0_api(prhs, (const mxArray**)plhs);
    }
    else if (nrhs == 5) {
         if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_create1:TooManyOutputArguments","Too many output arguments for entry-point 'crs_create1'.");
        /* Call the API function. */
        crs_create1_api(prhs, (const mxArray**)plhs);
    }
    else
        mexErrMsgIdAndTxt("crs_create:WrongNumberOfInputs","Incorrect number of input variables for entry-point 'crs_create'.");
}