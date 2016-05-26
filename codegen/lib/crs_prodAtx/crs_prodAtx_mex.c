/*
 * codegen/lib/crs_prodAtx/crs_prodAtx_mex.c
 *
 * Auxiliary code for mexFunction of crs_prodAtx
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf]), coder.typeof(int32(1), [1,1], [1,0])} crs_prodAtx_ser -args {crs_matrix, coder.typeof(0, [inf,inf])} crs_prodAtx_ser1 -args {crs_matrix, coder.typeof(0, [inf,inf]), coder.typeof(0, [inf,inf])}  enableInline=1 withOMP=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "crs_prodAtx.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void crs_prodAtx_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_real_T      x;
    emxArray_real_T      b;
    emxArray_int32_T     nthreads;

    struct0_T            A;
    mxArray              *_sub_mx1;

    /* Marshall in function inputs */

    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[0])!=5)
        mexErrMsgIdAndTxt("crs_prodAtx:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[0], 0, "row_ptr");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "col_ind");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "val");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputStruct",
            "Input argument A does not have the field val.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "nrows");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[0], 0, "ncols");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument x has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(prhs[1], (emxArray__common *)&x, "x", 2);
    plhs[0] = mxDuplicateArray(prhs[2]);
    if (mxGetNumberOfElements(plhs[0]) && mxGetClassID(plhs[0]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument b has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(plhs[0], (emxArray__common *)&b, "b", 2);
    if (mxGetNumberOfElements(prhs[3]) && mxGetClassID(prhs[3]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx:WrongInputType",
            "Input argument nthreads has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(prhs[3], (emxArray__common *)&nthreads, "nthreads", 1);

    /* Invoke the target function */
    crs_prodAtx_initialize();

    crs_prodAtx(&A, &x, &b, &nthreads);

    crs_prodAtx_terminate();

    /* Marshall out function outputs */
    if (b.canFreeData) plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&b, mxDOUBLE_CLASS);

    /* Free temporary variables */
    free_emxArray((emxArray__common*)A.val); mxFree(A.val);
    free_emxArray((emxArray__common*)A.col_ind); mxFree(A.col_ind);
    free_emxArray((emxArray__common*)A.row_ptr); mxFree(A.row_ptr);

    free_emxArray((emxArray__common*)&x);
    free_emxArray((emxArray__common*)&b);
    free_emxArray((emxArray__common*)&nthreads);
}
void crs_prodAtx_ser_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_real_T      x;
    emxArray_real_T      b;

    struct0_T            A;
    mxArray              *_sub_mx1;

    /* Marshall in function inputs */

    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[0])!=5)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[0], 0, "row_ptr");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "col_ind");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "val");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputStruct",
            "Input argument A does not have the field val.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "nrows");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[0], 0, "ncols");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser:WrongInputType",
            "Input argument x has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(prhs[1], (emxArray__common *)&x, "x", 2);

    /* Preallocate output variables */
    init_emxArray((emxArray__common*)&b, 2);

    /* Invoke the target function */
    crs_prodAtx_initialize();

    crs_prodAtx_ser(&A, &x, &b);

    crs_prodAtx_terminate();

    /* Marshall out function outputs */
    plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&b, mxDOUBLE_CLASS);

    /* Free temporary variables */
    free_emxArray((emxArray__common*)A.val); mxFree(A.val);
    free_emxArray((emxArray__common*)A.col_ind); mxFree(A.col_ind);
    free_emxArray((emxArray__common*)A.row_ptr); mxFree(A.row_ptr);

    free_emxArray((emxArray__common*)&x);
    free_emxArray((emxArray__common*)&b);
}
void crs_prodAtx_ser1_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_real_T      x;
    emxArray_real_T      b;

    struct0_T            A;
    mxArray              *_sub_mx1;

    /* Marshall in function inputs */

    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument A has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[0])!=5)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:InputStructWrongFields",
            "Input argument A has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongSizeOfInputArg",
            "Argument A must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[0], 0, "row_ptr");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputStruct",
            "Input argument A does not have the field row_ptr.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument A.row_ptr has incorrect data type. int32 is expected.");
    *(void**)&A.row_ptr = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.row_ptr, "A.row_ptr", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "col_ind");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputStruct",
            "Input argument A does not have the field col_ind.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument A.col_ind has incorrect data type. int32 is expected.");
    *(void**)&A.col_ind = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.col_ind, "A.col_ind", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "val");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputStruct",
            "Input argument A does not have the field val.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument A.val has incorrect data type. double is expected.");
    *(void**)&A.val = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)A.val, "A.val", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "nrows");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputStruct",
            "Input argument A does not have the field nrows.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument A.nrows has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongSizeOfInputArg",
            "Argument A.nrows should be a scalar.");
    A.nrows = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[0], 0, "ncols");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputStruct",
            "Input argument A does not have the field ncols.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument A.ncols has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongSizeOfInputArg",
            "Argument A.ncols should be a scalar.");
    A.ncols = *(int32_T*)mxGetData(_sub_mx1);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument x has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(prhs[1], (emxArray__common *)&x, "x", 2);
    plhs[0] = mxDuplicateArray(prhs[2]);
    if (mxGetNumberOfElements(plhs[0]) && mxGetClassID(plhs[0]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_prodAtx_ser1:WrongInputType",
            "Input argument b has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(plhs[0], (emxArray__common *)&b, "b", 2);

    /* Invoke the target function */
    crs_prodAtx_initialize();

    crs_prodAtx_ser1(&A, &x, &b);

    crs_prodAtx_terminate();

    /* Marshall out function outputs */
    if (b.canFreeData) plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&b, mxDOUBLE_CLASS);

    /* Free temporary variables */
    free_emxArray((emxArray__common*)A.val); mxFree(A.val);
    free_emxArray((emxArray__common*)A.col_ind); mxFree(A.col_ind);
    free_emxArray((emxArray__common*)A.row_ptr); mxFree(A.row_ptr);

    free_emxArray((emxArray__common*)&x);
    free_emxArray((emxArray__common*)&b);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 4) {
        if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_prodAtx:TooManyOutputArguments","Too many output arguments for entry-point crs_prodAtx.");
        /* Call the API function. */
        crs_prodAtx_api(prhs, (const mxArray**)plhs);
    }
    else if (nrhs == 2) {
        if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_prodAtx_ser:TooManyOutputArguments","Too many output arguments for entry-point crs_prodAtx_ser.");
        /* Call the API function. */
        crs_prodAtx_ser_api(prhs, (const mxArray**)plhs);
    }
    else if (nrhs == 3) {
        if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_prodAtx_ser1:TooManyOutputArguments","Too many output arguments for entry-point crs_prodAtx_ser1.");
        /* Call the API function. */
        crs_prodAtx_ser1_api(prhs, (const mxArray**)plhs);
    }
    else
        mexErrMsgIdAndTxt("crs_prodAtx:WrongNumberOfInputs","Incorrect number of input variables for entry-point crs_prodAtx.");
}