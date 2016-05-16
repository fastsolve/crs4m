/*
 * codegen/lib/crs_sort/crs_sort_mex.c
 *
 * Auxiliary code for mexFunction of crs_sort
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]), coder.typeof(0, [inf,1])} crs_sort0 -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf, 1])}  enableInline=1 withOMP=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "crs_sort.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void crs_sort_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_int32_T     row_ptr;
    emxArray_int32_T     col_ind;
    emxArray_real_T      val;

    /* Marshall in function inputs */
    if (mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_sort:WrongInputType",
            "Input argument row_ptr has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(prhs[0], (emxArray__common *)&row_ptr, "row_ptr", 1);
    plhs[0] = mxDuplicateArray(prhs[1]);
    if (mxGetNumberOfElements(plhs[0]) && mxGetClassID(plhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_sort:WrongInputType",
            "Input argument col_ind has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(plhs[0], (emxArray__common *)&col_ind, "col_ind", 1);
    plhs[1] = mxDuplicateArray(prhs[2]);
    if (mxGetNumberOfElements(plhs[1]) && mxGetClassID(plhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("crs_sort:WrongInputType",
            "Input argument val has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(plhs[1], (emxArray__common *)&val, "val", 1);

    /* Invoke the target function */
    crs_sort_initialize();

    crs_sort(&row_ptr, &col_ind, &val);

    crs_sort_terminate();

    /* Marshall out function outputs */
    if (col_ind.canFreeData) plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&col_ind, mxINT32_CLASS);
    if (val.canFreeData) plhs[1] = move_emxArray_to_mxArray((emxArray__common*)&val, mxDOUBLE_CLASS);

    /* Free temporary variables */
    free_emxArray((emxArray__common*)&row_ptr);
    free_emxArray((emxArray__common*)&col_ind);
    free_emxArray((emxArray__common*)&val);
}
void crs_sort0_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_int32_T     row_ptr;
    emxArray_int32_T     col_ind;

    /* Marshall in function inputs */
    if (mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_sort0:WrongInputType",
            "Input argument row_ptr has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(prhs[0], (emxArray__common *)&row_ptr, "row_ptr", 1);
    plhs[0] = mxDuplicateArray(prhs[1]);
    if (mxGetNumberOfElements(plhs[0]) && mxGetClassID(plhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_sort0:WrongInputType",
            "Input argument col_ind has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(plhs[0], (emxArray__common *)&col_ind, "col_ind", 1);

    /* Invoke the target function */
    crs_sort_initialize();

    crs_sort0(&row_ptr, &col_ind);

    crs_sort_terminate();

    /* Marshall out function outputs */
    if (col_ind.canFreeData) plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&col_ind, mxINT32_CLASS);

    /* Free temporary variables */
    free_emxArray((emxArray__common*)&row_ptr);
    free_emxArray((emxArray__common*)&col_ind);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Temporary copy for mex outputs. */
    mxArray *outputs[2];
    int i;
    int nOutputs = (nlhs < 1 ? 1 : nlhs);

    if (nrhs == 3) {
        if (nlhs > 2)
            mexErrMsgIdAndTxt("crs_sort:TooManyOutputArguments","Too many output arguments for entry-point crs_sort.");
        /* Call the API function. */
        crs_sort_api(prhs, (const mxArray**)outputs);
    }
    else if (nrhs == 2) {
        if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_sort0:TooManyOutputArguments","Too many output arguments for entry-point crs_sort0.");
        /* Call the API function. */
        crs_sort0_api(prhs, (const mxArray**)outputs);
    }
    else
        mexErrMsgIdAndTxt("crs_sort:WrongNumberOfInputs","Incorrect number of input variables for entry-point crs_sort.");

    /* Copy over outputs to the caller. */
    for (i = 0; i < nOutputs; ++i) {
        plhs[i] = outputs[i];
    }
}