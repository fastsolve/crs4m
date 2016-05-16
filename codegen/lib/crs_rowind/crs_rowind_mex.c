/*
 * codegen/lib/crs_rowind/crs_rowind_mex.c
 *
 * Auxiliary code for mexFunction of crs_rowind
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0),[inf,1])}  enableInline=1 withOMP=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "crs_rowind.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void crs_rowind_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_int32_T     row_ptr;
    emxArray_int32_T     col_ind;
    emxArray_int32_T     row_ind;

    /* Marshall in function inputs */
    if (mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_rowind:WrongInputType",
            "Input argument row_ptr has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(prhs[0], (emxArray__common *)&row_ptr, "row_ptr", 1);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("crs_rowind:WrongInputType",
            "Input argument col_ind has incorrect data type. int32 is expected.");
    alias_mxArray_to_emxArray(prhs[1], (emxArray__common *)&col_ind, "col_ind", 1);

    /* Preallocate output variables */
    init_emxArray((emxArray__common*)&row_ind, 1);

    /* Invoke the target function */
    crs_rowind_initialize();

    crs_rowind(&row_ptr, &col_ind, &row_ind);

    crs_rowind_terminate();

    /* Marshall out function outputs */
    plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&row_ind, mxINT32_CLASS);

    /* Free temporary variables */
    free_emxArray((emxArray__common*)&row_ptr);
    free_emxArray((emxArray__common*)&col_ind);
    free_emxArray((emxArray__common*)&row_ind);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 2) {
        if (nlhs > 1)
            mexErrMsgIdAndTxt("crs_rowind:TooManyOutputArguments","Too many output arguments for entry-point crs_rowind.");
        /* Call the API function. */
        crs_rowind_api(prhs, (const mxArray**)plhs);
    }
    else
        mexErrMsgIdAndTxt("crs_rowind:WrongNumberOfInputs","Incorrect number of input variables for entry-point crs_rowind.");
}