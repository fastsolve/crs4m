/*
 * cuSynchronize_mex.c
 *
 * Auxiliary code for mexFunction of cuSynchronize
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {CuStreamHandle}  enableInline=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "cuSynchronize.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void cuSynchronize_api(const mxArray ** prhs, const mxArray **plhs) {

    struct0_T            stm;
    mxArray              *_sub_mx1;

    int32_T              *errCode;
    boolean_T            *toplevel;

    /* Marshall in function inputs */

    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputType",
            "Input argument stm has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[0])!=3)
        mexErrMsgIdAndTxt("cuSynchronize:InputStructWrongFields",
            "Input argument stm has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("cuSynchronize:WrongSizeOfInputArg",
            "Argument stm must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[0], 0, "data");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputStruct",
            "Input argument stm does not have the field data.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxUINT8_CLASS)
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputType",
            "Input argument stm.data has incorrect data type. uint8 is expected.");
    *(void**)&stm.data = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)stm.data, "stm.data", 1);
    _sub_mx1 = mxGetField(prhs[0], 0, "type");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputStruct",
            "Input argument stm does not have the field type.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxCHAR_CLASS)
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputType",
            "Input argument stm.type has incorrect data type. char is expected.");
    *(void**)&stm.type = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)stm.type, "stm.type", 2);
    _sub_mx1 = mxGetField(prhs[0], 0, "nitems");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputStruct",
            "Input argument stm does not have the field nitems.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuSynchronize:WrongInputType",
            "Input argument stm.nitems has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cuSynchronize:WrongSizeOfInputArg",
            "Argument stm.nitems should be a scalar.");
    stm.nitems = *(int32_T*)mxGetData(_sub_mx1);

    /* Preallocate output variables */
    {mwSize l_size[] = {1, 1};
    *(void **)&errCode = prealloc_mxArray((mxArray**)&plhs[0], mxINT32_CLASS, 2, l_size); }
    {mwSize l_size[] = {1, 1};
    *(void **)&toplevel = prealloc_mxArray((mxArray**)&plhs[1], mxLOGICAL_CLASS, 2, l_size); }

    /* Invoke the target function */
    cuSynchronize_initialize();

    cuSynchronize(&stm, errCode, toplevel);

    cuSynchronize_terminate();

    /* Marshall out function outputs */
    /* Nothing to do for plhs[0] */
    /* Nothing to do for plhs[1] */

    /* Free temporary variables */
    free_emxArray((emxArray__common*)stm.type); mxFree(stm.type);
    free_emxArray((emxArray__common*)stm.data); mxFree(stm.data);

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Temporary copy for mex outputs. */
    mxArray *outputs[2];
    int i;
    int nOutputs = (nlhs < 1 ? 1 : nlhs);

    if (nrhs == 1) {
        if (nlhs > 2)
            mexErrMsgIdAndTxt("cuSynchronize:TooManyOutputArguments","Too many output arguments for entry-point cuSynchronize.");
        /* Call the API function. */
        cuSynchronize_api(prhs, (const mxArray**)outputs);
    }
    else
        mexErrMsgIdAndTxt("cuSynchronize:WrongNumberOfInputs","Incorrect number of input variables for entry-point cuSynchronize.");

    /* Copy over outputs to the caller. */
    for (i = 0; i < nOutputs; ++i) {
        plhs[i] = outputs[i];
    }
}