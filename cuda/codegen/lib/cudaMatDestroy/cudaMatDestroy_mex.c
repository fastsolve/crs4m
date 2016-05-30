/*
 * cudaMatDestroy_mex.c
 *
 * Auxiliary code for mexFunction of cudaMatDestroy
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {CudaMat}  enableInline=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "cudaMatDestroy.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void cudaMatDestroy_api(const mxArray ** prhs, const mxArray **plhs) {

    struct0_T            mat;
    mxArray              *_sub_mx1;

    int32_T              *errCode;
    boolean_T            *toplevel;

    /* Marshall in function inputs */

    if (!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputType",
            "Input argument mat has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[0])!=3)
        mexErrMsgIdAndTxt("cudaMatDestroy:InputStructWrongFields",
            "Input argument mat has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongSizeOfInputArg",
            "Argument mat must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[0], 0, "data");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputStruct",
            "Input argument mat does not have the field data.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputType",
            "Input argument mat.data has incorrect data type. uint64 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongSizeOfInputArg",
            "Argument mat.data should be a scalar.");
    mat.data = *(uint64_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[0], 0, "type");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputStruct",
            "Input argument mat does not have the field type.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputType",
            "Input argument mat.type has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongSizeOfInputArg",
            "Argument mat.type should be a scalar.");
    mat.type = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[0], 0, "dims");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputStruct",
            "Input argument mat does not have the field dims.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongInputType",
            "Input argument mat.dims has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 2)
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongSizeOfInputArg",
            "Argument mat.dims must contain 2 numbers.");
    copy_mxArray_to_array(_sub_mx1, mat.dims, 2);

    /* Preallocate output variables */
    {mwSize l_size[] = {1, 1};
    *(void **)&errCode = prealloc_mxArray((mxArray**)&plhs[0], mxINT32_CLASS, 2, l_size); }
    {mwSize l_size[] = {1, 1};
    *(void **)&toplevel = prealloc_mxArray((mxArray**)&plhs[1], mxLOGICAL_CLASS, 2, l_size); }

    /* Invoke the target function */
    cudaMatDestroy_initialize();

    cudaMatDestroy(&mat, errCode, toplevel);

    cudaMatDestroy_terminate();

    /* Marshall out function outputs */
    /* Nothing to do for plhs[0] */
    /* Nothing to do for plhs[1] */

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Temporary copy for mex outputs. */
    mxArray *outputs[2];
    int i;
    int nOutputs = (nlhs < 1 ? 1 : nlhs);

    if (nrhs == 1) {
        if (nlhs > 2)
            mexErrMsgIdAndTxt("cudaMatDestroy:TooManyOutputArguments","Too many output arguments for entry-point cudaMatDestroy.");
        /* Call the API function. */
        cudaMatDestroy_api(prhs, (const mxArray**)outputs);
    }
    else
        mexErrMsgIdAndTxt("cudaMatDestroy:WrongNumberOfInputs","Incorrect number of input variables for entry-point cudaMatDestroy.");

    /* Copy over outputs to the caller. */
    for (i = 0; i < nOutputs; ++i) {
        plhs[i] = outputs[i];
    }
}