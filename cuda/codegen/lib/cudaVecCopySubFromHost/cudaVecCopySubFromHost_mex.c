/*
 * cudaVecCopySubFromHost_mex.c
 *
 * Auxiliary code for mexFunction of cudaVecCopySubFromHost
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {m2c_int, m2c_vec, m2c_int, m2c_int, CudaVec, m2c_int}  enableInline=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "cudaVecCopySubFromHost.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void cudaVecCopySubFromHost_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_real_T      vec;

    struct0_T            cuVec;
    mxArray              *_sub_mx1;

    int32_T              n;
    int32_T              istart;
    int32_T              inc;
    int32_T              incCuVec;
    int32_T              *errCode;
    boolean_T            *toplevel;

    /* Marshall in function inputs */
    if (mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument n has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument n should be a scalar.");
    n = *(int32_T*)mxGetData(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument vec has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(prhs[1], (emxArray__common *)&vec, "vec", 1);
    if (mxGetNumberOfElements(prhs[2]) && mxGetClassID(prhs[2]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument istart has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[2]) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument istart should be a scalar.");
    istart = *(int32_T*)mxGetData(prhs[2]);
    if (mxGetNumberOfElements(prhs[3]) && mxGetClassID(prhs[3]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument inc has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument inc should be a scalar.");
    inc = *(int32_T*)mxGetData(prhs[3]);

    if (!mxIsStruct(prhs[4]))
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument cuVec has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[4])!=3)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:InputStructWrongFields",
            "Input argument cuVec has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[4]) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument cuVec must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[4], 0, "data");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputStruct",
            "Input argument cuVec does not have the field data.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument cuVec.data has incorrect data type. uint64 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument cuVec.data should be a scalar.");
    cuVec.data = *(uint64_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[4], 0, "type");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputStruct",
            "Input argument cuVec does not have the field type.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument cuVec.type has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument cuVec.type should be a scalar.");
    cuVec.type = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[4], 0, "len");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputStruct",
            "Input argument cuVec does not have the field len.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument cuVec.len has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument cuVec.len should be a scalar.");
    cuVec.len = *(int32_T*)mxGetData(_sub_mx1);
    if (mxGetNumberOfElements(prhs[5]) && mxGetClassID(prhs[5]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongInputType",
            "Input argument incCuVec has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[5]) != 1)
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongSizeOfInputArg",
            "Argument incCuVec should be a scalar.");
    incCuVec = *(int32_T*)mxGetData(prhs[5]);

    /* Preallocate output variables */
    {mwSize l_size[] = {1, 1};
    *(void **)&errCode = prealloc_mxArray((mxArray**)&plhs[0], mxINT32_CLASS, 2, l_size); }
    {mwSize l_size[] = {1, 1};
    *(void **)&toplevel = prealloc_mxArray((mxArray**)&plhs[1], mxLOGICAL_CLASS, 2, l_size); }

    /* Invoke the target function */
    cudaVecCopySubFromHost_initialize();

    cudaVecCopySubFromHost(n, &vec, istart, inc, &cuVec, incCuVec, errCode, toplevel);

    cudaVecCopySubFromHost_terminate();

    /* Marshall out function outputs */
    /* Nothing to do for plhs[0] */
    /* Nothing to do for plhs[1] */

    /* Free temporary variables */
    free_emxArray((emxArray__common*)&vec);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Temporary copy for mex outputs. */
    mxArray *outputs[2];
    int i;
    int nOutputs = (nlhs < 1 ? 1 : nlhs);

    if (nrhs == 6) {
        if (nlhs > 2)
            mexErrMsgIdAndTxt("cudaVecCopySubFromHost:TooManyOutputArguments","Too many output arguments for entry-point cudaVecCopySubFromHost.");
        /* Call the API function. */
        cudaVecCopySubFromHost_api(prhs, (const mxArray**)outputs);
    }
    else
        mexErrMsgIdAndTxt("cudaVecCopySubFromHost:WrongNumberOfInputs","Incorrect number of input variables for entry-point cudaVecCopySubFromHost.");

    /* Copy over outputs to the caller. */
    for (i = 0; i < nOutputs; ++i) {
        plhs[i] = outputs[i];
    }
}