/*
 * cuMatCopySubFromGPU_mex.c
 *
 * Auxiliary code for mexFunction of cuMatCopySubFromGPU
 *
 * C source code generated by m2c.
 * %#m2c options: codegenArgs=-args {m2c_int, m2c_int, CuMat, m2c_mat} cuMatCopySubFromGPU_async -args {m2c_int, m2c_int, CuMat, m2c_mat, CuStreamHandle}  enableInline=1
 *
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#else
#include "matrix.h"
#endif
/* Include the C file generated by codegen in lib mode */
#include "cuMatCopySubFromGPU.h"
#include "m2c.c"
/* Include declaration of some helper functions. */
#include "lib2mex_helper.c"

void cuMatCopySubFromGPU_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_real_T      mat;

    struct0_T            cuMat;
    mxArray              *_sub_mx1;

    int32_T              nrows;
    int32_T              ncols;
    int32_T              *errCode;
    boolean_T            *toplevel;

    /* Marshall in function inputs */
    if (mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument nrows has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongSizeOfInputArg",
            "Argument nrows should be a scalar.");
    nrows = *(int32_T*)mxGetData(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument ncols has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongSizeOfInputArg",
            "Argument ncols should be a scalar.");
    ncols = *(int32_T*)mxGetData(prhs[1]);

    if (!mxIsStruct(prhs[2]))
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument cuMat has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[2])!=3)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:InputStructWrongFields",
            "Input argument cuMat has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[2]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongSizeOfInputArg",
            "Argument cuMat must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[2], 0, "data");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputStruct",
            "Input argument cuMat does not have the field data.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument cuMat.data has incorrect data type. uint64 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongSizeOfInputArg",
            "Argument cuMat.data should be a scalar.");
    cuMat.data = *(uint64_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[2], 0, "type");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputStruct",
            "Input argument cuMat does not have the field type.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument cuMat.type has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongSizeOfInputArg",
            "Argument cuMat.type should be a scalar.");
    cuMat.type = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[2], 0, "dims");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputStruct",
            "Input argument cuMat does not have the field dims.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument cuMat.dims has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 2)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongSizeOfInputArg",
            "Argument cuMat.dims must contain 2 numbers.");
    copy_mxArray_to_array(_sub_mx1, cuMat.dims, 2);
    plhs[0] = mxDuplicateArray(prhs[3]);
    if (mxGetNumberOfElements(plhs[0]) && mxGetClassID(plhs[0]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongInputType",
            "Input argument mat has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(plhs[0], (emxArray__common *)&mat, "mat", 2);

    /* Preallocate output variables */
    {mwSize l_size[] = {1, 1};
    *(void **)&errCode = prealloc_mxArray((mxArray**)&plhs[1], mxINT32_CLASS, 2, l_size); }
    {mwSize l_size[] = {1, 1};
    *(void **)&toplevel = prealloc_mxArray((mxArray**)&plhs[2], mxLOGICAL_CLASS, 2, l_size); }

    /* Invoke the target function */
    cuMatCopySubFromGPU_initialize();

    cuMatCopySubFromGPU(nrows, ncols, &cuMat, &mat, errCode, toplevel);

    cuMatCopySubFromGPU_terminate();

    /* Marshall out function outputs */
    if (mat.canFreeData) plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&mat, mxDOUBLE_CLASS);
    /* Nothing to do for plhs[1] */
    /* Nothing to do for plhs[2] */

    /* Free temporary variables */
    free_emxArray((emxArray__common*)&mat);
}
void cuMatCopySubFromGPU_async_api(const mxArray ** prhs, const mxArray **plhs) {

    emxArray_real_T      mat;

    struct0_T            cuMat;
    struct1_T            strm;
    mxArray              *_sub_mx1;

    int32_T              nrows;
    int32_T              ncols;
    int32_T              *errCode;
    boolean_T            *toplevel;

    /* Marshall in function inputs */
    if (mxGetNumberOfElements(prhs[0]) && mxGetClassID(prhs[0]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument nrows has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[0]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument nrows should be a scalar.");
    nrows = *(int32_T*)mxGetData(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) && mxGetClassID(prhs[1]) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument ncols has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument ncols should be a scalar.");
    ncols = *(int32_T*)mxGetData(prhs[1]);

    if (!mxIsStruct(prhs[2]))
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument cuMat has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[2])!=3)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:InputStructWrongFields",
            "Input argument cuMat has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[2]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument cuMat must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[2], 0, "data");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputStruct",
            "Input argument cuMat does not have the field data.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxUINT64_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument cuMat.data has incorrect data type. uint64 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument cuMat.data should be a scalar.");
    cuMat.data = *(uint64_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[2], 0, "type");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputStruct",
            "Input argument cuMat does not have the field type.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument cuMat.type has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument cuMat.type should be a scalar.");
    cuMat.type = *(int32_T*)mxGetData(_sub_mx1);
    _sub_mx1 = mxGetField(prhs[2], 0, "dims");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputStruct",
            "Input argument cuMat does not have the field dims.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument cuMat.dims has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 2)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument cuMat.dims must contain 2 numbers.");
    copy_mxArray_to_array(_sub_mx1, cuMat.dims, 2);
    plhs[0] = mxDuplicateArray(prhs[3]);
    if (mxGetNumberOfElements(plhs[0]) && mxGetClassID(plhs[0]) != mxDOUBLE_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument mat has incorrect data type. double is expected.");
    alias_mxArray_to_emxArray(plhs[0], (emxArray__common *)&mat, "mat", 2);

    if (!mxIsStruct(prhs[4]))
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument strm has incorrect data type. struct is expected.");
    if (mxGetNumberOfFields(prhs[4])!=3)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:InputStructWrongFields",
            "Input argument strm has incorrect number of fields.");
    if (mxGetNumberOfElements(prhs[4]) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument strm must contain 1 items.");

    _sub_mx1 = mxGetField(prhs[4], 0, "data");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputStruct",
            "Input argument strm does not have the field data.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxUINT8_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument strm.data has incorrect data type. uint8 is expected.");
    *(void**)&strm.data = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)strm.data, "strm.data", 1);
    _sub_mx1 = mxGetField(prhs[4], 0, "type");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputStruct",
            "Input argument strm does not have the field type.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxCHAR_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument strm.type has incorrect data type. char is expected.");
    *(void**)&strm.type = mxCalloc(1, sizeof(emxArray__common));
    alias_mxArray_to_emxArray(_sub_mx1, (emxArray__common*)strm.type, "strm.type", 2);
    _sub_mx1 = mxGetField(prhs[4], 0, "nitems");
    if (_sub_mx1==NULL)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputStruct",
            "Input argument strm does not have the field nitems.");
    if (mxGetNumberOfElements(_sub_mx1) && mxGetClassID(_sub_mx1) != mxINT32_CLASS)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongInputType",
            "Input argument strm.nitems has incorrect data type. int32 is expected.");
    if (mxGetNumberOfElements(_sub_mx1) != 1)
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:WrongSizeOfInputArg",
            "Argument strm.nitems should be a scalar.");
    strm.nitems = *(int32_T*)mxGetData(_sub_mx1);

    /* Preallocate output variables */
    {mwSize l_size[] = {1, 1};
    *(void **)&errCode = prealloc_mxArray((mxArray**)&plhs[1], mxINT32_CLASS, 2, l_size); }
    {mwSize l_size[] = {1, 1};
    *(void **)&toplevel = prealloc_mxArray((mxArray**)&plhs[2], mxLOGICAL_CLASS, 2, l_size); }

    /* Invoke the target function */
    cuMatCopySubFromGPU_initialize();

    cuMatCopySubFromGPU_async(nrows, ncols, &cuMat, &mat, &strm, errCode, toplevel);

    cuMatCopySubFromGPU_terminate();

    /* Marshall out function outputs */
    if (mat.canFreeData) plhs[0] = move_emxArray_to_mxArray((emxArray__common*)&mat, mxDOUBLE_CLASS);
    /* Nothing to do for plhs[1] */
    /* Nothing to do for plhs[2] */

    /* Free temporary variables */
    free_emxArray((emxArray__common*)&mat);
    free_emxArray((emxArray__common*)strm.type); mxFree(strm.type);
    free_emxArray((emxArray__common*)strm.data); mxFree(strm.data);

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Temporary copy for mex outputs. */
    mxArray *outputs[3];
    int i;
    int nOutputs = (nlhs < 1 ? 1 : nlhs);

    if (nrhs == 4) {
        if (nlhs > 3)
            mexErrMsgIdAndTxt("cuMatCopySubFromGPU:TooManyOutputArguments","Too many output arguments for entry-point cuMatCopySubFromGPU.");
        /* Call the API function. */
        cuMatCopySubFromGPU_api(prhs, (const mxArray**)outputs);
    }
    else if (nrhs == 5) {
        if (nlhs > 3)
            mexErrMsgIdAndTxt("cuMatCopySubFromGPU_async:TooManyOutputArguments","Too many output arguments for entry-point cuMatCopySubFromGPU_async.");
        /* Call the API function. */
        cuMatCopySubFromGPU_async_api(prhs, (const mxArray**)outputs);
    }
    else
        mexErrMsgIdAndTxt("cuMatCopySubFromGPU:WrongNumberOfInputs","Incorrect number of input variables for entry-point cuMatCopySubFromGPU.");

    /* Copy over outputs to the caller. */
    for (i = 0; i < nOutputs; ++i) {
        plhs[i] = outputs[i];
    }
}