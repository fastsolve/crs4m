% Build script for crs_transp
if ~isnewer( ['../../../crs_transp.' mexext], 'crs_transp_mex.c', 'crs_transp.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_transp_mex.c  ' output ' ../../../crs_transp '  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_transp/mex_crs_transp.m');
    eval(cmd);
else
    fprintf( 'crs_transp.%s is up to date.\n', mexext);
end
