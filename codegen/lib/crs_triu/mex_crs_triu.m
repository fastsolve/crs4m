% Build script for crs_triu
if ~isnewer( ['../../../crs_triu.' mexext], 'crs_triu_mex.c', 'crs_triu.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_triu_mex.c  ' output ' ../../../crs_triu '  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_triu/mex_crs_triu.m');
    eval(cmd);
else
    fprintf( 'crs_triu.%s is up to date.\n', mexext);
end
