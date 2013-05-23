% Build script for crs_create
if ~isnewer( ['../../../crs_create.' mexext], 'crs_create_mex.c', 'crs_create.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_create_mex.c  ' output ' ../../../crs_create '  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_create/mex_crs_create.m');
    eval(cmd);
else
    fprintf( 'crs_create.%s is up to date.\n', mexext);
end
