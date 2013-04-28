% Build script for crs_tril
if ~isnewer( ['../../../crs_tril.' mexext], 'crs_tril_mex.c', 'crs_tril.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);

    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' ompcflag ' ' dbopts ' -I"' dir 'include" -I. '  'crs_tril_mex.c  ' output ' ../../../crs_tril '  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_tril/mex_crs_tril.m');
    eval(cmd);
else
    fprintf( 'crs_tril.%s is up to date.\n', mexext);
end
