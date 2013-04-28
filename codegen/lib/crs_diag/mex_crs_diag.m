% Build script for crs_diag
if ~isnewer( ['../../../crs_diag.' mexext], 'crs_diag_mex.c', 'crs_diag.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);

    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' ompcflag ' ' dbopts ' -I"' dir 'include" -I. '  'crs_diag_mex.c  ' output ' ../../../crs_diag '  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_diag/mex_crs_diag.m');
    eval(cmd);
else
    fprintf( 'crs_diag.%s is up to date.\n', mexext);
end
