% Build script for crs_diag
if ~isnewer( ['../../../crs_diag.' mexext], 'crs_diag_mex.c', 'crs_diag.c')
    if ~exist('dbopts.m', 'file'); dbopts = ''; end
    dir = which('m2c.m'); dir=dir(1:end-6);


    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' dbopts ' -I"' dir '" -I. '  'crs_diag_mex.c  ' output ' ../../../crs_diag ' ];
    disp( 'run codegen/lib/crs_diag/mex_crs_diag.m');
    eval(cmd);
else
    fprintf( 'crs_diag.%s is up to date.\n', mexext);
end
