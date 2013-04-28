% Build script for crs_prodAB
if ~isnewer( ['../../../crs_prodAB.' mexext], 'crs_prodAB_mex.c', 'crs_prodAB.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);

    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' ompcflag ' ' dbopts ' -I"' dir 'include" -I. '  'crs_prodAB_mex.c  ' output ' ../../../crs_prodAB '  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_prodAB/mex_crs_prodAB.m');
    eval(cmd);
else
    fprintf( 'crs_prodAB.%s is up to date.\n', mexext);
end
