% Build script for crs_prodAx
if ~isnewer( ['../../../crs_prodAx.' mexext], 'crs_prodAx_mex.c', 'crs_prodAx.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('m2c.m'); dir=dir(1:end-6);
    [mpicflag, mpildflag] = mpiflags;
    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' mpicflag  ' ' ompcflag ' ' dbopts ' -I"' dir '" -I. '  'crs_prodAx_mex.c  ' output ' ../../../crs_prodAx '  ' ' mpildflag  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_prodAx/mex_crs_prodAx.m');
    eval(cmd);
else
    fprintf( 'crs_prodAx.%s is up to date.\n', mexext);
end
