% Build script for crs_prodAtx
if ~isnewer( ['../../../crs_prodAtx.' mexext], 'crs_prodAtx_mex.c', 'crs_prodAtx.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('lib2mex'); dir=dir(1:end-10);
    [mpicflag, mpildflag] = mpiflags;
    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' mpicflag  ' ' ompcflag ' ' dbopts ' -I"' dir '" -I. '  'crs_prodAtx_mex.c  ' output ' ../../../crs_prodAtx '  ' ' mpildflag  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_prodAtx/mex_crs_prodAtx.m');
    eval(cmd);
else
    fprintf( 'crs_prodAtx.%s is up to date.\n', mexext);
end
