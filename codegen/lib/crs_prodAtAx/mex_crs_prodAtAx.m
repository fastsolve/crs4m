% Build script for crs_prodAtAx
if ~isnewer( ['../../../crs_prodAtAx.' mexext], 'crs_prodAtAx_mex.c', 'crs_prodAtAx.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [mpicflag, mpildflag] = mpiflags;
    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' mpicflag  ' ' ompcflag ' ' dbopts ' -I"' dir 'include" -I. '  'crs_prodAtAx_mex.c  ' output ' ../../../crs_prodAtAx '  ' ' mpildflag  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_prodAtAx/mex_crs_prodAtAx.m');
    eval(cmd);
else
    fprintf( 'crs_prodAtAx.%s is up to date.\n', mexext);
end
