% Build script for crs_prodAtx
if ~isnewer( ['../../../crs_prodAtx.' mexext], 'crs_prodAtx_mex.c', 'crs_prodAtx.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [mpicflag, mpildflag] = mpiflags;    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' mpicflag  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_prodAtx_mex.c  ' output ' ../../../crs_prodAtx '  ' ' mpildflag  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_prodAtx/mex_crs_prodAtx.m');
    eval(cmd);
else
    fprintf( 'crs_prodAtx.%s is up to date.\n', mexext);
end
