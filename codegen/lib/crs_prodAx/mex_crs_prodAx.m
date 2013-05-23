% Build script for crs_prodAx
if ~isnewer( ['../../../crs_prodAx.' mexext], 'crs_prodAx_mex.c', 'crs_prodAx.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [mpicflag, mpildflag] = mpiflags;    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' mpicflag  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_prodAx_mex.c  ' output ' ../../../crs_prodAx '  ' ' mpildflag  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_prodAx/mex_crs_prodAx.m');
    eval(cmd);
else
    fprintf( 'crs_prodAx.%s is up to date.\n', mexext);
end
