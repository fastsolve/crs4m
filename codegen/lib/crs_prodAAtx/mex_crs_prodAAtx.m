% Build script for crs_prodAAtx
if ~isnewer( ['../../../crs_prodAAtx.' mexext], 'crs_prodAAtx_mex.c', 'crs_prodAAtx.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [mpicflag, mpildflag] = mpiflags;    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' mpicflag  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_prodAAtx_mex.c  ' output ' ../../../crs_prodAAtx '  ' ' mpildflag  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_prodAAtx/mex_crs_prodAAtx.m');
    eval(cmd);
else
    fprintf( 'crs_prodAAtx.%s is up to date.\n', mexext);
end
