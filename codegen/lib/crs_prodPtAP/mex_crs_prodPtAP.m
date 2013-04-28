% Build script for crs_prodPtAP
if ~isnewer( ['../../../crs_prodPtAP.' mexext], 'crs_prodPtAP_mex.c', 'crs_prodPtAP.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);

    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' ompcflag ' ' dbopts ' -I"' dir 'include" -I. '  'crs_prodPtAP_mex.c  ' output ' ../../../crs_prodPtAP '  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_prodPtAP/mex_crs_prodPtAP.m');
    eval(cmd);
else
    fprintf( 'crs_prodPtAP.%s is up to date.\n', mexext);
end
