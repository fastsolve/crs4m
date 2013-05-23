% Build script for crs_rowind
if ~isnewer( ['../../../crs_rowind.' mexext], 'crs_rowind_mex.c', 'crs_rowind.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);
    [accCFLAG, accLDFLAG] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' accCFLAG ' ' dbopts ' -I"' dir 'include" '  'crs_rowind_mex.c  ' output ' ../../../crs_rowind '  ' ' accLDFLAG ];
    disp( 'run codegen/lib/crs_rowind/mex_crs_rowind.m');
    eval(cmd);
else
    fprintf( 'crs_rowind.%s is up to date.\n', mexext);
end
