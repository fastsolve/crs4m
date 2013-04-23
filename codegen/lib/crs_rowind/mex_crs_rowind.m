% Build script for crs_rowind
if ~isnewer( ['../../../crs_rowind.' mexext], 'crs_rowind_mex.c', 'crs_rowind.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('lib2mex'); dir=dir(1:end-10);

    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' ompcflag ' ' dbopts ' -I"' dir '" -I. '  'crs_rowind_mex.c  ' output ' ../../../crs_rowind '  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_rowind/mex_crs_rowind.m');
    eval(cmd);
else
    fprintf( 'crs_rowind.%s is up to date.\n', mexext);
end
