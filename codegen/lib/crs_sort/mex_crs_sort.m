% Build script for crs_sort
if ~isnewer( ['../../../crs_sort.' mexext], 'crs_sort_mex.c', 'crs_sort.c')
    if ~exist('dbopts.m', 'file'); dbopts = ' -O'; end
    dir = which('opaque_obj.m'); dir=dir(1:end-12);

    [ompcflag, ompldflag] = ompflags;
    if exist('octave_config_info', 'builtin'); output = '-o'; else output = '-largeArrayDims -output'; end
    cmd = ['mex'  ' ' ompcflag ' ' dbopts ' -I"' dir 'include" -I. '  'crs_sort_mex.c  ' output ' ../../../crs_sort '  ' ' ompldflag ];
    disp( 'run codegen/lib/crs_sort/mex_crs_sort.m');
    eval(cmd);
else
    fprintf( 'crs_sort.%s is up to date.\n', mexext);
end
