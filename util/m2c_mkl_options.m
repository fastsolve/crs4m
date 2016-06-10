function opts = m2c_mkl_options(version)
% Try to locate MKL automatically using MKLROOT.
% The output is either {}, {'-mkl'}, or {'-mkl', {MKLPATH}}, which
% can be passed as arguments to m2c using m2c(opts{:}, ...)

if nargin==0
    version = '';
end

if getenv('MKLROOT')
    opts = {['-mkl' version]};
else
    if exist('/opt/intel/mkl/lib/libmkl_core.a', 'file')
        opts = {['-mkl' version], '/opt/intel/mkl'};
    else
        opts = {};
    end
end
