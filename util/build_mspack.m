function build_mspack(varargin)
%build_msp Build script for MSP

mspackroot = fileparts(which('startup_mspack'));
curpath = pwd;
cd(mspackroot);

try
    cuda_opts = m2c_cuda_options;
    build_cuda(varargin{:});
    
    mkl_opts = m2c_mkl_options;
    if isempty(mkl_opts)
        warning('build_mspack:NoMKL', 'MKLROOT was not set. Compiling without MKL.\n');
    end
    
    %grep_pattern('Mat/*.m', '\n%#codegen -args') ...
    % First, compile test scripts
    lines = [grep_pattern('Vec/*.m', '\n%#codegen -args') ...
        grep_pattern('Mat/BCRS/*.m', '\n%#codegen -args') ...
        grep_pattern('KSP/*.m', '\n%#codegen -args') ...
        grep_pattern('PC/*.m', '\n%#codegen -args')];
    files = regexp(lines, '(\w+.m):', 'tokens');
    
    for j=1:length(files)
        file = files{j}{1};
        m2c('-omp', '-mex', mkl_opts{:}, cuda_opts{:}, '-noinf', '-O', ...
            ['-I' mspackroot '/include'], varargin{:}, file);
    end
catch ME
    cd(curpath);
    rethrow(ME)
end

cd(curpath);
