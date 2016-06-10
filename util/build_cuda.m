function build_cuda(varargin)
%build_msp Build script for Cuda submodule

mspackroot = fileparts(which('startup_mspack'));
curpath = pwd;
cd(mspackroot);

opts = m2c_cuda_options;
if isempty(opts)
    warning('build_cuda:NoCUDA_PATH', 'CUDA_PATH was not set and nvcc is not in path. Compiling without CUDA.\n');
end

try
    % First, compile test scripts
    lines = [grep_pattern('cuda/*.m', '\n%#codegen -args') ...
        grep_pattern('cuda/sparse/*.m', '\n%#codegen -args')];
    files = regexp(lines, '(\w+.m):', 'tokens');
    
    for j=1:length(files)
        file = files{j}{1};
        m2c(opts{:}, '-mex', '-noinf', '-O', ['-I' mspackroot '/include'], ...
            varargin{:}, file);
    end
catch ME
    cd(curpath);
    rethrow(ME)
end

cd(curpath);
