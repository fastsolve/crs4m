function build_cuda(varargin)
%build_msp Build script for Cuda submodule

mspackroot = fileparts(which('startup_mspack'));
curpath = pwd;
cd(mspackroot);

try
    % First, compile test scripts
    lines = grep_pattern('cuda/*.m', '\n%#codegen -args');
    files = regexp(lines, '(\w+.m):', 'tokens');
    
    for j=1:length(files)
        file = files{j}{1};
        m2c('-cuda', '-mex', '-noinf', '-O', '-cppflags', ...
            ['{''-I' mspackroot '/include''}'], varargin{:}, file);
    end
catch ME
    cd(curpath);
    rethrow(ME)
end

cd(curpath);
