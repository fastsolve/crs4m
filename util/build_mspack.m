function build_mspack(varargin)
%build_msp Build script for MSP

mspackroot = fileparts(which('startup_mspack'));
curpath = pwd;
cd(mspackroot);

try
    % First, compile test scripts
    lines = [grep_pattern('Mat/*.m', '\n%#codegen -args') ...
        grep_pattern('Vec/*.m', '\n%#codegen -args') ...
        grep_pattern('KSP/*.m', '\n%#codegen -args') ...
        grep_pattern('PC/*.m', '\n%#codegen -args')];
    files = regexp(lines, '(\w+.m):', 'tokens');
    
    for j=1:length(files)
        file = files{j}{1};
        m2c('-omp', '-mex', '-noinf', '-O', '-cppflags', ...
            ['{''-I' mspackroot '/include''}'], varargin{:}, file);
    end
catch ME
    cd(curpath);
    rethrow(ME)
end

cd(curpath);
