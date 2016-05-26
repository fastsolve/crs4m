%Test script for MSPACK

mspackroot = fileparts(which('startup_mspack'));
curpath = pwd;
cd(mspackroot);

try
    % First, compile test scripts
    lines = [grep_pattern('Mat/*.m', '\n%!test') ...
        grep_pattern('Vec/*.m', '\n%!test') ...
        grep_pattern('KSP/*.m', '\n%!test') ...
        grep_pattern('PC/*.m', '\n%!test')];
    files = regexp(lines, '([\.\/\\\w]+.m):', 'tokens');
    
    for j=1:length(files)
        file = files{j}{1};
        fprintf('Testing %s... ', file);
        test_mcode(file);
    end
catch ME
    cd(curpath);
    rethrow(ME)
end

cd(curpath);
