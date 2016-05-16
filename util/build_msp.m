function build_msp(varargin)
%build_msp Build script for MSP

% First, compile test scripts
lines = grep_pattern('crs*.m', '\n%#codegen -args');
files = regexp(lines, '(\w+.m):', 'tokens');

for j=1:length(files)
    file = files{j}{1};
    compile('-omp', '-noinf', '-O', varargin{:}, file);
end
