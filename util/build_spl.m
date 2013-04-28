function build_spl( varargin)
%build_spl Build script for SpaLab

% First, compile test scripts
lines = grep_pattern( '*.m', '\n%#codegen -args');
files = regexp( lines, '(\w+.m):', 'tokens');

for j=1:length(files)
    file = files{j}{1};
    compile('-acc -noinf -O', file, varargin{:});
end
