%Test script for palSpa

% First, compile test scripts
if exist( './test_psp.m', 'file')
    lines = grep_pattern( '../*.m', '\n%!test');
else
    lines = grep_pattern( '*.m', '\n%!test');
end
files = regexp( lines, '([\.\/\\\w]+.m):', 'tokens');

for j=1:length(files)
    file = files{j}{1};
    fprintf('Testing %s... ', file);
    test_mcode( file);
end
