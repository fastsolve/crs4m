% Startup script of SpaLab

if ~exist('./util/build_sl.m', 'file')
    error('You must run the startup script in the SpaLab''s root directory.');
end

% Add M2C
if ~exist('m2c', 'file')
    try
        run('../M2C/startup.m');
    catch %#ok<CTCH>
        error(['Could not find M2C in the path. Please run M2C''s ' ...
            'startup script manually.']);
    end
end

% Add MACC
if ~exist('MACC_begin_parallel', 'file')
    try
        run('../MACC/startup.m');
    catch %#ok<CTCH>
        error(['Could not find MACC in the path. Please run MACC''s  ' ...
            'startup script manually.']);
    end
end

% Add MMPI
if ~exist('MMPI_Init', 'file')
    try
        run('../MMPI/startup.m');
    catch %#ok<CTCH>
        error(['Could not find MMPI in the path. Please run MMPI''s  ' ...
            'startup script manually.']);
    end
end

% Add LinaLab
if ~exist('build_lal', 'file')
    try
        run('../LinaLab/startup.m');
    catch %#ok<CTCH>
        error(['Could not find LinaLab in the path. Please run LinaLab''s  ' ...
            'startup script manually.']);
    end
end

addpath(pwd); %#ok<*MCAP>
addpath([pwd '/util']); %#ok<*MCAP>

if ~exist(['dotprod.' mexext], 'file')
    build_sl
end
