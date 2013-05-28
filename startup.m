% Startup script of palSpa

if ~exist('./util/build_psp.m', 'file')
    error('You must run the startup script in the palSpa''s root directory.');
end

% Add palCoder
if ~exist('palc', 'file')
    try
        run('../palCoder/startup.m');
    catch %#ok<CTCH>
        error(['Could not find palCoder in the path. Please run palCoder''s ' ...
            'startup script manually.']);
    end
end

% Add palACC
if ~exist('pACC_begin_parallel', 'file')
    try
        run('../palACC/startup.m');
    catch %#ok<CTCH>
        error(['Could not find palACC in the path. Please run palACC''s  ' ...
            'startup script manually.']);
    end
end

% Add palMPI
if ~exist('pMPI_Init', 'file')
    try
        run('../palMPI/startup.m');
    catch %#ok<CTCH>
        error(['Could not find palMPI in the path. Please run palMPI''s  ' ...
            'startup script manually.']);
    end
end

% Add palLina
if ~exist('build_pla', 'file')
    try
        run('../palLina/startup.m');
    catch %#ok<CTCH>
        error(['Could not find palLina in the path. Please run palLina''s  ' ...
            'startup script manually.']);
    end
end

addpath(pwd); %#ok<*MCAP>
addpath([pwd '/util']); %#ok<*MCAP>

if ~isinmpi
    try
        build_psp
    catch
        warning('Could not build palSpa.');
    end
end
