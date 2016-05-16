% Startup script of MSPACK

addpath(pwd); %#ok<*MCAP>
addpath([pwd '/util']); %#ok<*MCAP>

if ~exist('./util/build_msp.m', 'file')
    error('You must run the startup script in the MSPACK''s root directory.');
end

% Add M2C
if ~exist('m2c', 'file')
    try
        if exist('../M2C','dir')==7
            run('../M2C/startup.m');
        elseif exist('../M2C_dist','dir')==7
            run('../M2C_dist/startup.m');
        else
            run('../CodeGen/startup.m');
        end
    catch %#ok<CTCH>
        error('Could not find m2c in the path. Please add M2C to your path.');
    end

    if usejava('jvm')
        build_msp
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
if ~exist('pMPI_Init', 'file')
    try
        run('../MMPI/startup.m');
    catch %#ok<CTCH>
        error(['Could not find MMPI in the path. Please run MMPI''s  ' ...
            'startup script manually.']);
    end
end

% Add LinaLab
if ~exist('build_pla', 'file')
    try
        run('../LinaLab/startup.m');
    catch %#ok<CTCH>
        error(['Could not find LinaLab in the path. Please run LinaLab''s  ' ...
            'startup script manually.']);
    end
end
