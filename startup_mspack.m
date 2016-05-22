% Startup script of MSPACK

msproot = fileparts(which('startup_mspack.m'));

addpath(msproot); %#ok<*MCAP>
addpath([msproot '/util']); %#ok<*MCAP>

if ~exist('./util/build_msp.m', 'file')
    error('You must run the startup script in the MSPACK''s root directory.');
end

% Add M2C
if ~exist('m2c', 'file')
    try
        run('../M2C/startup.m');
    catch %#ok<CTCH>
        error(['Could not find m2c in the path. Please run M2C''s ' ...
            'startup script.']);
    end

    % Add MACC_
    if ~exist('MACC_begin_parallel', 'file')
        try
            run('../MACC_/startup.m');
        catch %#ok<CTCH>
            error(['Could not find MACC_ in the path. Please run MACC_''s  ' ...
                   'startup script.']);
        end
    end

    % Add MMPI
    if ~exist('mpi_Init', 'file')
        try
            run('../MMPI/startup.m');
        catch %#ok<CTCH>
            error(['Could not find MMPI in the path. Please run MMPI''s  ' ...
                   'startup script.']);
        end
    end
    
    if usejava('jvm')
        build_msp;
    end
end
