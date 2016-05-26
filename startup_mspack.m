% Startup script of MSPACK

msproot = fileparts(which('startup_mspack.m'));

addpath(msproot); %#ok<*MCAP>
addpath([msproot '/Mat']); %#ok<*MCAP>
addpath([msproot '/Vec']); %#ok<*MCAP>
addpath([msproot '/KSP']); %#ok<*MCAP>
addpath([msproot '/PC']); %#ok<*MCAP>
addpath([msproot '/util']); %#ok<*MCAP>

if ~exist('./util/build_mspack.m', 'file')
    error('You must run the startup script in the MSPACK''s root directory.');
end

% Add M2C
if ~exist('m2c', 'file')
    run('../M2C/startup.m');
    run('../MOMP/startup.m');
    
    if ~exist(['crs_create.' mexext], 'file')
        fprintf(1, 'To build MSPACK, run command build_mspack in MATLAB/Octave.\n');
    end
end
