% Startup script of MSPACK

msproot = fileparts(which('startup_mspack.m'));

addpath(msproot);
addpath([msproot '/cuda']);
addpath([msproot '/cuda/sparse']);
addpath([msproot '/Vec']);
addpath([msproot '/Mat']);
addpath([msproot '/Mat/BSR']);
addpath([msproot '/Mat/CRS']);
addpath([msproot '/KSP']);
addpath([msproot '/PC']);
addpath([msproot '/util']);

% Add M2C
if ~exist('m2c', 'file')
    run([fileparts(msproot) '/M2C/startup.m']);
    run([fileparts(mspoot) '../MOMP/startup.m']);
    run([fileparts(mspoot) '../MCUDA/startup.m']);
    
    if ~exist(['crs_create.' mexext], 'file')
        fprintf(1, 'To build MSPACK, run command build_mspack in MATLAB/Octave.\n');
    end
end
