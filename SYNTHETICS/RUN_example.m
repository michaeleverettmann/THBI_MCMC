close all
clear all
%% Setup
run('../a0_STARTUP_BAYES.m')

proj = struct('name', 'SYNTHETICS'); % bb2021.08.04 changed from EXAMPLE because I don't have the example data files. %,'EXAMPLE');
paths = getPaths(); 
proj.STAinversions = paths.STAinversions; 
proj.dir = [fileparts(mfilename('fullpath')),'/'];
% proj.dir = paths.STAinversions; 
proj.STAinversions = paths.STAinversions; ; % [proj.dir,'inversion_results/'];
save([proj.dir,'project_details.mat'],'proj')


wd = pwd; addpath(wd);
cd(proj.dir);


%% specify details of this run
generation = 0; % generation of solution and data processing
gc = '';
BWclust = '';
STAMP = 'examplerun';
% STAMP = 'adding_sediment_pt1'; 

onesta = '';

%% put parameters in place 
global run_params

% you can obviously adapt these (likely in some loop over stations etc.) 
% to be appropriate for your dataset
run_params.projname = proj.name; % from above
run_params.gc = gc; % great circle distance of body wave data, if relevant
run_params.BWclust = BWclust; % cluster of BW data, if relevant
run_params.datN = generation; % processing iteration, if relevant
run_params.STAMP = STAMP; % NEED - some identifier for this inversion run
run_params.overwrite = 1; % do you want to overwrite previous results?
if ~ (exist('network_manual', 'var') && exist('station_manual', 'var')) ; 
    network_manual = 'testnwk'; 
    station_manual = 'crat_2mld'; 
    fprintf('\nReseting to %s.%s\n',network_manual,station_manual)
end
run_params.sta = station_manual; % name of station
run_params.nwk = network_manual; % name of network

%% Run it
run_all_chains;
