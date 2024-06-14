clc; 
close all; 
clear; 
restoredefaultpath; 

%% Setup
% These define what model we are running. Make a list of all desired options. Check SYNTHETICS/bayes_inv_parms.m for more options. 
STAMP_all =          {'standard_temp'  };
network_manual_all = {'testnwk'        }; 
station_manual_all = {'simple_layers_1'}; %

start_dir = pwd(); 
for istamp = 1:length(STAMP_all); % Loop over each synthetic model. 

    % Handling of directories, paths, etc. 
    cd(start_dir); 
    STAMP = STAMP_all{istamp};
    network_manual = network_manual_all{istamp}; 
    station_manual = station_manual_all{istamp}; %

    run('../a0_STARTUP_BAYES.m') % Setup paths and stuff
    
    proj = struct('name', 'SYNTHETICS'); 
    paths = getPaths(); 
    proj.STAinversions = paths.STAinversions; 
    proj.dir = [fileparts(mfilename('fullpath')),'/'];
    save([proj.dir,'project_details.mat'],'proj')
    
    wd = pwd; 
    addpath(wd);
    cd(proj.dir); 
    
    %% specify details of this run
    generation = 0; % generation of solution and data processing
    gc = '';
    BWclust = '';
    onesta = '';
    
    %% put parameters in place 
    global run_params
    
    % you can adapt these (likely in some loop over stations etc.) to be appropriate for your dataset
    run_params.projname = proj.name; % from above
    run_params.gc = gc; % great circle distance of body wave data, if relevant
    run_params.BWclust = BWclust; % cluster of BW data, if relevant
    run_params.datN = generation; % processing iteration, if relevant
    run_params.STAMP = STAMP; % NEED - some identifier for this inversion run
    run_params.overwrite = 1; % do you want to overwrite previous results?
    run_params.sta = station_manual; % name of station
    run_params.nwk = network_manual; % name of network
    
    %% Run it
    fprintf('Starting synthetic test %s %s %s\n', network_manual, station_manual, STAMP)
    execute_run_all_chains; % Put "run_all_chains.m" in a function so that "clear" won't erase values like STAMP_all 
end

function execute_run_all_chains
    run_all_chains; % Put "run_all_chains.m" in a function, so that maybe "clear" won't erase values like STAMP_all 
end