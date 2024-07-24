% close all 
% clear all % Cannot clear all, because then we allways reset network_manual and station_manual below. 
% These define what we are running. Make a list of all desired options. 

% % Loop over various synthetic tests. stamp defines the data used. station determines what the synthetic structure will be. 
% STAMP_all =          {'standard'       , 'all_no_hv'       , 'all_no_sp'       , 'standard'       , 'standard'         , };
% network_manual_all = {'testnwk'        , 'testwk'          , 'testnwk'         , 'testnwk'        , 'testnwk'          , }; 
% station_manual_all = {'cont_EProt-s1m1', 'cont_EProt-s1m1' , 'cont_EProt-s1m1' , 'cont_EProt-s1'  , 'cont_EProt-s1m1m2', }; %


STAMP_all = {'more_anis'};
network_manual_all = {'testnwk'}; 
station_manual_all = {'cont_EProt-s1m1m2'}; %

start_dir = pwd(); 
for istamp = 1:length(STAMP_all); 
    cd(start_dir); % Don't know if directories get messed up, but change to here just in case.
    STAMP = STAMP_all{istamp};
    network_manual = network_manual_all{istamp}; 
    station_manual = station_manual_all{istamp}; %


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
    % if ~ (exist('network_manual', 'var') && exist('station_manual', 'var')) ; 
    fprintf('\nReseting to %s.%s\n',network_manual,station_manual)
    % end
    run_params.sta = station_manual; % name of station
    run_params.nwk = network_manual; % name of network
    
    %% Run it
    fprintf('Starting synthetic test %s %s %s\n', network_manual, station_manual, STAMP)
    execute_run_all_chains; % Put "run_all_chains.m" in a function, so that maybe "clear" won't erase values like STAMP_all 
end

function execute_run_all_chains
    run_all_chains;% Put "run_all_chains.m" in a function, so that maybe "clear" won't erase values like STAMP_all 
end