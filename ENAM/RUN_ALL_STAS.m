close all
clear all
%% Setup
run('../a0_STARTUP_BAYES.m')
evdata1_database; disp('Temporary - downloading station info in RUN_ALL_STAS.m') 

paths = getPaths(); 
proj = struct('name','ENAM');
proj.dir = [paths.THBIpath '/' proj.name];
wd = pwd; addpath(wd);
cd(proj.dir);

%% load project, station, and request details and request details
try
    load([proj.dir,'/project_details.mat']);
catch
    run([proj.dir,'/project_details.m']);
end
load([proj.infodir,'stations.mat']); % bb2021.09.28 Get this using evdata1_database.m

%% specify details of this run
generation = 1; % generation of solution and data processing
gc = 1;
BWclust = 1;
STAMP = 'ENAM_trial_US_CEH';

onesta = '';

overwrite = true;

%% put parameters in place for running all stations
global run_params

run_params.projname = proj.name;
run_params.gc = gc;
run_params.BWclust = BWclust;
run_params.datN = generation;
run_params.STAMP = STAMP;
run_params.overwrite = overwrite;

tempSta = and(string({stations_IRIS.NetworkCode})=='US', string({stations_IRIS.StationCode})=='CEH'); 
tempStaInd = find(tempSta); 

%% ==================  LOOP OVER STATIONS IN DB  ================== 
% for is = 1:stainfo.nstas % done 30-end, need to do before and after...
for is = tempStaInd; disp('Only doing 1 station right now!!!!')
    
    if exist('onesta') && ~isempty(onesta)
        if ~strcmp(stainfo.stas{is},onesta), continue; end
    end
    
    fprintf('\n'); for i=1:3, for j=1:40, fprintf('**'); end; fprintf('*\n'); end; fprintf('\n');
    
    fprintf('STATION: %s\n',stainfo.stas{is})
    fprintf('NETWORK: %s\n\n',stainfo.nwk{is})
    run_params.sta = stainfo.stas{is};
    run_params.nwk = stainfo.nwk{is};
    
    % do the work (and make all the Mineos files) in a workdir
    if exist('workdir','dir')~=7, mkdir('workdir'); end
%     cd('workdir')
    cd(proj.dir)
    
    execute_MASTER_par
    
    cd(proj.dir)
end

function execute_MASTER_par
    global run_params
%     try
    MASTER_par;
%     catch
%          1;
%     end
end