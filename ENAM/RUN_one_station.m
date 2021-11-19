% function RUN_one_station(net, sta) % Can comment out the function part of this and run as normal script, if you define net and sta. But have to move run_param
%

%% setup data and stuff. 
% evdata1_database; disp('Temporary - downloading station info in RUN_ALL_STAS.m') 
% RUN_prep_data; warning('bb2021.11.15. Doing RUN_prep_data in RUN_one_station - but this only needs to be ran once, not for every station.') % Gets event data and data paths on this computer. 
run('../a0_STARTUP_BAYES.m')


%%%%% Important! Must define network and station before runnig this! 
% network_manual = 'US' ; disp('resetting to US.CEH')
% station_manual = 'CEH'; disp('resetting to US.CEH')
% network = input('Network? (use quotes):')
disp('Network and station') 
disp(network_manual) 
disp(station_manual)

% Paths stuff
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
STAMP = 'ENAM_trial';

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

tempSta = and(string(stainfo.nwk)==network_manual, string(stainfo.stas)==station_manual); 

tempStaInd = find(tempSta); 

%% ==================  LOOP OVER STATIONS IN DB  ================== 
% No, actually just find the right station and run it. 
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
    MASTER_par;
end
% end