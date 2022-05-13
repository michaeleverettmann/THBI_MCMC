%% specify details of this run
generation = 1; % generation of solution and data processing
gc = 1;
BWclust = 1;
onesta = '';
overwrite = true;

% STAMP = 'ENAM_trial';
% STAMP_all = {...
%       'adding_sediment_pt1',...
%     'h_kappa_tests/degree_of_freedom/4',...
%     'h_kappa_tests/degree_of_freedom/9',...
%     'h_kappa_tests/degree_of_freedom/12',...
%     'h_kappa_tests/degree_of_freedom/15',...
%     'h_kappa_tests/degree_of_freedom/19',... don
%     'h_kappa_tests/degree_of_freedom/25',...
%     'h_kappa_tests/degree_of_freedom/30',...
%     'h_kappa_tests/degree_of_freedom/40',...
%     'h_kappa_tests/h_kappa_solo/sigma_constant/0.3'
%     'h_kappa_tests/observe_convergence/invert_sigma/0.3', ...
%     'h_kappa_tests/observe_convergence/invert_sigma_false/0.3', ...
%     'h_kappa_tests/observe_convergence/invert_sigma/large_cooldown/0.3', ...
%     'delayed_rejection_WRONG'
%     'real_all_data',...
%     'burnin_push_toward_minimum',...
%     'quarter_root_error',...
%     'higher_min_sigma',...
%     'prior_sigma_0.2_min_sigma_0.1',...
%     'prior_sigma_1_min_sigma_0.5',...
%     'prior_sigma_1_min_sigma_0.1_err',...
%      'HKappa_001',...
%      'HKappa_002',...
%       'HKappa_003',...
%     'HKappa_004',...
%     'HKappa_005',...
%     'HKappa_006',...
%     'HKappa_007',...
%     'all_002',...
%     'all_no_hv',...
%     'prior_sigma_10_min_sigma_3',...
%     'SSA_2022_v2',...
%     'SW_Ray_phV_only',...  % { --- Start test to see influence of different data types independently
%     'SW_Lov_phV_only',...  % ...
%     'RF_Sp_ccp_only' ,...  % ...
%     'HKstack_P_only' ,...  :q!
%     'SW_HV_only'     ,...  % } --- End   test to see influence of different data types independently
%     'SW_Ray_phV_only_one_chain',...  % { --- Start test to see influence of different data types independently
%     'SW_Lov_phV_only_one_chain',...  % ...
%     'RF_Sp_ccp_only_one_chain' ,...  % ...
%     'HKstack_P_only_one_chain' ,...  :q!
%     'SW_HV_only_one_chain'     ,...  % } --- End   test to see influence of different data types independently
%     'all_demo'       ,...
%     'HK_fast'   ,...
%     'HK_faster',...
% }; % This determines which tests we want to run now. They will run sequentially, not in parallel (each station only has one ram folder. )

% Copy what you want to run here. 
STAMP_all = {...
    'HK_faster'     ,...  % } --- End   test to see influence of different data types independently
}; % This determines which tests we want to run now. They will run sequentially, not in parallel (each station only has one ram folder. )




if exist('STAMP', 'var'); % If we defined this already in a bash script, then only run that stamp. Otherwise, go through the list of stamps above. 
    STAMP_all = {STAMP}; 
end

for istamp = [1:length(STAMP_all)]; 
    STAMP = STAMP_all{istamp}; 

    %% setup data and stuff. 
    % evdata1_database; disp('Temporary - downloading station info in RUN_ALL_STAS.m') 
    % RUN_prep_data; warning('bb2021.11.15. Doing RUN_prep_data in RUN_one_station - but this only needs to be ran once, not for every station.') % Gets event data and data paths on this computer. 
    run('../a0_STARTUP_BAYES.m')


    %%%%% Important! Must define network and station before runnig this! 
    % If we have not defined network and station, use default US.CEH
    if ~ (exist('network_manual', 'var') && exist('station_manual', 'var')) ; 
        network_manual = 'TA'; 
        station_manual = 'O53A'; 
        fprintf('\nReseting to %s.%s\n',network_manual,station_manual)
    end
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
    catch e 
        fprintf('\n%s\n',getReport(e)); 
        run([proj.dir,'/project_details.m']);
    end
    load([proj.infodir,'stations.mat']); % bb2021.09.28 Get this using evdata1_database.m

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
    for is = tempStaInd; 

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
        cd(proj.dir)

        execute_MASTER_par

        cd(proj.dir)
    end

end

function execute_MASTER_par
%     try
    MASTER_par;
%     catch e
%         fprintf('Cant execute master_par.m. %s\n%s','(todo put stamp here)',getReport(e)) 
%     end
end
% end