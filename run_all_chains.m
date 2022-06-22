clear % clear all to make sure we use values below
close all

global run_params
paths = getPaths(); 

if isempty(run_params)
    projname = 'AK'; % SYNTHETICS, LAB_tests, or US, for now
    sta = 'BPAW';
    nwk = 'AK';
    gc = 'all'; % will search for gcarcs +/-3 of this value; % 'all' to do all gcs
    BWclust = 0;
    datN = 30; % generation of data processing
    STAMP=[sta,datestr(now,'_yyyymmddHHMM_pll')];
    overwrite = true;
    cd AK
else
    projname = run_params.projname;
    sta = run_params.sta;
    nwk = run_params.nwk;
    gc = run_params.gc;
    BWclust = run_params.BWclust;
    datN = run_params.datN;
    STAMP = run_params.STAMP;
    overwrite = run_params.overwrite;
end



notes = ['']; 

%% ------------------------- START -------------------------
global projdir TRUEmodel
projdir = [paths.THBIpath,'/',projname,'/'];
cd(projdir);

% load paths which might change from one to another computer
% pathsToChange = loadPathsToChange(); bb2021.09.13 Started this but it's going to be a PITA

run([paths.THBIpath,'/a0_STARTUP_BAYES']);
load('project_details'); %TODO_STATION_NETWORK bb2021.11.12
addpath([proj.dir,'matguts/']);



%% PARMS
run parms/bayes_inv_parms
[par, inv] = update_bayes_inv_parms(par, STAMP); % Modify this function to make different tests. 

% if exist('external_data_types', 'var') && external_data_types; 
%     % For this block, external_data_types and this_data_type have to be
%     % created before running run_all_chains.m 
%     % They should only serve to manually modify par.inv.datatypes after
%     % running bayes_inv_parms. 
%     par.inv.datatypes = this_data_type; % bb2021.12.07 Thinking of way to loop over each data type in solo. 
% end

if strcmp(projname,'SYNTHETICS')
% % %     bb2021.12.07 Removing the re-definitions of sta. I want synthetic test charactaristic, not noise, to define sta. 
% % %     if isfield(par.synth,'noisetype') && strcmp(par.synth.noisetype,'real'); 
% % %         sta=['SYNTH_',sta]; 
% % %     else
% % %         sta = 'SYNTH'; 
% % %     end
    par.stadeets = struct('sta',sta','nwk',nwk'); 

	% noise details, if "real"
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise

elseif strcmp(projname,'LAB_tests')
	zsed = 0;
	zmoh = 45;
	zlab = 130;
	wlab = 10;
	flab = 0.05;
    par.synth.model = struct('zsed',zsed,'zmoh',zmoh,'zlab',zlab,'wlab',wlab,'flab',flab);
	dtps = {'BW_Ps','BW_Sp','BW_Sp_lo','BW_Ps_lo','SW_Ray_phV','SW_Lov_phV','SW_HV'};

	% noise details, if "real"
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise

	% naming convention
	dtpstr='_';
	if any(strcmp(dtps,'BW_Ps')), dtpstr=[dtpstr,'Ps']; end
	if any(strcmp(dtps,'BW_Ps_lo')), dtpstr=[dtpstr,'Pslo']; end
	if any(strcmp(dtps,'BW_Ps_cms')), dtpstr=[dtpstr,'Pscms']; end
	if any(strcmp(dtps,'BW_Sp')), dtpstr=[dtpstr,'Sp']; end
	if any(strcmp(dtps,'BW_Sp_lo')), dtpstr=[dtpstr,'Splo']; end
	if any(strcmp(dtps,'SW_Ray_phV')), dtpstr=[dtpstr,'Ray']; end
	if any(strcmp(dtps,'SW_Lov_phV')), dtpstr=[dtpstr,'Lov']; end
	if any(strcmp(dtps,'SW_HV')), dtpstr=[dtpstr,'HV']; end

	sta = ['LAB_s',num2str(zsed),'_m',num2str(zmoh),'_z',num2str(zlab),'_w',num2str(wlab),'_f',num2str(100*flab),dtpstr];
    nwk = 'LAB_test';
end

if strcmp(projname,'SYNTHETICS') || strcmp(projname,'LAB_tests')
    par.synth.noise_sta_deets = struct('datadir',['/Volumes/data/THBI/US/STAsinv/',noisesta,'_dat20/'],...
                         'sta',noisesta,'nwk',noisenwk,'gc',noisegcarcs,'noiseup',noiseup,'noiseshape',noiseshape);
end

par.inv.BWclust = BWclust;
ifsavedat = false;

%% get saving things ready
par.proj = proj;
avardir = sprintf('%s%s_%s_dat%.0f/',par.proj.STAinversions,sta,nwk,datN);
resdir = [avardir,STAMP];
if ~exist(resdir,'dir'), try mkdir(resdir); catch, error('Looks like no path to output directory - is it mounted?'); end, end
fid = fopen([resdir,'/notes.txt'],'w'); fprintf(fid,notes); fclose(fid);

par.data = struct('stadeets',struct('sta',sta,'nwk',nwk,'Latitude',[],'Longitude',[]),...
                  'gc',gc,'datN',datN,'avardir',avardir);

par.res.STAMP = STAMP;
par.res.resdir= resdir;
par.res = orderfields(par.res,{'STAMP','resdir','zatdep'});

par_ORIG = par;

% save par
if exist([resdir,'/misfits_perchain_orig.mat'],'file')==2
%     yn=input('Station already has results - overwrite? (y/[n]) ','s');
%     if ~strcmp(yn,'y'), return; end
    fprintf('Station already has results...')
    if ~overwrite, fprintf(' SKIPPING\n'); return; else, fprintf(' OVERWRITING\n'); end
end
save([resdir,'/par'],'par');

% save par and copy the parms.m file to the results directory
eval(sprintf('! cp parms/bayes_inv_parms.m %s',resdir))

allpdytp = parse_dtype_all(par);

%% ========================  LOAD + PREP DATA  ========================
[trudata,par] = a2_LOAD_DATA(par);
plot_TRU_WAVEFORMS(trudata);
% check_data(trudata,par)

%% ===========================  PRIOR  ===========================
% see if prior exists already
if exist([proj.dir,'prior.mat'],'file') 
    a = load([proj.dir,'prior.mat']);
    % only do prior if par has changed
    redoprior = par_equiv_prior(par,a.par);
    prior = a.prior;
else
    redoprior = true;
end

if redoprior
    fprintf('  > Building prior distribution from %.0f runs\n',max([par.inv.niter,1e5]))
    prior = a3_BUILD_EMPIRICAL_PRIOR(par,max([par.inv.niter,1e5]),14,par.res.zatdep);
    plot_MODEL_SUMMARY(prior,par,1,[resdir,'/prior_fig.pdf']);
    save([proj.dir,'/prior'],'prior','par');
end

%% ---------------------------- INITIATE ----------------------------

% ===== Prepare for parallel pool =====
delete(gcp('nocreate'));
myCluster = parcluster('local');
maxWorkers = myCluster.NumWorkers; 
% % Use the following two lines if you want more parallel chains than you have cores. 
% myCluster.NumWorkers = max(maxWorkers, par.inv.nchains); 
% maxWorkers = myCluster.NumWorkers; 

parpool(min([par.inv.nchains,maxWorkers])); % TODOcomp May need to change based on computer. I set so that we only use as many workers as the local parallel profile will allow.  
TD = parallel.pool.Constant(trudata); PR = parallel.pool.Constant(par);
% (((( If not parallel: ))))
% TD(1). Value = trudata; par.inv.verbose = 0;

%% START DIFFERENT MARKOV CHAINS IN PARALLEL
model0_perchain = cell(par.inv.nchains,1);
misfits_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);
SWs_perchain = cell(par.inv.nchains,1);

%% ========================================================================
%% ========================================================================
fprintf('\n ============== STARTING CHAIN(S) ==============\n')
%% ========================================================================
%% ========================================================================
t = now;
% mkdir([resdir,'/chainout']);
% parfor iii = 1:par.inv.nchains

if profileRun;  % Start profiling parfor iterations, where most calculations happen. 
    mpiprofile on; 
end

% mainDir = [paths.execPath '/' nwk '_' sta]; % Keep track of where the main folder is, where we want to return after changing directory back from ram drive. 
% TODO might cause problems to change to new directory because of prior.mat (which loads from absolute directory though, maybe ok) and project_details.mat) which might be different for different stations? 
% if ~ exist(mainDir); mkdir(mainDir); end % This is where we will cd to for final processing, and save final results. using exist here is ok, because we only do it once per stations.  NOTE don't need if exists, but it's REALLY important to not accidentally overwrite the whole folder. 
cd(paths.ramDrive); % Execute everything from a folder in ram for major speedup. 
mkdir([nwk '_' sta]); cd([nwk '_' sta]); % Go to station specific folder to keep things clean . TODO just to cd once. 



%% % % % % % parfor iii = 1:par.inv.nchains % TODO Will need to change between for and parfor, depending on circumstance. for is needed if wanting to do debuging. 
for iii = 1:par.inv.nchains % TODO Will need to change between for and parfor, depending on circumstance. for is needed if wanting to do debuging. 
warning('BB2021.11.22 Not in parallel!!!')
par = PR.Value; 
trudata = TD.Value; 

[ model0_perchain{iii}, misfits_perchain{iii},...
    allmodels_perchain{iii}, SWs_perchain{iii} ]=...
    run_one_chain(par, trudata, nwk, sta, iii)

end % parfor loop
%%
[ram_copy_stats] = ram_to_HD(paths, resdir, nwk, sta); % Copy final results from ram to hard disk. 
[~,duRam] = system(sprintf('du -h %s',paths.ramDrive)); 
fprintf('\nDisk usage of ram drive after moving files off of it: \n%s\n', duRam )
cd(resdir); % Get back out of ram and go to stations hard drive folder 

if profileRun; % Get results from profiling. 
    mpiprofile off; 
    mpiStats = mpiprofile('info'); 
    save(['mpiProfileData_' nwk '_' sta], 'mpiStats');  % Save profile results. Can transfer from HPC and bring to local computer for viewing. 
    fprintf('\n\nSaved mpiProfileData....mat to %s\n\n',resdir); % Use this to see the results: load('mpiProfileData'); mpiprofile('viewer', mpiStats);
end 

if par.inv.niter > 2000
	delete(gcp('nocreate'));
else
    warning('brb2022.03.16 Not deleting parallel pool, niter<2000 so I assume we are debugging. ')
end

%% ========================================================================
%% ========================================================================
%% ----------------------- End loop on chains  ----------------------------
%% ========================================================================
%% ========================================================================
c0_SAVE_OUTPUT(resdir,misfits_perchain,allmodels_perchain);

fprintf('Duration of entire run: %.0f s\n',(now-t)*86400)
plot_invtime(misfits_perchain,[resdir,'/invTime.pdf']);

%% Process results
fprintf('  > Processing results\n')
% misfits_perchain = misfits_perchain_original;
% allmodels_perchain = allmodels_perchain_original;
[misfits_perchain,allmodels_perchain,goodchains,...
 misfits_perchain_original,...
 allmodels_perchain_original,...
 allmodels_collated] ...
 = c1_PROCESS_RESULTS( misfits_perchain,allmodels_perchain,par,1,[resdir,'/modMisfits']);

[ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain,[resdir,'/hyperparmtrend.pdf'] );
plot_KNOT_TRENDS( allmodels_perchain,par,[resdir,'/knottrends']  )

posterior = c2_BUILD_POSTERIOR(allmodels_collated,par,par.res.zatdep);

fprintf('  > Plotting posterior\n')
plot_MODEL_SUMMARY(posterior,par,1,[resdir,'/posterior.pdf'])

fprintf('  > Plotting prior vs. posterior\n')
plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[resdir,'/prior2posterior.pdf'])
% plot_P2P_recovery(prior,posterior,TRUEmodel,par,1,[resdir,'/mparm_recovery_p2p.pdf'])

fprintf('  > Plotting model suite\n')
[ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par );
plot_SUITE_of_MODELS( suite_of_models,posterior,1,[resdir,'/suite_of_models.png'],[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
plot_HEATMAP_ALLMODELS(suite_of_models,par,1,[resdir,'/heatmap_of_models.pdf']);

%% Save some things
fprintf('  > Saving misfits, allmods, posterior, model suite\n')
save([resdir,'/misfits_perchain'],'misfits_perchain');
save([resdir,'/allmodels_perchain'],'allmodels_perchain');
save([resdir,'/posterior'],'posterior');
% save([resdir,'/allmods_collated'],'allmodels_collated');
save([resdir,'/mod_suite'],'suite_of_models');
save([resdir,'/goodchains'],'goodchains');
try
save([resdir,'/SWs_pred'],'SWs_perchain');
end

%% Final interpolated model with errors
% If you have to few iterations/chains, there aren't enough models, and you
% will get an error in c4_FINAL_MODEL. 
final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[resdir,'/final_model']);
plot_FINAL_MODEL( final_model,posterior,1,[resdir,'/final_model.pdf'],true,[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);

%% predict data with the final model, and calculate the error!
[ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,trudata,posterior );

% distribute data for different processing (e.g. _lo, _cms)
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    pdt = parse_dtype( dtype );
    if strcmp(pdt{1},'BW') && (~strcmp(pdt{3},'def') || ~strcmp(pdt{4},'def'))
        if any(strcmp(par.inv.datatypes,['BW_',pdt{2}])) % only if there IS a standard!
            disp(['replacing ',dtype,' with ',[pdt{1},'_',pdt{2}]])
            final_predata.(dtype) = final_predata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
        end
    end
end
% window, filter data
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ final_predata ] = predat_process( final_predata,dtype,par);
end
save([resdir,'/final_predata'],'final_predata');


[ final_misfit ] = b4_CALC_MISFIT( trudata,final_predata,par,0, 'plotRFError',true );
[ final_log_likelihood,final_misfit ] = b5_CALC_LIKELIHOOD( final_misfit,trudata,final_model.hyperparms,par );
plot_TRUvsPRE( trudata,final_predata,1,[resdir,'/final_true_vs_pred_data.pdf'], allmodels_collated);
plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,1,[resdir,'/final_true_vs_pred_data_wavs.png']);
% plot_TRUvsPRE_WAVEFORMS_indiv_figs(trudata,final_predata,1,[resdir,'/final_true_vs_pred_data_wavs.pdf']);
save([resdir,'/final_misfit'],'final_misfit');


plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,1,[resdir,'/fig2_FIT_MODEL.pdf']);

% did we save the data?
if ifsavedat
    savedat.gdmods = find([allmodels.bestmods]');
    savedat.gdmods(savedat.gdmods==0) = [];
    save([resdir,'/savedat'],'savedat');
% plot_DATAFITS(trudata,savedat,par,1)
    plot_FIG1_FIT_DATA( trudata,savedat,par,1,[resdir,'/fig1_FIT_DATA.pdf'])
end

% clear('TRUEmodel')
% return
