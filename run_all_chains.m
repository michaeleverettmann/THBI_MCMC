% This script does: 
% 1. Load data and prepare for inversion. 
% 2. Do a parallel loop over each MCMC chain. 
% 3. Plot and store results. 

clear % clear all to make sure we use values below. No longer needed maybe? 
close all

time_start_inversion = now() * 24; % In hours. 

global run_params
paths = getPaths(); 

if isempty(run_params) % Some default values. 
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

run([paths.THBIpath,'/a0_STARTUP_BAYES']); % Crucial: setup all paths. 
load('project_details'); 
addpath([proj.dir,'matguts/']); 

%% PARMS
run parms/bayes_inv_parms
[par, inv] = update_bayes_inv_parms(par, STAMP); % Modify inversion parameters depending on your STAMP. Useful for doing various synthetic tests. 

if any(strcmp(projname,{'SYNTHETICS','matlab_to_mineos_vJBR', 'example_github'})); 

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

%% Get some directories ready. 
% Switch to execution folder, to make synthetic data. 
prev_dir = pwd(); 
cd(paths.ramDrive); % Execute everything from a folder in ram for major speedup. 
mkdir([nwk '_' sta]); cd([nwk '_' sta]); % Go to station specific folder to keep things clean.  

%% ========================  LOAD + PREP DATA  ========================
[trudata,par] = a2_LOAD_DATA(par, 'nwk', nwk, 'sta', sta);
plot_TRU_WAVEFORMS(trudata);
cd(prev_dir); % Not sure if need to change back to whatever directory we were previously in, but I will for safety. 
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

% Estimate prior PDFs. For comparison to the inverted PDFs. 
if ~par.mod.force_no_new_prior && redoprior
    fprintf('  > Building prior distribution from %.0f runs\n',max([par.inv.niter,1e5]))
    prior = a3_BUILD_EMPIRICAL_PRIOR(par,max([par.inv.niter,1e5]),14,par.res.zatdep);
    plot_MODEL_SUMMARY(prior,par,1,[resdir,'/prior_fig.pdf']);
    save([proj.dir,'/prior'],'prior','par');
end

%% ---------------------------- INITIATE PARALLEL -------------------------

% ===== Prepare for parallel pool =====
delete(gcp('nocreate'));
myCluster = parcluster('local');
maxWorkers = myCluster.NumWorkers; 

% Matlab stores parallel job metadata in some folder. If running many parpools in parallel, they might overwrite each others metadata, giving a corrupt file error message. Give each parallel pool its own folder. 
fprintf('\nChecking storage availability on /tmp/ folder\n\n'); 
!du -h /tmp/
!df -h /tmp/
JobStorageLocation = sprintf('%s/job_info_scratch/%s_%s_%s',...
    paths.THBIpath, nwk, sta, STAMP); 
JobStorageLocation_old = [JobStorageLocation '_old']; 
% JobStorageLocation = sprintf('/tmp/mcmcthbi_%s_%s_%s/',... % Use tmp for local disk. % https://researchcomputing.princeton.edu/support/knowledge-base/matlab
mkdir(JobStorageLocation);
mkdir(JobStorageLocation_old); 
system(sprintf('mv %s/* %s',JobStorageLocation, JobStorageLocation_old)); 
pause_time = 3; 
fprintf('Pausing for %1.3f seconds to allow job storage folder to be created\n',pause_time); 
pause(pause_time);
myCluster.JobStorageLocation = JobStorageLocation; 

% % Use the following two lines if you want more parallel chains than you have cores. 
% myCluster.NumWorkers = max(maxWorkers, par.inv.nchains); 
% maxWorkers = myCluster.NumWorkers; 

fprintf('Starting parpool with JobStorageLocation %s\n',JobStorageLocation); 
tstartppool = tic; 
try 
    myCluster.parpool(min([par.inv.nchains,maxWorkers])); % TODOcomp May need to change based on computer. I set so that we only use as many workers as the local parallel profile will allow.  
catch error_parpool
    fprintf(['Problem opening parpool! Error displayed below. ',...
        '\n >>>>>>>>>>>>>>>> \n %s \n <<<<<<<<<<<<<<< \n',...
        'Trying to start parpool one more time!'],...
        getReport(error_parpool) ); 
    delete(gcp('nocreate')); 
    system(sprintf('mv %s/* %s',JobStorageLocation, JobStorageLocation_old)); % Remove possibly corrupted .mat files. 
    pause_time = rand(1)*3; 
    fprintf('Pausing for %1.3f seconds before starting pool\n',pause_time); 
    pause(pause_time);
    myCluster.parpool(min([par.inv.nchains,maxWorkers])); % TODOcomp May need to change based on computer. I set so that we only use as many workers as the local parallel profile will allow.  
end
fprintf('Time to start parpool: \n'); 
toc(tstartppool)

    
TD = parallel.pool.Constant(trudata); PR = parallel.pool.Constant(par);
% (((( If not parallel: ))))
% TD(1). Value = trudata; par.inv.verbose = 0;

%% START DIFFERENT MARKOV CHAINS IN PARALLEL
model0_perchain = cell(par.inv.nchains,1);
misfits_perchain = cell(par.inv.nchains,1);
allmodels_perchain = cell(par.inv.nchains,1);
SWs_perchain = cell(par.inv.nchains,1);

t = now;

% Option to profile code. 
if profileRun;  % Start profiling parfor iterations, where most calculations happen. 
    mpiprofile on; 
end

cd(paths.ramDrive); % Execute everything from a folder in ram for major speedup. 
mkdir([nwk '_' sta]); cd([nwk '_' sta]); % Go to station specific folder to keep things clean . TODO just to cd once. 

% fprintf('Debug: if this split definition if the Seizmo one and not Matlab builtin, the code will break:\n'); 
% help split
% fprintf('END Some info on split:\n'); 

%% ========================================================================
%% ========================================================================
fprintf('\n ============== STARTING CHAIN(S) ==============\n')
%% ========================================================================
%% ========================================================================

% Probably do not run in parallel until you know it works in serial for loop. 
% parfor iii = 1:par.inv.nchains % TODO Will need to change between for and parfor, depending on circumstance. for is needed if wanting to do debuging. 
parfor iii = 1:par.inv.nchains % TODO Will need to change between for and parfor, depending on circumstance. for is needed if wanting to do debuging. 
    par = PR.Value; 
    trudata = TD.Value; 
    [ model0_perchain{iii}, misfits_perchain{iii},...
        allmodels_perchain{iii}, SWs_perchain{iii} ]=...
        run_one_chain(par, trudata, nwk, sta, iii)
end 

[ram_copy_stats] = ram_to_HD(paths, resdir, nwk, sta); % Copy final results from ram to hard disk. Also remove the ram drive for this station, and change directory to main results folder for this statino. 

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
c0_SAVE_OUTPUT(resdir,misfits_perchain,allmodels_perchain,par);

fprintf('Duration of entire run: %.0f s\n',(now-t)*86400)
plot_invtime(misfits_perchain,[resdir,'/invTime.pdf']);

%% Process results
fprintf('  > Processing results\n')
[misfits_perchain,allmodels_perchain,goodchains,...
 misfits_perchain_original,...
 allmodels_perchain_original,...
 allmodels_collated] ...
 = c1_PROCESS_RESULTS( misfits_perchain,allmodels_perchain,par,1,[resdir,'/modMisfits']);

plot_corrplot(par, allmodels_collated, [resdir,'/xi_correlation_matrix.pdf'] ); 

[ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain,[resdir,'/hyperparmtrend.pdf'] );
plot_MISFIT_TRENDS(par,allmodels_perchain,misfits_perchain_original,resdir );
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

%% Save some things
fprintf('  > Saving misfits, allmods, posterior, model suite\n')
save([resdir,'/posterior'],'posterior');
save([resdir,'/goodchains'],'goodchains');
try
save([resdir,'/SWs_pred'],'SWs_perchain');
end

take_lots_of_space = false; 
if take_lots_of_space; 
    save([resdir,'/misfits_perchain'],'misfits_perchain'); % Already saved the "orig" version of this. 
    save([resdir,'/allmodels_perchain'],'allmodels_perchain'); % Already saved the "orig" version of this. Not savingthis now. We save the original version earlier. This takes a lot of space, so I'm prioritizing the unmodified misfits which we can re-use to get back to this point without running the inversion. 
%     save([resdir,'/allmods_collated'],'allmodels_collated');
    save([resdir,'/mod_suite'],'suite_of_models');
end


%% Final interpolated model with errors
% If you have to few iterations/chains, there aren't enough models. You may get an error in c4_FINAL_MODEL. 
final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[resdir,'/final_model']);
plot_FINAL_MODEL( final_model,posterior,1,[resdir,'/final_model.pdf'],true,[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
plot_HEATMAP_ALLMODELS(suite_of_models,final_model,par,1,[resdir,'/heatmap_of_models.pdf']);
plot_HEATMAP_ALLMODELS_shallow(suite_of_models,final_model,par,1,[resdir,'/heatmap_of_models_shallow.pdf']);

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
plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,posterior,par,1,[resdir,'/final_true_vs_pred_data_wavs.png']);
% plot_TRUvsPRE_WAVEFORMS_indiv_figs(trudata,final_predata,1,[resdir,'/final_true_vs_pred_data_wavs.pdf']);
save([resdir,'/final_misfit'],'final_misfit');


plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,1,[resdir,'/fig2_FIT_MODEL.pdf']);

%Plot kernels for final model. 
plot_all_sensitivity_kernels(final_model,trudata,par,resdir)

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
if ~ take_lots_of_space % Clear some space while we are at it and saving things. 
    fprintf('\nRemoving velocity profile files (allmodels should still have those)\n'); 
    system(sprintf('rm ./*%s.%s*vel_profile',nwk,sta)); 
    fprintf('\nRemoving invState files\n'); 
    system(sprintf('rm ./%s.%s*invState',nwk,sta)); 
end

time_end_inversion = now() * 24; % In hours. 
fprintf('\nInversion took a total of %2.5f hours. \n',time_end_inversion-time_start_inversion); 
