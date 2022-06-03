% brb2022.06.03. Look at how the model pdfs will change if we use more/less
% chains. Pick some station where you have results you like. Then loop
% through iNumChains, allowing between 1 and 16 chains to be utilized in
% estimating the final model. Then make some plots. We could also expand to
% look at how the estimated data changes. 

%% Startup. Path definitions
clc; clear; 
run('../../a0_STARTUP_BAYES.m');

pathToStaResult = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/O53A_TA_dat1/all_002'; 
resdir_main = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/figures/number_of_iterations'; 
prior_path = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/prior.mat'; 

%% Load files from the inversion. 
mat_files = ls([pathToStaResult '/*.mat']); 
mat_files = regexp(mat_files,'\n','split'); 
for i_mat_file = 1:length(mat_files); 
    if isempty(mat_files{i_mat_file}); continue; end; % we get an empty value after last \n when splitting string into multiple cells. 
    fprintf('\nLoading %s\n',mat_files{i_mat_file}); 
    load(mat_files{i_mat_file}); 
end

fprintf('\n Loading prior \n'); 
load(prior_path); 

misfits_perchain_orig = misfits_perchain; 
allmodels_perchain_orig = allmodels_perchain; 

%% Do some tests. This is modified from MASTER_par.m
for iNumChains = [1, 2, 3, 4, 6, 9, 12, 16]; 
    fprintf('\n----------------------------------\nStarting iNumChains = %1.0f\n',iNumChains); 
    resdir = sprintf('%s/nchain_%1.0f', resdir_main, iNumChains); 
    par.res.resdir = resdir; % For saving files to new location
    par.inv.synthTest = false; % Get an error without this. 
    mkdir(resdir); 
    
    [misfits_perchain,allmodels_perchain,goodchains,...
         misfits_perchain_original,...
         allmodels_perchain_original,...
         allmodels_collated] ...
         = c1_PROCESS_RESULTS( misfits_perchain_orig,allmodels_perchain_orig,par,1,[resdir,'/modMisfits'],...
         [1:iNumChains]')

    posterior = c2_BUILD_POSTERIOR(allmodels_collated,par,par.res.zatdep);
    plot_MODEL_SUMMARY(posterior,par,1,[resdir,'/modMisfits.pdf']); 
    plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[resdir,'/prior2posterior.pdf'])
    fprintf('  > Plotting model suite\n')
    [ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par );
    plot_SUITE_of_MODELS( suite_of_models,posterior,1,[resdir,'/suite_of_models.png'],[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
    plot_HEATMAP_ALLMODELS(suite_of_models,par,1,[resdir,'/heatmap_of_models.pdf']);
    final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[resdir,'/final_model']);
    plot_FINAL_MODEL( final_model,posterior,1,[resdir,'/final_model.pdf'],true,[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
end
