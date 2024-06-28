%% Compare the results of using the older version of Mineos (Zach's, was using before 2024, and Josh's, using after). 
clc; clear; 
run('../../a0_STARTUP_BAYES.m');

version = 'new'; % Old version is Zach's. New is Josh's. 
if strcmp(version, 'old'); 
    addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_mineos'); 
    rmpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/transition_to_russell_mineos/new_code'); 
    rmpath('/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/MINEOS_synthetics/run_MINEOS'); 
elseif strcmp(version, 'new'); 
    addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/transition_to_russell_mineos/new_code'); 
    addpath('/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/MINEOS_synthetics/run_MINEOS'); 
    rmpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_mineos'); 
end

resdir_data = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_collate/CEH_US_dat1/standard'; 
resdir_fig = [pwd() '/fig_out_' version '_mineos']; 
mkdir(resdir_fig); 
prior_path  = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/prior.mat' ; 
mkdir('working_directory'); 
cd('working_directory'); 

%% Load model files from an earlier inversion. 
mat_files = ls([resdir_data '/*.mat']); 
mat_files = regexp(mat_files,'\n','split'); 
fprintf('\n Loading prior \n'); 
load(prior_path); % Load this before other things - prior.mat also has par which might be outdated. 

for i_mat_file = 1:length(mat_files); 
    if isempty(mat_files{i_mat_file}); continue; end; % we get an empty value after last \n when splitting string into multiple cells. 
    fprintf('\nLoading %s\n',mat_files{i_mat_file}); 
    load(mat_files{i_mat_file}); 
end

par.res.resdir = resdir_fig; % For saving files to new location
par.synth.propmat_or_telewavesim = 'propmat'; % This is expected in newer version of mcmc. 
delete final_predata % Make sure we are recalculating the pre_data

%% Forward modelling and making plots. 
plot_FINAL_MODEL( final_model,posterior,1,[resdir_fig,'/final_model.pdf'],true,[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
[ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,trudata,posterior );
plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,posterior,par,1,[resdir_fig,'/final_true_vs_pred_data_wavs.png']);
kbase = plot_all_sensitivity_kernels(final_model,trudata,par,resdir_fig) % Get kernels for final model. 
save([resdir_fig, '/dat_and_kbase.mat'], 'final_predata', 'kbase'); % Make and save kernels. 