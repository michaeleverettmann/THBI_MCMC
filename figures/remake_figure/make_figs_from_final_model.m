% brb2022.06.03. Look at how the model pdfs will change if we use more/less
% chains. Pick some station where you have results you like. Then loop
% through iNumChains, allowing between 1 and 16 chains to be utilized in
% estimating the final model. Then make some plots. We could also expand to
% look at how the estimated data changes. 

%% Startup. Path definitions
clc; clear; 
run('../../a0_STARTUP_BAYES.m');

% resdir_data = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/O53A_TA_dat1/many_sw_authors'; 
% resdir_data = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/R54A_TA_dat1/add_sediment_try2'; 
resdir_data = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_cnsi/S57A_TA_dat1/layerise_normal'; 
% resdir_data = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/S57A_TA_dat1/sage_gage'; 
resdir_fig = '/Users/brennanbrunsvik/Documents/temp/remake_thbi_figures/S57A_TA_dat1/layerise_normal/original_parent_pulse'; 
prior_path = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/prior.mat' ; 

%% Load files from the inversion. 
mat_files = ls([resdir_data '/*.mat']); 
mat_files = regexp(mat_files,'\n','split'); 
fprintf('\n Loading prior \n'); 
load(prior_path); % Load this before other things - prior.mat also has par which might be outdated. 

% mat_files = {[resdir_data '/misfits_perchain_orig.mat'],...
%     [resdir_data '/allmodels_perchain_orig.mat'],...
%     [resdir_data '/trudata_USE.mat'],[resdir_data '/par.mat']}; % These variables are all that are totally required. 

for i_mat_file = 1:length(mat_files); 
    if isempty(mat_files{i_mat_file}); continue; end; % we get an empty value after last \n when splitting string into multiple cells. 
    fprintf('\nLoading %s\n',mat_files{i_mat_file}); 
    load(mat_files{i_mat_file}); 
end

% % % par.inv.datatypes = {'SW_Ray_phV_eks', 'SW_Ray_phV_dal', ...
% % %         'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}

misfits_perchain_orig = misfits_perchain; 
allmodels_perchain_orig = allmodels_perchain; 
% [par, ~] = update_bayes_inv_parms(par, 'add_sediment_try2'); 

%%% Only temporary things here! Things to make your specific files run. 
% par.inv.datatypes = {'SW_Ray_phV_eks', 'SW_Ray_phV_dal', ...
%     'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; warning('Updated par inv datatypes'); 
% par.inv.datatypes = {'RF_Sp_ccp'}; 
% par.inv.datatypes = {'SW_Ray_phV_eks', 'SW_Ray_phV_dal', ...
%     'SW_Lov_phV', 'SW_HV'}; 

% par = update_bayes_inv_parms(par, 'RF_Sp_ccp_only'); 

%% Do some tests. This is modified from run_all_chains.m

par.res.resdir = resdir_fig; % For saving files to new location
par.inv.synthTest = false; % Get an error without this. 
mkdir(resdir_fig); 
    
goodChainManual = logical([ones(12,1)]); warning('brb2022.07.06: Setting good chains manual'); 
% goodChainManual(2:end,:)=false; 
% goodChainManual = logical([zeros(12,1)]); warning('brb2022.07.06: Setting good chains manual'); 
% goodChainManual(3)=true; 
% goodChainManual = []; 

[misfits_perchain,allmodels_perchain,goodchains,...
     misfits_perchain_original,...
     allmodels_perchain_original,...
     allmodels_collated] ...
     = c1_PROCESS_RESULTS( misfits_perchain_orig,allmodels_perchain_orig,par,1,[resdir_fig,'/modMisfits'],goodChainManual)

[ hypparm_trends ] = plot_HYPERPARAMETER_TRENDS( allmodels_perchain,[resdir_fig,'/hyperparmtrend.pdf'] );
plot_MISFIT_TRENDS(par,allmodels_perchain,misfits_perchain,resdir_fig );

posterior = c2_BUILD_POSTERIOR(allmodels_collated,par,par.res.zatdep);
plot_MODEL_SUMMARY(posterior,par,1,[resdir_fig,'/modMisfits.pdf']); 
plot_PRIORvsPOSTERIOR(prior,posterior,par,1,[resdir_fig,'/prior2posterior.pdf'])
fprintf('  > Plotting model suite\n')
[ suite_of_models ] = c3_BUILD_MODEL_SUITE(allmodels_collated,par );
plot_SUITE_of_MODELS( suite_of_models,posterior,1,[resdir_fig,'/suite_of_models.png'],[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
final_model = c4_FINAL_MODEL(posterior,allmodels_collated,par,1,[resdir_fig,'/final_model']);
plot_FINAL_MODEL( final_model,posterior,1,[resdir_fig,'/final_model.pdf'],true,[par.data.stadeets.Latitude,par.data.stadeets.Longitude]);
plot_HEATMAP_ALLMODELS(suite_of_models,final_model,par,1,[resdir_fig,'/heatmap_of_models.pdf']);
plot_HEATMAP_ALLMODELS_shallow(suite_of_models,final_model,par,1,[resdir_fig,'/heatmap_of_models_shallow.pdf']);

% % par.datprocess.HKappa = struct(              ...
% %                        'min_error', 0.002,           ... % Add this much "error" to h-kappa stacks (error of 0 can result in sigma inverting improperly)
% %                        'scale_error', 1,           ... % Multiply h-kappa error by this constant. Sigma needs to be scaled accordingly. If using 100, we can think of it like percent. 
% %                        'weightDistanceMax', 0,   ... % At start of burnin, gives 0 to 1 weight toward the (scaled) Euclidian distance from HKappa energy maximum. In otherwords, this tends toward disregarding the actual energy value, and pays attention to its position. 
% %                        'weightDistanceMin', 0); 
% par.forc.mindV = 0.05; warning('Changing propmat resolution');  
% par.forc.synthperiod = 2.5; warning('Changing propmat Period');  
% par.forc.mindV = 0.075; 
% par.datprocess.CCP.simple_parent_pulse = false; warning('Setting simple parent pulse = true'); 
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
% save([resdir_fig,'/final_predata'],'final_predata');

[ final_misfit ] = b4_CALC_MISFIT( trudata,final_predata,par,0, 'plotRFError',true );
[ final_log_likelihood,final_misfit ] = b5_CALC_LIKELIHOOD( final_misfit,trudata,final_model.hyperparms,par );
plot_TRUvsPRE( trudata,final_predata,1,[resdir_fig,'/final_true_vs_pred_data.pdf'], allmodels_collated);
% plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,1,[resdir_fig,'/final_true_vs_pred_data_wavs.png']);
plot_TRUvsPRE_WAVEFORMS( trudata,final_predata,1,[resdir_fig,'/final_true_vs_pred_data_wavs.png']);

% plot_TRUvsPRE_WAVEFORMS_indiv_figs(trudata,final_predata,1,[resdir,'/final_true_vs_pred_data_wavs.pdf']);
% save([resdir_fig,'/final_misfit'],'final_misfit');
plot_FIG2_FIT_MODEL( final_model,posterior,prior,par,1,[resdir_fig,'/fig2_FIT_MODEL.pdf'])

%Plot kernels for final model. 
plot_all_sensitivity_kernels(final_model,trudata,par,resdir_fig)