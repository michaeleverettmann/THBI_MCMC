function [Kbase, KbasePrev, predata, log_likelihood, misfit, Pm_prior, nchain, ifpass] = ...
    compare_kernel_full_calc(...
        model, Kbase, predata, trudata, ID, par, fail_chain, ...
        ii, ifpass, misfit, ptbnorm, log_likelihood, ifaccept); 
   

[KbaseNew,predataNew] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,1);

if ifforwardfail(predataNew,par) % Sometimes kernel reset gives nans, but doesn't explicitly fail? 
    fail_chain=fail_chain+1; 
    ifpass=0;
    error('Failed making new kernels! Fail_chain %.0f, ii=%1.0f\n',fail_chain,ii)
else; 
    ifpass = 1; %brb2022.07.15 Not really sure if this is appropriate. 
end

KbasePrev = Kbase; % Now that we knot b7 worked, we can define Kbaseprev to be the last Kbase we had. 

[log_likelihoodNew,misfitNew] = b8_LIKELIHOOD_RESET(par,predataNew,trudata,KbaseNew,model.datahparm);

%%% Evaluate how much the data has changed from using kernels to using the full calculation. 
if ifaccept; % Only do comparison if last model was accepted. Otherwise predata doesn't correspond to model. 
rst_str = '\n>>>>>>>>>>>>>>>>>>>\nChanges between kernel predicted and full forward modelled error: \n'; 
rst_str_l = ''; 
fns = fieldnames(misfitNew.E2); 
for ifn = 1:length(fns); 
    
    %%% Get info on change in model/data fit. 
    fn = fns{ifn}; 
    E2New   = misfitNew.E2.(fn); 
    E2      = misfit   .E2.(fn); 
    p_change = (E2New - E2) / E2; 
    
    logLNew = misfitNew.logL_indivdat.(fn); 
    logL    = misfit   .logL_indivdat.(fn); 
    logL_diff = logLNew - logL; 
    %%% End model/data fit. 
    
    %%% building strings that collectively show each data fit. 
    rst_str = [rst_str ...
        sprintf('E2     %10.3f -> %10.3f. | dE/E %10.3f : %s\\n', ...
        E2, E2New, p_change, fn)]; 
    
    rst_str_l = [rst_str_l ...
        sprintf('l_like %10.3f -> %10.3f. | Diff  %10.3f : %s\\n', ...
        logL, logLNew, logL_diff, fn)]; 
    %%% End strings for each data fit. 
% % % 
% % %     lik_drop_mak_fig = -inf; % if lik_drop_mak_fig < 5; warning('Making extra HV kernel plots'); end;  % Temporary. Should make this very high. 
% % %     if (abs(logL_diff) > lik_drop_mak_fig) && strcmp(fn,'SW_HV') ; % This shouldn't happen. Make some plots. 
% % %         % Compare HV from previous Kbase, new Kbase, and with the current model based on old Kbase + kernels. 
% % %         LW = 1.5; 
% % %         figure(3); clf; hold on; set(gcf, 'pos', [-1215 280 894 638]); 
% % %         tiledlayout(1,3,'TileSpacing','compact'); 
% % % 
% % %         %%% HV ratio
% % %         nexttile(1, [1,1]); cla; hold on; box on; 
% % %         set(gca, 'ydir', 'reverse', 'LineWidth', 1.5);  
% % %         xlim([0.6, 1.2]); 
% % %         ylim([min(trudata.SW_HV.periods-1.5), max(trudata.SW_HV.periods+5)]); 
% % %         set(gca, 'yscale', 'log'); 
% % %         grid off; grid on;  % Log scale requires turning grid off before on
% % %         
% % %         plot(trudata.SW_HV.HVr, trudata.SW_HV.periods, ...
% % %             'DisplayName', 'Goal', 'LineWidth', LW*3/2, ...
% % %             'Color', 'Green'); 
% % %         
% % %         sw_hv = Kbase.HV; 
% % %         plot(sw_hv.HVr, sw_hv.periods, ...
% % %             'DisplayName', 'Old HV (full calc)', 'LineWidth', LW, ...
% % %             'Color', 'k');
% % %         
% % %         sw_hv = KbaseNew.HV; 
% % %         plot(sw_hv.HVr, sw_hv.periods, ...
% % %             'DisplayName', 'New HV (full calc)', 'LineWidth', LW, ...
% % %             'Color', 'blue');
% % %         
% % %         sw_hv = predata.(fn); 
% % %         plot(sw_hv.HVr, sw_hv.periods, ...
% % %             'DisplayName', 'New HV (from kernels)', 'LineWidth', LW, ...
% % %             'Color', 'r');
% % %        
% % %         legend('Location', 'Best'); 
% % %         xlabel('H/V ratio'); ylabel('Period');
% % %         
% % %         mismatch = rms(predata.(fn).HVr - KbaseNew.HV.HVr); % Mistmach between kernel and new full forward model
% % %         dC_est   = rms(predata.(fn).HVr - Kbase   .HV.HVr); % dC predicted from kernel
% % %         dC_full  = rms(Kbase   .HV.HVr  - KbaseNew.HV.HVr); % dC based on full calculations 
% % %         %%% End HV ratio
% % % 
% % %         %%% Models
% % %         nexttile(2, [1,2]); hold on; box on; 
% % %         set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
% % %         plot([Kbase.modelk.VS, Kbase.modelk.VP], Kbase.modelk.z, ...
% % %             'DisplayName', 'Previous kernel model', 'LineWidth', LW, ...
% % %             'Color', 'k');
% % %         plot([KbaseNew.modelk.VS,KbaseNew.modelk.VP], KbaseNew.modelk.z, ...
% % %             'DisplayName', 'New kernel model', 'LineWidth', LW, ...
% % %             'Color', 'blue');
% % %         leg=legend('Location', 'Best'); 
% % %         xlabel('V'); ylabel('Depth (km)'); 
% % %         ylim([-5, max(Kbase.modelk.z)]); 
% % %         %%% End models
% % % 
% % %         sgtitle(sprintf(...
% % %             ['HV kernel versus full calculation.\nIteration=%6.0f. ',...
% % %             'Model pertubation norm = %3.2f. ',...
% % %             'Change in log-likelihood = %3.5f\n',...
% % %             'RMS of: Kernel error: %1.4f -- dC kernel: %1.4f -- dC full model: %1.4f.'],... 
% % %             par.ii, ptbnorm, logL_diff,...
% % %             mismatch, dC_est, dC_full)); 
% % %         
% % % 
% % %         exportgraphics(gcf, sprintf('%s/%s_ii_%1.0f_likelihood_drop_hv.pdf',...
% % %             par.res.resdir, par.res.chainstr, ii ) ); % %!%! Temporary 
% % % 
% % %         %%% Plot sensitivity kernels from old and new Kbase. 
% % %         % If they are very different then we need to update kernels more often. 
% % %         plot_sensitivity_kernel_hv(Kbase.HV.KHV   , 'dat', Kbase.HV   ,... 
% % %             'model',Kbase.modelk   ,...
% % %             'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_old.pdf', ...
% % %             par.res.resdir, par.res.chainstr, par.ii) ); % Old kernels
% % %         
% % %         plot_sensitivity_kernel_hv(KbaseNew.HV.KHV, 'dat', KbaseNew.HV, ...
% % %             'model',KbaseNew.modelk,...
% % %             'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_new.pdf',...
% % %             par.res.resdir, par.res.chainstr, par.ii) ); % New kernels
% % %         %%% End sensitivity kernels.
% % %         
% % %         %%%
% % % % % %         plot_sensitivity_kernel_hv_grid(Kbase.HV.KHV   , 'dat', Kbase.HV   ,... 
% % % % % %             'model',Kbase.modelk   ,...
% % % % % %             'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_old_grid.pdf', ...
% % % % % %             par.res.resdir, par.res.chainstr, par.ii) ); % Old kernels
% % % % % %         
% % %         %%%
% % %         
% % %         save( sprintf('%s/%s_ii_%1.0f_likelihood_drop_hv_workspace.mat',...
% % %             par.res.resdir, par.res.chainstr, ii ) ); 
% % %     if (abs(logL_diff) > lik_drop_mak_fig) && strcmp(fn,'SW_Ray_phV_eks') ; % % Temporary
% % %         % Compare HV from previous Kbase, new Kbase, and with the current model based on old Kbase + kernels. 
% % %         LW = 1.5; 
% % %         figure(4); clf; hold on; set(gcf, 'pos', [-1215 280 894 638]); 
% % %         tiledlayout(1,3,'TileSpacing','compact'); 
% % %         
% % %         ray = Kbase.Ray ; 
% % %         rayNew = KbaseNew.Ray; 
% % % 
% % %         %%% Rayleigh wave phase velocities. 
% % %         nexttile(1, [1,1]); cla; hold on; box on; 
% % %         set(gca, 'ydir', 'reverse', 'LineWidth', 1.5);  
% % %         xlim([2.5, 4.5]); 
% % %         ylim([min(ray.periods-1.5), max(ray.periods+5)]); 
% % %         set(gca, 'yscale', 'log'); 
% % %         grid off; grid on;  % Log scale requires turning grid off before on
% % %                
% % %         plot(ray.phV, ray.periods, ...
% % %             'DisplayName', 'Old phv (full calc)', 'LineWidth', LW, ...
% % %             'Color', 'k');
% % %         
% % %         plot(rayNew.phV, rayNew.periods, ...
% % %             'DisplayName', 'New phv (full calc)', 'LineWidth', LW, ...
% % %             'Color', 'blue');
% % %         
% % %         for ifn2 = 1:length(fns); 
% % %             fn2 = fns{ifn2}; 
% % %             if ~contains(fn2, 'SW_Ray_phV'); 
% % %                 continue
% % %             end
% % %             plot(predata.(fn2).phV, predata.(fn2).periods, ...
% % %                 'DisplayName', 'New phV (from kernels)', 'LineWidth', LW, ...
% % %                 'Color', 'r');
% % %         end
% % %         legend('Location', 'Best'); 
% % %         xlabel('Rayleigh phase velocity (km/s)'); ylabel('Period');
% % %     end
end

d_log = log_likelihoodNew-log_likelihood; 

tot_change = sprintf('l_like %10.3f -> %10.3f. | Diff: %10.3f : Combined',...
    log_likelihood, log_likelihoodNew, d_log);

ptb_str = sprintf('ptbnorm was %1.5f\n',ptbnorm); 

if d_log < -5;
	tot_change = [tot_change '\n\nLog likelihood dropped a lot after resetting kernels!\n',...
     '    Check which data type is the problem and \n    make sure the kernels are resetting often enough!']; 
end

% Combine all the strings from above and finally print them. 
fprintf([rst_str '\n' ...
         rst_str_l '\n\n' ...
         ptb_str '\n' ...
         tot_change '\n<<<<<<<<<<<<<<<<<<<<\n']); 
end
%%% End evaluating influence of kernels. 

% Accept the new misfits and model.  
misfit           = misfitNew; 
log_likelihood   = log_likelihoodNew; 
Kbase            = KbaseNew; 
predata          = predataNew; 

Pm_prior = calc_Pm_prior(model,par);
nchain = 0;

end