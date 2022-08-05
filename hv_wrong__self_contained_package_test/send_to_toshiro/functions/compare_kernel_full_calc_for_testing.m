% function [Kbase, KbasePrev, predata, log_likelihood, misfit, Pm_prior, nchain] = ...
%     compare_kernel_full_calc(...
%         model, Kbase, predata, trudata, ID, par, fail_chain, ...
%         ii, ifpass, misfit, ptbnorm, log_likelihood, ifaccept); 
%    
    
% % % % Temp. Remake old model. 
% % % model0 = Kbase.modelk; 
% % % [Kbase0,predata0] = b7_KERNEL_RESET(model0,Kbase,predata,ID,ii,par,1);
% % % % [Kbase0.HV.HVr, Kbase.HV.HVr] % < were same 

% Interpolate onto old model basis. 
% % % %%% Let's do some shinanigains to model, the new one, and figure out what's
% % % %%% wrong. 
% % % parm_change = {'VS', 'VP', 'rho', 'Sanis', 'Panis'}; % 'z', 'z0', 
% % % model_copy = model; 
% % % for iparm = 1:length(parm_change); 
% % %     this_parm = parm_change{iparm}; 
% % %     model.(this_parm) = ...
% % %         linterp(model_copy.z, model_copy.(this_parm), model0.z); 
% % % end
% % % model.z = model0.z ; 
% % % model.z0= model0.z0; 
% % %     
    
% Interpolate onto new model basis. 
model0 = Kbase.modelk; 
% [Kbase0.HV.HVr, Kbase.HV.HVr] % < were same 
%%% Let's do some shinanigains to model, the new one, and figure out what's
%%% wrong. 
parm_change = {'VS', 'VP', 'rho', 'Sanis', 'Panis'}; % 'z', 'z0', 
model_copy = model0; 
% Use just one of these interpolation loops. 
% % % for iparm = 1:length(parm_change); 
% % %     this_parm = parm_change{iparm}; 
% % %     model0.(this_parm) = ...
% % %         linterp(model0.z, model0.(this_parm), model.z);     
% % % % %     %%% Make two models equivalent in some depth range
% % %     ngauss = 2; 
% % %     imean = 1; 
% % %     gwin = gausswin(ngauss); 
% % %     gwin(1:ceil(ngauss/2)) = 1; 
% % %     model.(this_parm)(imean:imean+ngauss-1) = ...
% % %         gwin .* model0.(this_parm)(imean:imean+ngauss-1) + ...
% % %         (1-gwin) .* model.(this_parm)(imean:imean+ngauss-1)
% % % end
% % % model0.z = model.z ; 
% % % model0.z0= model.z0; 
%%%% model0.Nz = length(model0.z); 

% % % % new_z = linspace(0, 300, 300)'; % Optional, can use at base of next for loop 
% % % new_z = [0; logspace(log10(0.1), log10(10), 50)'; [11:300]'];
% % % model0.z (model0.z (1:end-1)==model0.z (2:end)) = model0.z (model0.z (1:end-1)==model0.z (2:end)) - 0.01; 
% % % model.z  (model .z (1:end-1)==model .z (2:end)) = model .z (model .z (1:end-1)==model .z (2:end)) - 0.01; 
% % % model0.z0(model0.z0(1:end-1)==model0.z0(2:end)) = model0.z0(model0.z0(1:end-1)==model0.z0(2:end)) - 0.01; 
% % % model.z0 (model .z0(1:end-1)==model .z0(2:end)) = model .z0(model .z0(1:end-1)==model .z0(2:end)) - 0.01; 
% % % 
% % % for iparm = 1:length(parm_change); 
% % %     this_parm = parm_change{iparm}; 
% % %     model0.(this_parm) = ...
% % %             interp1(model0.z, model0.(this_parm), new_z, 'pchip');    
% % %     model.(this_parm) = ...
% % %         interp1(model.z, model.(this_parm), new_z, 'pchip');     
% % % end
% % % model0.z = new_z ; model0.z0 = new_z ; 
% % % model.z = new_z ; model.z0 = new_z; 
% % % model0.Nz = length(new_z); model.Nz = length(new_z); 

%%% Yet another fake model 1. 
model = model0; 
ind1 = 1; ind2 = 2; 
model.VS (ind1:ind2) = model.VS (ind1:ind2) + 1; 
model.VP (ind1:ind2) = model.VP (ind1:ind2) + 1; 
model.rho(ind1:ind2) = model.rho(ind1:ind2) + 1; 

%%%

figure(1); clf; hold on; plot(model0.VS, model0.z); plot(model.VS, model.z); set(gca, 'ydir', 'reverse'); 
Kbase.modelk = model0; 
[Kbase,predata] = b7_KERNEL_RESET(model0,Kbase,predata,ID,ii,par,1);
predata = b3_FORWARD_MODEL_SW_kernel(model,Kbase,par,predata );
% [model0.VS, model.VS, model0.z, model.z]
%%%

% [predata1,laymodel1] = b3__INIT_PREDATA(model,par,trudata,0 );
% [predata1,par] = b3_FORWARD_MODEL_BW(model,laymodel1,par,predata1,ID,0,predata);
% predata1 = b3_FORWARD_MODEL_RF_ccp(   model,laymodel1,par,predata1,ID,0 );
% predata1 = b3_FORWARD_MODEL_SW_kernel(model,Kbase0,par,predata1 );
% [predata.SW_HV.HVr, predata1.SW_HV.HVr] % < were same

[KbaseNew,predataNew] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,1);

if ifforwardfail(predataNew,par) % Sometimes kernel reset gives nans, but doesn't explicitly fail? 
    fail_chain=fail_chain+1; ifpass=0;
    error('Failed making new kernels! Fail_chain %.0f, ii=%1.0f\n',fail_chain,ii)
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

    lik_drop_mak_fig = -inf; % if lik_drop_mak_fig < 5; warning('Making extra HV kernel plots'); end;  % Temporary. Should make this very high. 
    if (abs(logL_diff) > lik_drop_mak_fig) && strcmp(fn,'SW_HV') ; % This shouldn't happen. Make some plots. 
        % Compare HV from previous Kbase, new Kbase, and with the current model based on old Kbase + kernels. 
        
        LW = 1.5; 
        figure(3); clf; hold on; set(gcf, 'pos', [-1215 280 894 638]); 
        tiledlayout(1,3,'TileSpacing','compact'); 

        %%% HV ratio
        nexttile(1, [1,1]); cla; hold on; box on; 
        set(gca, 'ydir', 'reverse', 'LineWidth', 1.5);  
        xlim([0.6, 1.2]); 
        ylim([min(trudata.SW_HV.periods-1.5), max(trudata.SW_HV.periods+5)]); 
        set(gca, 'yscale', 'log'); 
        grid off; grid on;  % Log scale requires turning grid off before on
        
        plot(trudata.SW_HV.HVr, trudata.SW_HV.periods, ...
            'DisplayName', 'Goal', 'LineWidth', LW*3/2, ...
            'Color', 'Green'); 
        
        sw_hv = Kbase.HV; 
        plot(sw_hv.HVr, sw_hv.periods, ...
            'DisplayName', 'Old HV (full calc)', 'LineWidth', LW, ...
            'Color', 'k');
        
        sw_hv = KbaseNew.HV; 
        plot(sw_hv.HVr, sw_hv.periods, ...
            'DisplayName', 'New HV (full calc)', 'LineWidth', LW, ...
            'Color', 'blue');
        
        sw_hv = predata.(fn); 
        plot(sw_hv.HVr, sw_hv.periods, ...
            'DisplayName', 'New HV (from kernels)', 'LineWidth', LW, ...
            'Color', 'r');
       
        legend('Location', 'Best'); 
        xlabel('H/V ratio'); ylabel('Period');
        
        mismatch = rms(predata.(fn).HVr - KbaseNew.HV.HVr); % Mistmach between kernel and new full forward model
        dC_est   = rms(predata.(fn).HVr - Kbase   .HV.HVr); % dC predicted from kernel
        dC_full  = rms(Kbase   .HV.HVr  - KbaseNew.HV.HVr); % dC based on full calculations 
        %%% End HV ratio

        %%% Models
        nexttile(2, [1,2]); hold on; box on; 
        set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
        plot([Kbase.modelk.VS, Kbase.modelk.VP], Kbase.modelk.z, ...
            'DisplayName', 'Previous kernel model', 'LineWidth', LW, ...
            'Color', 'k');
        plot([KbaseNew.modelk.VS,KbaseNew.modelk.VP], KbaseNew.modelk.z, ...
            'DisplayName', 'New kernel model', 'LineWidth', LW, ...
            'Color', 'blue');
        leg=legend('Location', 'Best'); 
        xlabel('V'); ylabel('Depth (km)'); 
        ylim([-5, max(Kbase.modelk.z)]); 
        %%% End models

        sgtitle(sprintf(...
            ['HV kernel versus full calculation.\nIteration=%6.0f. ',...
            'Model pertubation norm = %3.2f. ',...
            'Change in log-likelihood = %3.5f\n',...
            'RMS of: Kernel error: %1.4f -- dC kernel: %1.4f -- dC full model: %1.4f.'],... 
            par.ii, ptbnorm, logL_diff,...
            mismatch, dC_est, dC_full)); 
        

        exportgraphics(gcf, sprintf('%s/%s_ii_%1.0f_likelihood_drop_hv.pdf',...
            par.res.resdir, par.res.chainstr, ii ) ); % %!%! Temporary 

        %%% Plot sensitivity kernels from old and new Kbase. 
        % If they are very different then we need to update kernels more often. 
        plot_sensitivity_kernel_hv(Kbase.HV.KHV   , 'dat', Kbase.HV   ,... 
            'model',Kbase.modelk   ,...
            'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_old.pdf', ...
            par.res.resdir, par.res.chainstr, par.ii) ); % Old kernels
        
        plot_sensitivity_kernel_hv(KbaseNew.HV.KHV, 'dat', KbaseNew.HV, ...
            'model',KbaseNew.modelk,...
            'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_new.pdf',...
            par.res.resdir, par.res.chainstr, par.ii) ); % New kernels
        %%% End sensitivity kernels. 
    end
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

% end