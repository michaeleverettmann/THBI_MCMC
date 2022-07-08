function [Kbase, KbasePrev, predata, log_likelihood, misfit, Pm_prior, nchain] = ...
    compare_kernel_full_calc(...
        model, Kbase, predata, trudata, ID, par, fail_chain, ...
        ii, ifpass, misfit, ptbnorm, log_likelihood, ifaccept); 
   

[KbaseNew,predataNew] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,1);

if ifforwardfail(predata,par) % Sometimes kernel reset gives nans, but doesn't explicitly fail? 
    fail_chain=fail_chain+1; ifpass=0;
    error('Failed making new kernels! Fail_chain %.0f, ii=%1.0f\n',fail_chain,ii)
end
KbasePrev = Kbase; % Now that we knot b7 worked, we can define Kbaseprev to be the last Kbase we had. 

[log_likelihoodNew,misfitNew] = b8_LIKELIHOOD_RESET(par,predataNew,trudata,Kbase,model.datahparm);

%%% Evaluate how much the data has changed from using kernels to using the full calculation. 
if ifaccept; % Only do comparison if last model was accepted. Otherwise predata doesn't correspond to model. 
rst_str = '\n>>>>>>>>>>>>>>>>>>>\nChanges between kernel predicted and full forward modelled error: \n'; 
rst_str_l = ''; 
fns = fieldnames(misfitNew.E2); 
for ifn = 1:length(fns); 
    fn = fns{ifn}; 
    E2New   = misfitNew.E2.(fn); 
    E2      = misfit   .E2.(fn); 
    p_change = (E2New - E2) / E2; 
    
    logLNew = misfitNew.logL_indivdat.(fn); 
    logL    = misfit   .logL_indivdat.(fn); 
    logL_diff = logLNew - logL; 
    
    rst_str = [rst_str ...
        sprintf('E2     %10.3f -> %10.3f. | dE/E %10.3f : %s\\n', ...
        E2, E2New, p_change, fn)]; 
    
    rst_str_l = [rst_str_l ...
        sprintf('l_like %10.3f -> %10.3f. | Diff  %10.3f : %s\\n', ...
        logL, logLNew, logL_diff, fn)]; 

    lik_drop_mak_fig = -inf; % if lik_drop_mak_fig < 5; warning('Making extra HV kernel plots'); end;  % Temporary. Should make this very high. 
    if (abs(logL_diff) > lik_drop_mak_fig) && strcmp(fn,'SW_HV') ; 
        LW = 1.5; 
        figure(3); clf; hold on; set(gcf, 'pos', [-1221 458 894 391]); 
        tiledlayout(1,3,'TileSpacing','compact'); 

        nexttile(1, [1,1]); hold on; box on; 
        set(gca, 'ydir', 'reverse', 'LineWidth', 1.5);  
        sw_hv = predata.(fn); 
        plot(sw_hv.HVr, sw_hv.periods, ...
            'DisplayName', 'From kernels', 'LineWidth', LW, ...
            'Color', 'r');
        sw_hv = predataNew.(fn); 
        plot(sw_hv.HVr, sw_hv.periods, ...
            'DisplayName', 'Full calculation', 'LineWidth', LW, ...
            'Color', 'k');
        legend('Location', 'Best'); 
        xlabel('H/V ratio'); ylabel('Period');

        nexttile(2, [1,2]); hold on; box on; 
        set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
        plot([Kbase.modelk.VS, Kbase.modelk.VP], Kbase.modelk.z, ...
            'DisplayName', 'Previous kernel model', 'LineWidth', LW, ...
            'Color', 'blue');
        plot([KbaseNew.modelk.VS,KbaseNew.modelk.VP], KbaseNew.modelk.z, ...
            'DisplayName', 'New kernel model', 'LineWidth', LW, ...
            'Color', 'k');
        leg=legend('Location', 'Best'); 
        xlabel('V'); ylabel('Depth (km)'); 
        ylim([-5, max(Kbase.modelk.z)]); 

        sgtitle(sprintf(...
            ['HV kernel versus full calculation.\nIteration=%6.0f. ',...
            'Model pertubation norm = %3.2f\n ',...
            'Change in log-likelihood = %3.5f'],...
            par.ii, ptbnorm, logL_diff )); 
        

        exportgraphics(gcf, sprintf('%s/likelihood_drop_hv_%s_ii_%1.0f.pdf',...
            par.res.resdir, par.res.chainstr, ii) ); % %!%! Temporary 

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