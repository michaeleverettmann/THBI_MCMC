function [Kbase, KbasePrev, predata, log_likelihood, misfit, Pm_prior, nchain, ifpass] = ...
    compare_kernel_full_calc(...
        model, Kbase, predata, trudata, ID, par, fail_chain, ...
        ii, ifpassIn, misfit, ptbnorm, log_likelihood, ifaccept); 
ifpass = ifpassIn; % Sometimes there is a weird thing taking in and returning a variable. 

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