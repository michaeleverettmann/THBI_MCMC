% par = PR.Value; 
% trudata = TD.Value

function [model0, misfits, allmodels, preSW ...
    ]=run_one_chain(par, trudata, nwk, sta, iii ); 

%% Some setup
paths = getPaths(); 

% Disable a bspline warning that doesn't seem to matter. Needs to be placed in parfor or else individual workers don't keep this warning off. ; 
warning('off', 'MATLAB:rankDeficientMatrix'); % This comes up when doing least squares inversion for spline weights. Be careful, the rankDeficientMatrix could be needed at another point in the inversion...    

% Time the inversion
absTimeIter = timeseries(zeros(par.inv.niter,1)); 
absTimeIter(1) = datetime(); 

% Folder/name for this chain. 
chainstr = [nwk '.' sta '_' mkchainstr(iii)];
par.res.chainstr = chainstr; 
par.res.chainID = mkchainstr(iii); % ID to keep track of things like where the velocity profiles are. brb2022.03.30
par.res.chainExecFold = [paths.ramDrive '/' nwk '_' sta '/' par.res.chainID]; % Explicitly determine executable folder, since we will use rm ./* and don't want to risk executing this from somewhere on the hard drive!!! brb2022.03.30
mkdir(par.res.chainExecFold); 
cd(par.res.chainExecFold);  % Speed up execution of code by reducing number of files in a chains folder (this should help in Matlab system() calls). brb2022.03.30

% Structure to keep track of accepting/rejecting models and the related conditions
accept_info = make_accept_info(par); 

%%% Diary file stuff
diaryFile = sprintf('diary_%s.txt', chainstr); % Use a diary file to keep track of parallel inversions seperately. %TODO_STATION_NETWORK bb2021.11.12
diary off; 
delete(diaryFile); 
pause(0.01); % Because delete seems to fully execute after turning diary on...
diary(diaryFile);
%%% End diary file stuff. 

%% Fail-safe to restart chain if there's a succession of failures
fail_chain=20;
fail_total = 0; 
iter_fail = []; % Keep track of all iterations where there was an error. Needs to be resisable array, and keep track of previous failures even after resetting a chain. 
while fail_chain>=20

%%% Initialize several variables. 
newK = false; % Maybe make true. 
SW_precise = []; 
KbasePrev = []; % When this is empty, we've only made kernels once. Can reset to this model if needed. 
laymodel1 = []; 
non_acceptk = 0; % How often have we rejected the current baseline model. 
predataPrev = []; 
timeStartIter = zeros(par.inv.niter,1); % Vector of starting times of each iteration. 
timeVeryStart = now(); % Starting time of inversion
tFact = (24 * 60 * 60); % Factor to convert from serial time to seconds. 
%%% End initializing several varibales. 


%% Prep posterior structure
[ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par);

%% initiate model
[ifpass, numFails, model0, model,...
    Pm_prior, par, Kbase] = ...
    initiate_model(...
        par, trudata, chainstr, fail_chain, iii); 
    
    
%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
fprintf('\n =========== STARTING ITERATIONS %s ===========\n',chainstr)
ii = 0;
time0 = now;

ptb = cell({});
nchain = 0;
fail_chain = 0; 
fail_reset = 0;
ifaccept=true;
if isfield(trudata,'SW_Ray')
    preSW = zeros(length(trudata.SW_Ray.periods),ceil(par.inv.niter./par.inv.saveperN));
end
% reset_likelihood;
log_likelihood = -Inf;
predata=[]; 
predat_save = []; 
misfit = [];
par.hkResetInfo = struct('timesReset',0, 'timesSWKernelsReset',0); % Structure to keep track of when to reset HK stacks. 

while ii < par.inv.niter
       
save_exit_if_broken = false; % if save_exit_if_broken; warning('Will save and exit chain if there is an error'); end 
if (save_exit_if_broken) && (fail_chain > 0); 
    save(sprintf('%s/chain_with_broken_something_%s.mat',...
        par.res.resdir,chainstr))
    fail_chain = -100; 
    warning('fail_chain > 0. EXITING INVERSION FOR CHAIN %s',chainstr); 
    warning('Turn off save_exit_if_broken if you dont want an error to stop the inversion. %s',chainstr); 
    break
end
    
ii = ii+1; par.ii = ii; % bb2022.10.18 Keep track of ii in par for easily plotting how inversion is changing with ii. 
accept_info(ii).fail_chain = fail_chain; 

% % % %%%
if fail_chain > 0; % Problem -- There are various ways of dealing with this.
%     fprintf('\nFail chain: %2.0f. ii = %6.0f',fail_chain, ii)
%%% TODO add in Kbase.numRewinds. If it gets to > 3 or something, go to
%%% prev Kbase.
    resetData = false; 
    failInfoStr = sprintf('\n%s ii=%1.0f Failure. fail_chain = %1.0f, fail_total = %1.0f. ',...
            chainstr, ii, fail_chain, fail_total); 
        
    %%% If early fail
    earlyFailNum = 75; 
    if ii < earlyFailNum; % If chain is giving problems this early, just reset it.
        warning('%sError within first %1.0f iterations. Reseting chain.\n', failInfoStr, earlyFailNum);  
%         warning([failInfoStr 'Error within first 10 iterations. Reseting chain.']); 
        fail_chain = 100; fail_total = fail_total + 1; break; 
    end 
    
    %%% If later fail
    if fail_chain > 15; 
        if ii - Kbase.itersave > 50; % Last Kbase is probably a decent model
            warning([failInfoStr 'High fail chain. Reset to last Kbase.']); 
            model = Kbase.modelk; ii = Kbase.itersave; % Rewind to last model    
            resetData = true; 
        else % Last Kbase might also be bad. Go back two models. 
            if ~isempty(KbasePrev); 
                warning([failInfoStr 'High fail chain. Reset to last KbasePrev (NOT Kbase).']); 
                model = KbasePrev.modelk; ii = KbasePrev.itersave; 
                resetData = true; 
            else; % We don't know if last Kbase is good, and we only have one Kbase. Just reset the chain: we must not be very far anyway.  
                warning([failInfoStr 'High fail chain after reset Kbase, but KbasePrev doesnt exist. Resetting chain.']); 
                fail_chain = 100; fail_total = fail_total + 1; break; 
            end
        end
    end
    
    fprintf('\nNow on iteration %1.0f',ii); 
    iter_fail(end+1) = ii; 
    
    if resetData; 
        predata = b3_FORWARD_MODEL_BW( model,laymodel,par,predata,ID,0,predataPrev); % brb2022.04.12 The arguments to forward_model_bw were in the wrong order. Probaly an old version of the code. 
        predata = b3_FORWARD_MODEL_RF_ccp( model,laymodel,par,predata,ID,0 );
        predata = b3_FORWARD_MODEL_SW_kernel( model,Kbase,par,predata );
        predataPrev = predata; % Keep track of last predata. To keep the previous complete HK stack. 
        % need to also reset likelihood and misfit to the new, precise data (likelihood may have been artificially high due to kernel forward  calc. approximation - if so, need to undo this, or chain will get stuck once we reset kernels).
        [log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,model.datahparm);
        Pm_prior = calc_Pm_prior(model,par);
        nchain = 0;
    end
    %!%! Code to break chain entirely if reset to many times.     
end
% % % %%%

% Maintenance in case folder is getting filled up during last iterations. 
[~,~]=clean_ram(ii,par,'clearRamInterval',200) % Prevent large amount of files from building up in ram, which can dramatically slow down matlab. If the ram is getting very full, then there is probably a problem!
if ~strcmp(pwd,par.res.chainExecFold); 
    warning('\nSomehow you left the ram folder! This can tremendously slow down code execute. Changing back (note that cd is slow)... Fix it. brb2022.04.02'); 
    cd(par.res.chainExecFold); 
end
% END folder maintenance

% Quickly check to make sure inversion isn't going way to slow. 
timeStartIter(ii) = (now() - timeVeryStart) * tFact;  % Seconds since the inversion very first started. 
littleIter=max(1,ii-100); smallTime = timeStartIter(littleIter); bigTime = timeStartIter(ii); 
avgTime = (bigTime-smallTime) / (ii-littleIter); % Average computational time over past min([1,100]) iterations. 
if ii ~= 1; absTimeIter(ii) = datetime(); end; % Dont replace first time. Want that in case fail_chain is a PITA
% END getting computational time. 


%% SAVE model every saveperN
if mod(ii,par.inv.saveperN)==0 && log_likelihood ~= -Inf
    [misfits,allmodels,savedat] = b9_SAVE_RESULT(ii,log_likelihood,misfit,model,Pm_prior,misfits,allmodels,predat_save,savedat,time0);
%     if isfield(trudata,'SW_Ray_phV')
%         preSW(:,misfits.Nstored) = predata.SW_Ray_phV.phV;
%     end
%     figure(22); hold on
%     plot([model.vpvs],[model.zmoh],'-ow')
%     figure(23);
%     plot(model.VS,model.z,'k')
end

%% SAVE inv state every Nsavestate iterations
if rem(ii,par.inv.Nsavestate)==0
    save_inv_state(par.res.resdir,chainstr,allmodels,misfits)
%     [ram_copy_stats] = ram_to_HD(paths, chainstr, mainDir, nwk, sta); % bb2021.12.07 this is time consuming if done often. Just do at end of inversion. Copy current results from ram to hard disk. 
end

try
%     if ii == 3; error('Fake'); end % One error at ii = 5
%     if (ii > 10) && (ii < 20); error('Fake'); end % 10 errors 20 to 30
%     if (ii > 80) && (ii < 90); error('Fake');  end % 10 errors 20 to 30
    
    %% Figure out if chain is too slow or unstable. 
    % TODO changed from 4 to 1 times par.inv.savepern. So fairly verbose. 
    if rem(ii,1*par.inv.saveperN)==0 || ii==1, fprintf(...
            'Sta %s nwk %s  Iteration %s %1.0f Average comp time (last 100 iter) %2.2f (all) %2.2f\n',...
            par.data.stadeets.sta,par.data.stadeets.nwk,chainstr,ii,...
            avgTime, timeStartIter(ii)/ii); end
    if par.inv.verbose, pause(0.05); end
    ifaccept=false;
    ifpass = false;
    newK = false; resetK = false;
    if fail_chain>19
        % if not enough saved in this chain, abort and restart
        if (ii - par.inv.burnin)/par.inv.saveperN < 200
            fprintf('\nBreak chain -- spot a1\n')
            break
        % if enough saved in chain, abort and keep the incomplete chain
        else
            fprintf('\nBreak chain -- spot a2\n')
            fail_chain = -fail_chain; break
        end
    end
    if fail_reset>5
        % if cannot reset kernels because current saved model is not viable
        if (ii - par.inv.burnin)/par.inv.saveperN < 200
            fprintf('\nBreak chain -- spot a3\n')
            fail_chain = 100;  break % high fail_chain will mean we restart chain
        % if enough saved in chain, abort and keep the incomplete chain
        else
            fprintf('\nBreak chain -- spot a4\n')
            fail_chain = -100; break
        end
    end
    
    %%% Kill chain if too slow. 
    if avgTime > 1; % Things are running really slow... why? Might want to restart, unless we are at the beginning of iterations. 
%         if ~feature('IsDebugMode'); % Don't break slow chains if we are debugging. 
        tScaleKill = 1; 
        tryKeeping = ((ii - par.inv.burnin)/par.inv.saveperN) >= 200; 
        toBreak = ( (avgTime > 20  * tScaleKill ) && (ii > 3   ) ) || ...
                  ( (avgTime > 7   * tScaleKill ) && (ii > 20  ) ) || ...
                  ( (avgTime > 3.5 * tScaleKill ) && (ii > 200 ) ) || ...
                  ( (avgTime > 2   * tScaleKill ) && (ii > 500 ) ) || ...
                  ( (avgTime > 1.2 * tScaleKill ) && (ii > 1000) ); % Things should be going smoothly by now and always take ~ .6 s. However, there might be a time where mineos runs a few times and makes a good chain seem slow: maybe averaging over just 100 iterations could make us reject an otherwise ok chain. 
        if toBreak && tryKeeping; 
            fprintf('\nExecution too slow. Aborting (and keeping chain) %s ii=%1.0f, avg time last 100 iter = %1.2f\n',...
                chainstr,ii,avgTime); 
            fail_chain = -100; break; 
        elseif toBreak; 
            fprintf('\nExecution to slow. Resetting %s ii=%1.0f, avg time last 100 iter = %1.2f\n',...
                chainstr,ii,avgTime);
            fail_chain =  100;  break; 
        end  
%         end
    end
    %%% End Figure out if chain is too slow or unstable. 
    %%
   
    % temperature - for perturbation scaling and likelihood increase
    temp = (par.inv.tempmax-1)*erfc(2*(ii-1)./par.inv.cooloff) + 1;
%     if round_level(temp,0.01)>1
%         if par.inv.verbose, fprintf('TEMPERATURE = %.2f\n',temp); end
%     end

    while ifpass == false % only keep calculating if model passes (otherwise save and move on)

%% ===========================  PERTURB  ===========================
    if ii==1; log_likelihood1 = -Inf; model1 = model; ptbnorm = nan; ...
            Pm_prior1k = Pm_prior; end
    if non_acceptk == 0; p_bd = 1; end; %!%!
    
    delay_reject_bool = true; if (mod(ii, 100)==0) && ~delay_reject_bool; warning('brb2022.04.12 non_acceptk=0. NO DELAYED REJECTION.'); end; % brb2022.04.12. 
    if ~ delay_reject_bool; 
        non_acceptk = 0; 
    end
    
    % Perturb the model. Perturb up to twice if delay_reject is true. 
    accept_info(ii).non_acceptk = non_acceptk; % What perturbation number were we at when we perturbed the model? 0 will mean perturbed once, 1 will mean perturbed twice, 2 will mean perturbed once
    [model1,ptbnorm,ifpass,p_bd,Pm_prior1,...
    ptb,modptb,nchain,breakTrue,non_acceptk]...
    = delay_reject(...
        model, Pm_prior, ptb, ii, par, temp, Kbase,nchain,...
        model1,ptbnorm,p_bd,Pm_prior1k,non_acceptk); 
%     fprintf('\nptbnorm=%1.3f, non_acceptk=%1.0f\n',ptbnorm, non_acceptk) % Temporary

% % %     %%% Test stuff
% % %     for itest = [1:10000]; 
% % %         fprintf('\n         nac %1.0f', non_acceptk)
% % %         [model1,ptbnorm,ifpass,p_bd,Pm_prior1,...
% % %             ptb,modptb,nchain,breakTrue,non_acceptk]...
% % %             = delay_reject(...
% % %                 model, Pm_prior, ptb, ii, par, temp, Kbase,nchain,...
% % %                 model1,ptbnorm,p_bd,Pm_prior1k,non_acceptk); 
% % %         fprintf('\n%4.0f nac=%1.0f breaktrue=%1.0f',itest,non_acceptk,breakTrue)
% % %     end
% % %     %%%
    
    if breakTrue; break; end; % Using break here instead of break is equivalent to saying that if (1) this model has zero prior, then we want to bake that into our inverted PDFs in some way. Using a continue means we would just try ignoring the zero prior model and who knows what consequences that would have. Then (2), if ifaccept==false, then we simply already know that it has zero prior(?) and there is no point in conducting the forward modelling. brb2022.05.16.  

%% ===========================  FORWARD MODEL  ===========================
	% don't re-calc if the only thing perturbed is the error, or if there is zero probability of acceptance!
    if ~strcmp('sig',ptb{ii}(1:3)) || isempty(predata) % brb2022.03.06 If not recalculating, we get errors later. 
        % make random run ID (to avoid overwrites in parfor)
		ID = [chainstr,num2str(ii,'%05.f'),num2str(randi(99),'_%02.f')];
        
        try
            [predata,laymodel1] = b3__INIT_PREDATA(model1,par,trudata,0 );
% % %             if ii > 5; 
% % %                 error('Fake error, testing saving fail_chain > 0 model')
% % %             end
            [predata,par] = b3_FORWARD_MODEL_BW(model1,laymodel1,par,predata,ID,0,predataPrev);
            predata = b3_FORWARD_MODEL_RF_ccp(   model1,laymodel1,par,predata,ID,0 );
            predata = b3_FORWARD_MODEL_SW_kernel(model1,Kbase,par,predata );
            predataPrev = predata; % Keep track of last predata. To keep the previous complete HK stack. 
        catch e
            fail_chain=fail_chain+1;
            fprintf('\nbrb Forward model error on ii=%1.0f, fail_chain %.0f Report below:\n',ii,fail_chain)  
            fprintf('\n\n%s\n\n',getReport(e))
            break; % Break instead of continue -- equivalent to saying that if forward modelling doesn't work, then that model is crap and deserved a zero prior. We arent sure the best approach here yet. However, the forward modelling should usually work, so this shouldn't matter much. brb2022.05.16. 
        end

        % continue if any Sp or PS inhomogeneous or nan or weird output
        if ifforwardfail(predata,par)
            fail_chain=fail_chain+1;
            fprintf('Forward model error, run_one_chain line 230 or so. Fail_chain %.0f\n',fail_chain)  
            break; % Break instead of continue -- equivalent to saying that if forward modelling doesn't work, then that model is crap and deserved a zero prior. We arent sure the best approach here yet. However, the forward modelling should usually work, so this shouldn't matter much. brb2022.05.16. 
        end

        % process predata - filter/taper etc.
        predata0=predata; % save orig.
        for idt = 1:length(par.inv.datatypes)
            predata = predat_process( predata,par.inv.datatypes{idt},par);
        end

        %%% brb2022.04.05 Big change here? 
        %%% Was executing this often, but then model would be rejected. So
        %%% the kernels would never actually get reset. 
        %%% The change: Let's not use SW_precise. 
        %%% Just set newK = true. Then the next accepted model will involve
        %%% a kernel reset (not that only the next accepted model would
        %%% ever get a kernel reset anyway, but now we ensure that this
        %%% happens and we don't waste time on sw_precise for failed
        %%% models.
		% Explicitly use mineos + Tanimoto scripts if ptb is too large
% % %         if ptbnorm > par.inv.kerneltolmed; % todo not sure about this value. 
% % %             % Could also use par.inv.kerneltolmax somehow...
% % %             newK=true;
% % %             fprintf('\n ii=%1.0f Resetting kernel on next accepted model. NOT running non_acceptk',ii) 
% % %         end
        if ptbnorm/par.inv.kerneltolmax > random('unif',par.inv.kerneltolmed/par.inv.kerneltolmax,1,1) % control chance of going to MINEOS
            newK = true; 
            fprintf('\nsw presice: ptbnorm=%1.2f, ii=%1.0f, nonAcK=%1.0f\n',...
                ptbnorm, ii, non_acceptk)
            [ predata,SW_precise ] = b3_FORWARD_MODEL_SW_precise( model1,par,predata,ID );
        end

    end % only redo data if model has changed

%      plot_TRUvsPRE(trudata,predata);

    % continue if any Sp or PS inhomogeneous or nan or weird output
    if ifforwardfail(predata,par)
        fail_chain=fail_chain+1; ifpass=0;
        fprintf('Forward model error, fail_chain %.0f, ii=%1.0f\n',fail_chain,ii)
        break;
    else
        fail_chain = 0;
    end
       

%% =========================  CALCULATE MISFIT  ===========================

    % SW weights, if applicable
    [ SWwt ] = make_SW_weight( par,Kbase,trudata );
    [ misfit1 ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors

%% =======================  CALCULATE LIKELIHOOD  =========================
    [ log_likelihood1,misfit1 ] = b5_CALC_LIKELIHOOD( misfit1,trudata,model1.datahparm,par);
%     fprintf('MISFITS: Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.SpRF,misfit.PsRF,misfit.SW)
%     fprintf('CHI2S:   Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.chi2_sp,misfit.chi2_ps,misfit.chi2_SW)

    fail_chain = 0;
    predat_save1 = predata0;

    end % while ifpass
       
%% ========================  ACCEPTANCE CRITERION  ========================
    [ ifaccept ] = b6_IFACCEPT( log_likelihood1,log_likelihood,temp,p_bd*ifpass,Pm_prior1,Pm_prior);
    
    
%%% Figure out if chain is stuck
    stuckIfNoChange = 20; % p of not changing 13 times in row is about 0.0001. if there is 50% chance of acceptance each time and nothing has gone wrong.  
    if (ii>stuckIfNoChange+1) 
        % Temporary coding. accept_info.ifaccept doesn't exist until after
        % first iteration. In either case, don't check if stuck until we
        % try at least a few iterations. 
        accptArr = [accept_info.ifaccept]; % In case this wasn't already an array - it shouldn't be necessary after initializing ifaccept as an array brb2022.05.17
        if all( ~accptArr(ii-stuckIfNoChange:ii-1) ); % last stuckIfNoChange have been rejected. Somethings wrong! 
% % %             ifaccept = true; % Just accept this model. Hopefully will get us away from the model region that breaks the code.  Temporary solution. 
% % %             newK=true; % We might have got stuck because kernels needed to be reset. Let's reset them. 
% % %             fprintf('\nbrb2022.04.06 Stuck! Didnt accept model for %1.0f iterations. Simply accepted a new model to get unstuck. ii=%1.0f\n',stuckIfNoChange,ii)
            fprintf(['\nbrb2022.04.06 Stuck! Didnt accept model for %1.0f iterations.\n  ',...
                     'What errors came before this?. ii=%1.0f\n  ',...
                     'log like prev and current: %f3.2f, -> %f3.2f\n',...
                     'fail_chain=%1.0f'],...
                     stuckIfNoChange,ii,log_likelihood,log_likelihood1,fail_chain)
            fprintf('Increasing nchain artificially IF fail_chain == 0\nnchain: %1.0f',nchain); 
            if fail_chain == 0; 
                nchain = nchain + 50; 
            end
            fprintf(' -> nchain now: %1.0f\n', nchain); 

% % %             logLikArr = [accept_info(:).log_likelihood]; 
% % %             iterLik = [ii-stuckIfNoChange:ii-2];
% % %             logLikArr=logLikArr(iterLik);
% % %             figure(1); clf; hold on; plot(iterLik,logLikArr); 
% % %             ylabel('Log lik'); xlabel('Iter'); 
        end
    end
    %!%!
%%% End figure out if chain is stuck. 
    
% % %     if (ifaccept && non_acceptk == 2); 
% % %         disp('Accepted the second perturbation') 
% % %     end
    
    if delay_reject_bool; % brb2022.04.12 
        if (~ifaccept); 
            if (non_acceptk == 1); % Didn't accept a model, but we only perturbed it once. Without increasing ii, perturb the again and try iterating again. 
                ii = ii - 1; % Step back an iteration. After "continue", ii = ii + 1, so this just keeps us on the same iteration. 
                continue; 
            elseif (non_acceptk == 2); 
                ; 
            end
        end
    end
    
    dLog = (log_likelihood1-log_likelihood); 
    if ((log_likelihood1 - log_likelihood) < -10) && ifaccept; 
        disp('Accepted a shitty model')
    end
    
    % ======== PLOT ========  if accept
    if ifaccept && par.inv.verbose && fail_chain==0
        plot_TRUvsPRE( trudata,predata);  pause(0.001);
        if strcmp(projname,'SYNTHETICS')
            plot_MOD_TRUEvsTRIAL( TRUEmodel, model1 ); pause(0.001);
        end
    end
    
%% ========== TEMPORARY ====== make some plots and get insight into why acceptance does/doesn't happen
    accept_info(ii).ifaccept = ifaccept; 
    accept_info(ii).ifpass = ifpass; 
    accept_info(ii).misfit = misfit1; 
    accept_info(ii).log_likelihood = log_likelihood1; 
    accept_info(ii).model = model1; 
    accept_info(ii).predat_save = predat_save1; 
    accept_info(ii).iter = ii; 
    accept_info(ii).temp = temp; 
    accept_info(ii).ptbnorm = ptbnorm; 
    accept_info(ii).Pm_prior1 = Pm_prior1; 
    accept_info(ii).p_bd = p_bd; 
    if any('HKstack_P'==string(par.inv.datatypes)); 
        accept_info(ii).hk_Emax_per_iter = max(max(predata.HKstack_P.Esum)); 
        accept_info(ii).sig_hk = model1.datahparm.sig_HKstack_P; 
    else; 
        accept_info(ii).hk_Emax_per_iter = nan; 
        accept_info(ii).sig_hk = nan ;
    end
        
    if ii == 1; accept_info(1) .trudata = trudata; end; % Inefficient, but this is just for testing. 
    %%% TODO remove these lines eventually. 
    
    
%% Plot some inversion progress stuff. 
    try
        plot_progress_chain(absTimeIter, par, accept_info, ptb); 
    catch e 
        fprintf('\nProblem with plot_progress_chain. \n%s\n',getReport(e))
    end

%% ========================  IF ACCEPT ==> CHANGE TO NEW MODEL  =========================
    if ifaccept
        if par.inv.verbose
            fprintf('  *********************\n  Accepting model! logL:  %.4e ==>  %.4e\n  *********************\n',...
                log_likelihood,log_likelihood1)
            fprintf('                   Pm_+prior:  %.4e ==>  %.4e\n  *********************\n',...
                Pm_prior,Pm_prior1)
        end

        % save new model!
%         sprintf('ACCEPTED, non_acceptk = %1.0f', non_acceptk)
        non_acceptk = 0; %!%! We accepted a new model. So we have now rejected the current model 0 times. 
        model = model1;
        log_likelihood = log_likelihood1;
        Pm_prior = Pm_prior1;
        misfit = misfit1;
        predat_save = predat_save1;



    %% UPDATE KERNEL if needed
        if newK==true
            fprintf('\nModel changed ptbnorm=%1.3f since last kernel set. Running !=sw_precise and ==kernel_reset at ii=%1.0f.\n',ptbnorm,ii) 
%             [ predata,SW_precise ] = b3_FORWARD_MODEL_SW_precise( model1,par,predata,ID ); % brb2022.04.06 Now need to do this here since I'm no longer doing it earlier? Check github to verify. 
            KbasePrev = Kbase; 
            [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,0,SW_precise);
            % need to also reset likelihood and misfit to the new, precise data (likelihood may have been artificially high due to kernel forward  calc. approximation - if so, need to undo this, or chain will get stuck once we reset kernels).
% % %             [log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,model.datahparm);
% % %             Pm_prior = calc_Pm_prior(model,par); % This block executes only if accepted model. So no need to recalculate prior or likelihood. 
            nchain = 0;
        end

    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
        if newK, delete_mineos_files(ID,'R'); end
        if newK, delete_mineos_files(ID,'L'); end
    end

    % restart-chain if immediate failure
    if isinf(log_likelihood), fail_chain=100; break; end

%% =========  reset kernel at end of burn in or after too many iter =======
% Very important: kernel reset might be already done. This addresses the following conditions: 
% - We did not accept a model, but burning or maxnkchain requires us to reset kernel. We will reset USING THE LAST ACCEPTED MODEL. 
% - We did accept a model, but newk = false, and maxnkchains or burnin requires us to reset things. 
% % %     resetK = false;   
% % %     if (chainIsStuck && (fail_chain == 0)); % Chain gets stuck sometimes and it could be due to surface wave kernels having wrong likelihood: then we want to reset kernels very soon. But if we are stuck just because of forward models failing, there isnt a point in resetting the kernels. 
% % %         fprintf(' \n Increasing nchain for %s: chain was stuck.\n',chainstr); 
% % %         nchain = nchain + 50; 
% % %     end
    mightReset = (newK && ifaccept) == false; % In all cases EXCEPT when we just made K for an accepted model, then there is a chance we need to reset kernels. 
    if mightReset && ii == par.inv.burnin % reset kernel at end of burn in (and we didn't just reset it tacitly)
            fprintf('\n RECALCULATING %s KERNEL - end of burn in\n',chainstr);
            resetK = true;
    end
    if mightReset && nchain > par.inv.maxnkchain % reset kernel if chain too long (and we didn't just reset it tacitly)
            fprintf('\n RECALCULATING %s KERNEL at iter %.0f - chain too long\n',chainstr,ii);
            resetK = true;
    end
    
    par.hkResetInfo.timesSWKernelsReset = par.hkResetInfo.timesSWKernelsReset ...
    + int16(resetK); % If we are resetting surface wave kernels, also remake HK stacks with our current models parameters... but maybe only every N times kernel resets... 
    
    if resetK
        try
            % reset the kernels using the current model
            KbasePrev = Kbase; 
            [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,1);
            fail_reset = 0;
        catch e 
            fprintf('\nbrb Kernel reset fail: %s\n',getReport(e)); 
            % rewind back to last Kbase model that worked!
            fprintf('Kernel reset BROKE at %s-%.0f... REWIND to %s-%.0f <<<<<<<<<<<<<<< \n',chainstr,ii,chainstr,Kbase.itersave)
            model = Kbase.modelk;
            ii = Kbase.itersave;

            warning('brb2022.04.12. Replaced b3_FORWARD_MODEL_BW arguments with the newer argument order. Not sure if it works yet or not... Does the code ever even get to this point?')
%             predata = b3_FORWARD_MODEL_BW( model,par,trudata,ID,0 );
            predata = b3_FORWARD_MODEL_BW( model,laymodel,par,predata,ID,0,predataPrev); % brb2022.04.12 The arguments to forward_model_bw were in the wrong order. Probaly an old version of the code. 
            predata = b3_FORWARD_MODEL_RF_ccp( model,laymodel,par,predata,ID,0 );
            predata = b3_FORWARD_MODEL_SW_kernel( model,Kbase,par,predata );
            predataPrev = predata; % Keep track of last predata. To keep the previous complete HK stack. 
            
            fail_reset = fail_reset+1;
        end
        % need to also reset likelihood and misfit to the new, precise data (likelihood may have been artificially high due to kernel forward  calc. approximation - if so, need to undo this, or chain will get stuck once we reset kernels).
        [log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,model.datahparm);
        Pm_prior = calc_Pm_prior(model,par);
        nchain = 0;
    end

    
catch e % e is an MException struct
    fail_chain = fail_chain+1;
    fprintf('\nbrb Main master par catch: %s\n',getReport(e)); 
    fprintf('\nfail_chain = %1.0f',fail_chain); 
end % on try-catch

if newK||resetK, delete_mineos_files(ID,'R'); end
if newK||resetK, delete_mineos_files(ID,'L'); end

end % on iterations
%% -------------------------- End iteration  ------------------------------
end % on the fail_chain while...
% ----------
fprintf('\n ================= ENDING ITERATIONS %s =================\n',chainstr)
% save([resdir,'/chainout/',chainstr],'model0','misfits','allmodels')

plot_hk_progress_bool = true; if ~plot_hk_progress_bool; warning('Purposfully not plotting HK stack inversion progress'); end 
if plot_hk_progress_bool 
    try; 
        plot_h_kappa_progress2(trudata, allmodels, par.res.resdir, iii, accept_info, ...
            par, trudata.HKstack_P.Esum)
    catch e 
        fprintf('\nError in plot_h_kappa_progress2. This is normal if you dont have HK data. \n  %s  \n',getReport(e)); 
    end
end


save_inv_state(par.res.resdir,chainstr,allmodels,misfits)
misfits_perchain{iii} = misfits;
allmodels_perchain{iii} = allmodels;
preSW = nan; 
if isfield(trudata,'SW_Ray')
    SWs_perchain{iii} = preSW;
end

% save('accept_info', 'accept_info.mat'); 

diary off

% Aggregate files from this chain's folder to the stations folder. ram_to_HD later takes things from stations folder (RAM) to hard disk. 
[~,~]=clean_ram(ii,par,'clearRamInterval',1) % Prevent large amount of files from building up in ram, which can dramatically slow down matlab. 
[mvStatus, mvMessage]=system('mv -v ./* ..'); 
fprintf('\nMoving files from chain folder (%s)\n  to network folder (..).  \n  Result:\n  %s\n',pwd,mvMessage)
if ~mvStatus==0; warning('Possible problem moving files from chain folder to station folder. brb2022.03.30'); end


end