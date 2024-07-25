function [model0, misfits, allmodels, preSW ...
    ]=run_one_chain(par, trudata, nwk, sta, iii ); 
% Function to execute one chain. Is called in a loop over all chains
% from run_all_chains.m. This function gets really complicated because
% Fortran codes (Mineos...) can break at so many different points in the
% inversion, requiring different tactics for keeping a chain going. 

%% Some setup
paths = getPaths(); 

% Disable a bspline warning that doesn't seem to matter. This comes up when doing least squares inversion for spline weights.     
warning('off', 'MATLAB:rankDeficientMatrix'); 

if strcmp(par.synth.propmat_or_telewavesim, 'telewavesim'); 
    pyenv('Version', '~/opt/anaconda3/envs/tws/bin/python', ... % Use anaconda environment where telewavesim is installed. This affects the entire Matlab session. TODO define this path in somewhere more obvious. 
        'ExecutionMode','OutOfProcess'); % ERROR ALERT Could not import numpy if using an anconda environment. Matlab would simply crash. However, setting executionMode=OutOfProcess fixed that for me. https://www.mathworks.com/matlabcentral/answers/502458-why-can-py-numpy-array-not-be-resolved
end

% Time the inversion across iterations. 
absTimeIter = timeseries(zeros(par.inv.niter,1)); 
absTimeIter(1) = datetime(); 

% Folder/name for this chain. 
chainstr = [nwk '.' sta '_' mkchainstr(iii)];
par.res.chainstr = chainstr; 
par.res.chainID = mkchainstr(iii); % ID to keep track of things like where the velocity profiles are. brb2022.03.30
par.res.chainExecFold = [paths.ramDrive '/' nwk '_' sta '/' par.res.chainID]; % Explicitly determine executable folder, since we will use rm ./* and don't want to risk executing this from somewhere on the hard drive!!! brb2022.03.30
mkdir(par.res.chainExecFold); 
cd(par.res.chainExecFold);  % Keep each chain in its own folder. 

% Structure to keep track of accepting/rejecting models and the related conditions
accept_info = make_accept_info(par); 

% Diary file for keeping the command output for each chain. 
diaryFile = sprintf('diary_%s.txt', chainstr); 
diary off; 
delete(diaryFile); 
pause(0.01); % Because delete seems to fully execute after turning diary on...
diary(diaryFile);

%% Fail-safe to restart chain if there's a succession of failures
fail_chain=20;
fail_total = 0; 
iter_fail = []; % Keep track of all iterations where there was an error. Needs to be resizable array, and keep track of previous failures even after resetting a chain. 

%% Start this chain. Try again if it fails. 
while fail_chain>=20

% Initialize several variables. 
newK = false; % Have a new kernel for use? TODO maybe make true. 
SW_precise = []; 
KbasePrev = []; % When this is empty, we've only made kernels once. Can reset to this model if needed. 
laymodel1 = []; 
non_acceptk = 0; % How often have we rejected the current baseline model. 
predataPrev = []; 
timeStartIter = zeros(par.inv.niter,1); % Vector of starting times of each iteration. 
timeVeryStart = now(); % Starting time of inversion
tFact = (24 * 60 * 60); % Factor to convert from serial time to seconds. 

%% Prep posterior structure
[ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par);

%% initiate model
[ifpass, numFails, model0, model,...
    Pm_prior, par, Kbase] = ...
    initiate_model(...
        par, trudata, chainstr, fail_chain, iii); 
[~,~]=clean_ram(1,par,'clearRamInterval',1,'verbose',true); 


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

ii_no_reverse = 0; 
while ii < par.inv.niter
fprintf('ii=%1.0f',ii) % Don't add newlines. This way, we just get one line that shows each iteration. 

earlyFailNum = 100; 
save_model_if_error = false; % Debug! Make this true, or something like (ii > earlyFailNum) if you want to get the models that are breaking the code.  
if save_model_if_error && (fail_chain > 0); 
    if exist('e', 'var'); thiserr = e; else; thiserr = '??????'; end
    save(sprintf('%s/model-caused-error_%s_ii-%1.0f_time-%s.mat',...
        par.res.resdir,par.res.chainstr,ii,string(datestr(now,'mm-dd-yyyy_HH-MM-SS'))),...
        'model', 'model1', 'trudata', 'par', 'thiserr', 'Kbase', 'KbasePrev'); % Should be everything needed to reproduce a forward model error. 
end
           
ii = ii+1; par.ii = ii; % bb2022.10.18 Keep track of ii in par for easily plotting how inversion is changing with ii. 
ii_no_reverse = ii_no_reverse + 1; % This will keep going up, even if we rewind ii. 
accept_info(ii).fail_chain = fail_chain; 
accept_info(ii).fail_total = fail_total; 

% There is obnoxiously complicated logic, basically due to the fact that
% many of the Fortran codes often simply break. We need to be able to
% rewind to models that don't break the Fortran code. 
% To prevent one chain from taking way too long, there is also logic
% to reduce the chances that we proceed with any terrible chains that
% frequently, but don't always, break the Fortran codes. 
if fail_chain > 0; % Problem -- There are various ways of dealing with this.
    % fprintf('\nFail chain: %2.0f. ii = %6.0f',fail_chain, ii)
    % TODO add in Kbase.numRewinds. If it gets to > 3 or something, go to prev Kbase.
    resetData = false; 
    failInfoStr = sprintf(['',...
        '%s ii=%1.0f Failure. fail_chain = %1.0f, fail_total previously = %1.0f. '],...
            chainstr, ii, fail_chain, fail_total); 
        
    % If early fail
    if ii < earlyFailNum; % If chain is giving problems this early, just reset it. Such a chain will probably give many errors later in the inversion, even if it stumbles across a model that doesn't break Fortran codes now. 
        warning('%sError within first %1.0f iterations. Reseting chain.\n', failInfoStr, earlyFailNum);  
        fail_chain = 100; 
        fail_total = fail_total + 1; break; 
    end 
    
    % If later fail
    if fail_chain > 15; 
        fprintf('\n :| :| :| :| :| :| :| :| :| :|   \n'); 
        if ii - Kbase.itersave > 50; % Last Kbase has to be a decent model - immediately after making Kbase (currently this is at lines 580 during the inversion) we checked for signs of forward model failure. Kbase is also made in initiate model around line 60. If that Kbase was broken, the inversion would reset at iteration 1. 
            warning([failInfoStr 'High fail chain. Reset to last Kbase.']); 
            KbaseUseReset = Kbase; 
            resetData = true; 
            fail_total = fail_total + 1; 
        else % Last Kbase might also be bad. Go back two models. 
            if ~isempty(KbasePrev); 
                warning([failInfoStr 'High fail chain. Reset to last KbasePrev (NOT Kbase).']); 
                KbaseUseReset = KbasePrev; 
                resetData = true; 
                fail_total = fail_total + 1; 
            else; % We don't know if last Kbase is good, and we only have one Kbase. Just reset the chain: we must not be very far anyway.  
                warning([failInfoStr 'High fail chain after reset Kbase, but KbasePrev doesnt exist. Resetting chain.']); 
                fail_chain = 100; 
                fail_total = fail_total + 1; break; 
            end
        end
    end
       
    if resetData; 
        try; 
            fprintf(['\n *********************************************** \n',...
                'Trying to reset data according to previous Kbase.\n'])
            
                [model, ii, predata, predataPrev, laymodel, ...
                    log_likelihood, misfit, Pm_prior, nchain] = ...
                    reverse_inversion_to_last_kbase(...
                        KbaseUseReset, par, trudata, predata, ID, predataPrev); 
                    
            fail_chain = 0; 
            if ifforwardfail(predata,par); 
                error('How is it possible to have forward fail now when we successfully did forward modelling with this model earlier? Ending inversion.'); 
            end
            
        catch e
            warning('\n-----------\nCHAIN BROKE. You should always be able to rewind to a model and run all forward modelling on it. A model wasnt accepted if we couldnt run the forward modelling. We simply have to exit the inversion. Report: \n%s\n--------------\n',getReport(e)); 
%             for asdfasdf=[1:30]; fprintf('\nREMOVE THIS SAVE LINE saving temporarily\n'); end % Uncomment this to save the workspace to help you debug. 
%                 save([par.res.resdir '/fail_rewind_chain_' num2str(iii) '.mat'])
%                 fail_chain = - 100; 
%             break;
        end
    end
    
    fprintf('\nNow on iteration %1.0f\n',ii); 
    iter_fail(end+1) = ii; 
    %!%! Code to break chain entirely if reset to many times.     
end

% Maintenance in case folder is getting filled up during last iterations. 
[~,~]=clean_ram(ii,par,'clearRamInterval',500) % Prevent large amount of files from building up in ram, which can dramatically slow down matlab. If the ram is getting very full, then there is probably a problem!
if ~strcmp(pwd,par.res.chainExecFold); 
    warning('\nSomehow you left the ram folder! This can tremendously slow down code execute. Changing back (note that cd is slow)... Fix it. brb2022.04.02'); 
    cd(par.res.chainExecFold); 
end

% Quickly check to make sure inversion isn't going way too slow. Take a moving average of inversion speed. 
timeStartIter(ii) = (now() - timeVeryStart) * tFact;  % Seconds since the inversion very first started. 
littleIter=max(1,ii-100); smallTime = timeStartIter(littleIter); bigTime = timeStartIter(ii); 
avgTime = (bigTime-smallTime) / (ii-littleIter); % Average computational time over past min([1,100]) iterations. 
if ii ~= 1; absTimeIter(ii) = datetime(); end; % Dont replace first time. Want that in case fail_chain is a PITA

%% SAVE model every saveperN
if mod(ii,par.inv.saveperN)==0 && log_likelihood ~= -Inf
    [misfits,allmodels,savedat] = b9_SAVE_RESULT(ii,log_likelihood,misfit,model,Pm_prior,misfits,allmodels,predat_save,savedat,time0);
end

%% SAVE inv state every Nsavestate iterations
if rem(ii,par.inv.Nsavestate)==0
    save_inv_state(par.res.resdir,chainstr,allmodels,misfits)
end

try  
    %% Figure out if chain is too slow or unstable. 
    % TODO changed from 4 to 1 times par.inv.savepern. So, fairly verbose. 
    if rem(ii,1*par.inv.saveperN)==0 || ii==1; 
            fprintf(...
            'Sta %s nwk %s  Iteration %s %1.0f, (%1.0f with rewinds), Average comp time (last 100 iter) %2.2f (all) %2.2f\n',...
            par.data.stadeets.sta,par.data.stadeets.nwk,chainstr,ii, ii_no_reverse,...
            avgTime, timeStartIter(ii)/ii); end
    if par.inv.verbose, pause(0.05); end % TODO not sure why this is needed. 
    newK = false; resetK = false;
    
    % Newer, simpler way of deciding when to kill chain. Rely on reversing chains. Only consider quitting early if we have done successfull+failed iterations > par.inv.niter. At that point, kill the chain if it fails at resetting X times IN A ROW
    if (ii_no_reverse > par.inv.niter) && (fail_reset>4) ; % If we did more iterations than we meant to, consider exiting the chain - Then if the kernels fail to reset X times in a row, call it quits. 
        fail_chain = -100; 
        break;
    end  
    if ii_no_reverse > (par.inv.niter * 22000/16000); % Due to reversing, most chains take a little less than 22k iterations to finish if I aimed for 16k. So just stop all chains at 22/16 of par.inv.niter. 
        fail_chain = -100; 
        break;
    end
    
% Below are some ways we previously determined to give up on a chain. They
% are left here as comments in case the current system is not suitable for
% a project. Maybe these can be a useful starting place. 
% % %     if fail_chain>19
% % %         fail_total = fail_total + 1;
% % %         % if not enough saved in this chain, abort and restart
% % %         if (ii - par.inv.burnin)/par.inv.saveperN < 200
% % %             fprintf('\nBreak chain -- spot a1\n') 
% % %             break
% % %         % if enough saved in chain, abort and keep the incomplete chain
% % %         else
% % %             fprintf('\nBreak chain -- spot a2\n')
% % %             fail_chain = -fail_chain; break
% % %         end
% % %     end
    
% % %     if fail_reset>5
% % %         fail_total = fail_total + 1; 
% % %         % if cannot reset kernels because current saved model is not viable
% % %         if (ii - par.inv.burnin)/par.inv.saveperN < 200
% % %             fprintf('\nBreak chain -- spot a3\n')
% % %             fail_chain = 100;  break % high fail_chain will mean we restart chain
% % %         % if enough saved in chain, abort and keep the incomplete chain
% % %         else
% % %             fprintf('\nBreak chain -- spot a4\n')
% % %             fail_chain = -100; break
% % %         end
% % %     end
    
% % %     %%% Kill chain if too slow. 
% % %     if avgTime > 1; % Things are running really slow... why? Might want to restart, unless we are at the beginning of iterations. 
% % % %         if ~feature('IsDebugMode'); % Don't break slow chains if we are debugging. 
% % %         tScaleKill = 1; 
% % %         tryKeeping = ((ii - par.inv.burnin)/par.inv.saveperN) >= 200; 
% % %         toBreak = ( (avgTime > 20  * tScaleKill ) && (ii > 3   ) ) || ...
% % %                   ( (avgTime > 7   * tScaleKill ) && (ii > 20  ) ) || ...
% % %                   ( (avgTime > 3.5 * tScaleKill ) && (ii > 200 ) ) || ...
% % %                   ( (avgTime > 2   * tScaleKill ) && (ii > 500 ) ) || ...
% % %                   ( (avgTime > 1.2 * tScaleKill ) && (ii > 1000) ); % Things should be going smoothly by now and always take ~ .6 s. However, there might be a time where mineos runs a few times and makes a good chain seem slow: maybe averaging over just 100 iterations could make us reject an otherwise ok chain. 
% % %         if toBreak && tryKeeping; 
% % %             fprintf('\nExecution too slow. Aborting (and keeping chain) %s ii=%1.0f, avg time last 100 iter = %1.2f\n',...
% % %                 chainstr,ii,avgTime);
% % %             fail_total = fail_total + 1;
% % %             fail_chain = -100; break; 
% % %         elseif toBreak; 
% % %             fprintf('\nExecution to slow. Resetting %s ii=%1.0f, avg time last 100 iter = %1.2f\n',...
% % %                 chainstr,ii,avgTime);
% % %             fail_total = fail_total + 1;
% % %             fail_chain =  100;  break; 
% % %         end  
% % % %         end
% % %     end
% % %     %%% End Figure out if chain is too slow or unstable. 
    %%
   
    % temperature - for perturbation scaling and likelihood increase
    temp = (par.inv.tempmax-1)*erfc(2*(ii-1)./par.inv.cooloff) + 1;

    ifpass = false;
    ifaccept=false;
    while ifpass == false % only keep calculating if model passes (otherwise save and move on)
        
        %% ===========================  PERTURB  ===========================
        if ii==1; log_likelihood1 = -Inf; model1 = model; ptbnorm = nan; ...
                Pm_prior1k = Pm_prior; end
        if non_acceptk == 0; p_bd = 1; end; 
        
        % Delayed rejection. If we didn't accept the last model due to decreased data fit, try perturbing it one more time. If that model doesn't work still, then reject the two models. 
        delay_reject_bool = true; 
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
        % fprintf('\nptbnorm=%1.3f, non_acceptk=%1.0f\n',ptbnorm, non_acceptk) % Temporary
        
        if breakTrue; break; end; % break while ifpass loop. Using break here instead of break is equivalent to saying that if (1) this model has zero prior, then we want to bake that into our inverted PDFs in some way. Using a continue means we would just try ignoring the zero prior model and who knows what consequences that would have. Then (2), if ifaccept==false, then we simply already know that it has zero prior(?) and there is no point in conducting the forward modelling. brb2022.05.16.  
    
        %% ===========================  FORWARD MODEL  ===========================
	    % don't re-calc if the only thing perturbed is the error, or if there is zero probability of acceptance!
        if ~strcmp('sig',ptb{ii}(1:3)) || isempty(predata) % brb2022.03.06 If not recalculating, we get errors later. 
            % make random run ID (to avoid overwrites in parfor)
		    ID = [chainstr,num2str(ii,'%05.f'),num2str(randi(99),'_%02.f')];
            
            try % Do forward modelling. 
                [predata,laymodel1] = b3__INIT_PREDATA(model1,par,trudata,0 );
                [predata,par] = b3_FORWARD_MODEL_BW(model1,laymodel1,par,predata,ID,0,predataPrev);
                predata = b3_FORWARD_MODEL_RF_ccp(   model1,laymodel1,par,predata,ID,0 );
                predata = b3_FORWARD_MODEL_SW_kernel(model1,Kbase,par,predata );
                predataPrev = predata; % Keep track of last predata. To keep the previous complete HK stack. 
            catch e % Forward modelling got broke. 
                fail_chain=fail_chain+1;
                fprintf('\nbrb Forward model error on ii=%1.0f, fail_chain %.0f Report below:\n',ii,fail_chain)  
                fprintf('\n\n%s\n\n',getReport(e))
                break; % Break instead of continue -- equivalent to saying that if forward modelling doesn't work, then that model is crap and deserved a zero prior. We arent sure the best approach here yet. However, the forward modelling should usually work after the burning, so this shouldn't matter much. brb2022.05.16. 
            end
    
            % break if any Sp or PS inhomogeneous or nan or weird output
            if ifforwardfail(predata,par)
                fail_chain=fail_chain+1;
                fprintf('Forward model error, run_one_chain line 329 or so. Fail_chain %.0f\n',fail_chain)  
                break; 
            end
    
            % process predata - filter/taper etc.
            predata0=predata; % save orig.
            for idt = 1:length(par.inv.datatypes)
                predata = predat_process( predata,par.inv.datatypes{idt},par);
            end
    
            newK = false; % brb2022.07.01. TODO might remove. I'm in the process of not using SW precise at this point in the inversion. 
        end % only redo data if model has changed
    
        % continue if any Sp or PS inhomogeneous or nan or weird output
        if ifforwardfail(predata,par)
            fail_chain=fail_chain+1; ifpass=0;
            fprintf('Forward model error, fail_chain %.0f, ii=%1.0f\n',fail_chain,ii)
            break;
        else
            fail_chain = 0; % Not sure this is necessary. brb2022.07.26  
        end   
    
        %% =========================  CALCULATE MISFIT  ===========================
        % SW weights, if applicable
        [ SWwt ] = make_SW_weight( par,Kbase,trudata );
        [ misfit1 ] = b4_CALC_MISFIT( trudata,predata,par,0,SWwt ); % misfit has structures of summed errors
    
        %% =======================  CALCULATE LIKELIHOOD  =========================
        [ log_likelihood1,misfit1 ] = b5_CALC_LIKELIHOOD( misfit1,trudata,model1.datahparm,par);
        % fprintf('MISFITS: Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.SpRF,misfit.PsRF,misfit.SW)
        % fprintf('CHI2S:   Sp %5.2e  Ps %5.2e  SW %5.2e\n',misfit.chi2_sp,misfit.chi2_ps,misfit.chi2_SW)
    
        if isinf(log_likelihood1) || isinf(-log_likelihood1); 
            fail_chain=fail_chain+1; ifpass=0;
            fprintf('\nFailure at ii=%1.0f. log_likelihood1=%f, chi2=%1.2f\n',ii, log_likelihood1, misfit1.chi2sum)
            break;
        end

        fail_chain = 0; % This chain succeeded, so we now have 0 failed chains in a row?  
        predat_save1 = predata0;
        
        %% ========================  ACCEPTANCE CRITERION  ========================
        [ ifaccept ] = b6_IFACCEPT( log_likelihood1,log_likelihood,temp,p_bd*ifpass,Pm_prior1,Pm_prior);
    
    end 

    %%% Figure out if chain is stuck
    stuckIfNoChange = 20; % Determine if stuck. For reference, probability of not changing 13 times in row is about 0.0001, assuming there is 50% chance of acceptance each time and nothing has gone wrong with codes/calculations.  
    if (ii>stuckIfNoChange+1) 
        accptArr = [accept_info.ifaccept]; % In case this wasn't already an array - it shouldn't be necessary after initializing ifaccept as an array brb2022.05.17
        if all( ~accptArr(ii-stuckIfNoChange:ii-1) ); % last stuckIfNoChange have been rejected. Somethings wrong! 
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
        end
    end
    
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
        disp('Accepted a shitty model. This should not happen very often. ')
    end
    
    % ======== PLOT ========  if accept
    if ifaccept && par.inv.verbose && fail_chain==0
        plot_TRUvsPRE( trudata,predata);  pause(0.001);
        % if strcmp(projname,'SYNTHETICS')
        %     plot_MOD_TRUEvsTRIAL( TRUEmodel, model1 ); pause(0.001);
        % end
    end
    
    %% Gather info to make some plots. Can help diagnose why there is never acceptance for some chain. 
    % TODO I think is is very inefficient. However, accept_info does get used in other places in the code. This might be a target for speeding the code up. Check with profiler.     % However, it might have just been inefficient if storing .model during each iteration, which eats ram. brb20240612
    accept_info(ii).ifaccept = ifaccept; 
    accept_info(ii).ifpass = ifpass; 
    accept_info(ii).misfit = misfit1; 
    accept_info(ii).log_likelihood = log_likelihood1; 
    % accept_info(ii).model = model1; warning('Saving all models') 
    accept_info(ii).iter = ii; 
    accept_info(ii).temp = temp; 
    accept_info(ii).ptbnorm = ptbnorm; 
    accept_info(ii).Pm_prior1 = Pm_prior1; 
    accept_info(ii).p_bd = p_bd; 
    if any('HKstack_P'==string(par.inv.datatypes)); 
        accept_info(ii).hk_Emax_per_iter = max(max(predata.HKstack_P.Esum)); 
        accept_info(ii).sig_hk = model1.datahparm.sig_HKstack_P; 
        accept_info(ii).predat_save = struct('HKstack_P', struct('H',predat_save1.HKstack_P.H,'K',predat_save1.HKstack_P.K, 'E_by_Emax',predat_save1.HKstack_P.E_by_Emax)); % Kind of a hack. Used to plot HK stack progress stuff in plot_h_kappa_progress2.m. Don't store the rest of predata thuogh, it takes WAY too much ram. 
    else; 
        accept_info(ii).hk_Emax_per_iter = nan; 
        accept_info(ii).sig_hk = nan;
        accept_info(ii).predat_save = nan; 
    end
        
    if ii == 1; accept_info(1).trudata = trudata; end;    
    
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
        % sprintf('ACCEPTED, non_acceptk = %1.0f', non_acceptk)
        non_acceptk = 0; %!%! We accepted a new model. So we have now rejected the current model 0 times. 
        model = model1;
        log_likelihood = log_likelihood1;
        Pm_prior = Pm_prior1;
        misfit = misfit1;
        predat_save = predat_save1;

        % Below was the old logic to update kernels. brb
        %% UPDATE KERNEL if needed
        % % % % %         if newK==true
        % % % % %             fprintf('\nModel changed ptbnorm=%1.3f since last kernel set. Running !=sw_precise and ==kernel_reset at ii=%1.0f.\n',ptbnorm,ii) 
        % % % % % %             [ predata,SW_precise ] = b3_FORWARD_MODEL_SW_precise( model1,par,predata,ID ); % brb2022.04.06 Now need to do this here since I'm no longer doing it earlier? Check github to verify. 
        % % % % %             KbasePrev = Kbase; 
        % % % % %             [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,0,SW_precise);
        % % % % %             % need to also reset likelihood and misfit to the new, precise data (likelihood may have been artificially high due to kernel forward  calc. approximation - if so, need to undo this, or chain will get stuck once we reset kernels).
        % % % % % % % %             [log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,model.datahparm);
        % % % % % % % %             Pm_prior = calc_Pm_prior(model,par); % This block executes only if accepted model. So no need to recalculate prior or likelihood. 
        % % % % %             nchain = 0;
        % % % % %         end

    else
        if par.inv.verbose, fprintf('  --FAIL--\n'); end
    end

    % restart-chain if immediate failure
    if isinf(log_likelihood), fail_chain=100; break; end

    %% =========  reset kernel at end of burn in or after too many iter =======
    mightReset = true; 
    if mightReset && ii == par.inv.burnin % reset kernel at end of burn in (and we didn't just reset it tacitly)
        fprintf('\n RECALCULATING %s KERNEL - end of burn in\n',chainstr);
        resetK = true;
    end
    if mightReset && nchain > par.inv.maxnkchain % reset kernel if chain too long (and we didn't just reset it tacitly)
        fprintf('\n RECALCULATING %s KERNEL at iter %.0f - chain too long\n',chainstr,ii);
        resetK = true;
    end
    if (ifaccept) && (ptbnorm>2); 
        fprintf('\n RECALCULATING %s KERNEL at iter %.0f - ptbnorm large = %1.3f\n',chainstr,ii,ptbnorm);
        resetK = true; 
    end
    
    % If we are resetting surface wave kernels, also remake HK stacks with our current models parameters. Resetting the HK stack is somewhat slow. 
    par.hkResetInfo.timesSWKernelsReset = par.hkResetInfo.timesSWKernelsReset + int16(resetK); 
    
    if resetK % Kernel recalculation
        try % Try making new kernels based on the last accepted model. 
            if par.inv.verbose; fprintf('\nTrying to reset K.\n'); end
            [Kbase, KbasePrev, predata, log_likelihood, misfit, Pm_prior, nchain, ifpass] = ...
                compare_kernel_full_calc(...
                    model, Kbase, predata, trudata, ID, par, fail_chain, ...
                    ii, ifpass, misfit, ptbnorm, log_likelihood, ifaccept); 
            fail_reset = 0; 
        catch e % If it doesn't work, rewind to last Kbase model that worked. 
            fprintf('\nbrb Kernel reset fail. Error was: \n --------\n %s\n --------\n',getReport(e)); 
            fprintf('Kernel reset BROKE at %s-%.0f... REWIND to %s-%.0f <<<<<<<<<<<<<<< \n',chainstr,ii,chainstr,Kbase.itersave)
            
            KbaseUseReset = Kbase; 
            [model, ii, predata, predataPrev, laymodel, ...
                log_likelihood, misfit, Pm_prior, nchain] = ...
                reverse_inversion_to_last_kbase(...
                    KbaseUseReset, par, trudata, predata, ID, predataPrev);             
            
            fail_reset = fail_reset+1;
            fail_chain = fail_chain+1;  
            fail_total = fail_total+1; 
            
            if ifforwardfail(predata,par); 
                warning('You should always be able to rewind to a model and run all forward modelling on it, but it just happened. A model wasnt accepted if we couldnt run the forward modelling. What is happening?. We simply have to exit the inversion.'); 
                fail_chain = - 100; 
                break;
            end
        end
    end
    
catch e % e is an MException struct
    fail_chain = fail_chain+1;
    fprintf('\nbrb Main master par catch: %s\n',getReport(e)); % try encapsulates the bulk of the chains code, so we have to report the errors manually. 
    fprintf('\nfail_chain = %1.0f',fail_chain); 
    fprintf('\n--------------------------------------------------\n\n\n'); 
end 

if resetK, delete_mineos_files(ID,'R'); end
if resetK, delete_mineos_files(ID,'L'); end

if resetK; % Remove Mineos files to keep ram available. 
    [~,~]=clean_ram(ii,par,'clearRamInterval',1,'verbose',true); 
end

end % on iterations
%% -------------------------- End iteration  ------------------------------
end % on the fail_chain while...
% ----------
fprintf('\n ================= ENDING ITERATIONS %s =================\n',chainstr)

% Plot to show the different HK values that were sampled. Helps debug. 
plot_hk_progress_bool = true; if ~plot_hk_progress_bool; warning('Purposfully not plotting HK stack inversion progress'); end 
if plot_hk_progress_bool && any(("HKstack_P" == string(par.inv.datatypes)))
    try 
        plot_h_kappa_progress2(trudata, allmodels, par.res.resdir, iii, accept_info, ...
            par, trudata.HKstack_P.Esum)
    catch e 
        fprintf('\nError in plot_h_kappa_progress2: \n ------ \n %s \n ------- \n',getReport(e)); 
    end
end

save_inv_state(par.res.resdir,chainstr,allmodels,misfits)
misfits_perchain{iii} = misfits;
allmodels_perchain{iii} = allmodels;
preSW = nan; 
if isfield(trudata,'SW_Ray')
    SWs_perchain{iii} = preSW;
end

diary off

% Aggregate files from this chain's folder to the stations folder. ram_to_HD later takes things from stations folder (RAM) to hard disk. 
[~,~]=clean_ram(ii,par,'clearRamInterval',1) % Prevent large amount of files from building up in ram, which can dramatically slow down matlab. 
[mvStatus, mvMessage]=system('mv -v ./* ..'); 
fprintf('\nMoving files from chain folder (%s)\n  to network folder (..).  \n  Result:\n  %s\n',pwd,mvMessage)
if ~mvStatus==0; warning('Possible problem moving files from chain folder to station folder. brb2022.03.30'); end

end