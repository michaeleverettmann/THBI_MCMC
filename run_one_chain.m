% par = PR.Value; 
% trudata = TD.Value

% run_one_chain(par, nwk, sta, iii)



function [model0, misfits, allmodels, preSW ...
    ]=run_one_chain(par, trudata, nwk, sta, iii ); 

% run_one_chain(par, trudata, nwk, sta, iii,model0_perchain, misfits_perchain, allmodels_perchain, SWs_perchain )
% % % 
% % % % For some reason matlab doesn't like returning a variable that was an input. 
% % % model0_perchain    = model0_perchainO    ;
% % % misfits_perchain   = misfits_perchainO   ;
% % % allmodels_perchain = allmodels_perchainO ;
% % % SWs_perchain       = SWs_perchainO       ;

paths = getPaths(); 

% Disable a bspline warning that doesn't seem to matter. Needs to be placed in parfor or else individual workers don't keep this warning off. ; 
warning('off', 'MATLAB:rankDeficientMatrix'); % This comes up when doing least squares inversion for spline weights. Be careful, the rankDeficientMatrix could be needed at another point in the inversion...    

absTimeIter = timeseries(zeros(par.inv.niter,1)); 
absTimeIter(1) = datetime(); 

accept_info = struct('ifaccept'          , zeros(par.inv.niter,1)*nan,...; ; % Initialize structure with acceptance info for each iteration. 
                     'misfit'            , zeros(par.inv.niter,1)*nan,...
                     'log_likelihood'    , zeros(par.inv.niter,1)*nan,...
                     'sig_hk'            , zeros(par.inv.niter,1)*nan,...
                     'model'             , zeros(par.inv.niter,1)*nan,...
                     'predat_save'       , zeros(par.inv.niter,1)*nan,...
                     'iter'              , zeros(par.inv.niter,1)*nan,...
                     'temp'              , zeros(par.inv.niter,1)*nan,...
                     'ptbnorm'           , zeros(par.inv.niter,1)*nan,...
                     'Pm_prior1'         , zeros(par.inv.niter,1)*nan,...
                     'p_bd'              , zeros(par.inv.niter,1)*nan,...
                     'non_acceptk'       , zeros(par.inv.niter,1)*nan,...
                     'hk_Emax_per_iter'  , zeros(par.inv.niter,1)*nan,...
                     'trudata'           , zeros(par.inv.niter,1)*nan); 

for iInit=1:par.inv.niter;
    accept_info(iInit).ifaccept          = []; ; % Initialize structure with acceptance info for each iteration. 
    accept_info(iInit).misfit            = []; % misfit1; 
    accept_info(iInit).log_likelihood    = []; % log_likelihood1; 
    accept_info(iInit).sig_hk            = []; % model1.datahparm.sig_HKstack_P; 
    accept_info(iInit).model             = []; % model1; 
    accept_info(iInit).predat_save       = []; % predat_save1; 
    accept_info(iInit).iter              = []; % ii; 
    accept_info(iInit).temp              = []; % temp; 
    accept_info(iInit).ptbnorm           = []; % ptbnorm; 
    accept_info(iInit).Pm_prior1         = []; % Pm_prior1; 
    accept_info(iInit).p_bd              = []; % p_bd; 
    accept_info(iInit).non_acceptk       = []; % non_acceptk; 
    accept_info(iInit).hk_Emax_per_iter  = []; % max(max(predata.HKstack_P.Esum)); 
    accept_info(iInit).trudata           = []; % trudata; end; % Inefficient, but this is just for testing. 
end


chainstr = [nwk '.' sta '_' mkchainstr(iii)];
par.res.chainstr = chainstr; 
par.res.chainID = mkchainstr(iii); % ID to keep track of things like where the velocity profiles are. brb2022.03.30
% mkdir(par.res.chainID); 
par.res.chainExecFold = [paths.ramDrive '/' nwk '_' sta '/' par.res.chainID]; % Explicitly determine executable folder, since we will use rm ./* and don't want to risk executing this from somewhere on the hard drive!!! brb2022.03.30
mkdir(par.res.chainExecFold); 
cd(par.res.chainExecFold);  % Speed up execution of code by reducing number of files in a chains folder (this should help in Matlab system() calls). brb2022.03.30

%%% Diary file stuff
diaryFile = sprintf('diary_%s.txt', chainstr); % Use a diary file to keep track of parallel inversions seperately. %TODO_STATION_NETWORK bb2021.11.12
diary off; 
delete(diaryFile); 
pause(0.01); % Because delete seems to fully execute after turning diary on...
diary(diaryFile);
%%% End diary file stuff. 

%% Fail-safe to restart chain if there's a succession of failures
fail_chain=20;
while fail_chain>=20

%%% Initialize several variables. 
newK = false; % Maybe make true. 
SW_precise = []; 
laymodel1 = []; 
non_acceptk = 0; % How often have we rejected the current baseline model. 
predataPrev = []; 
% accept_info = struct(); % bb2022.01.05 Temporary tests for h-kappa inversion. 
timeStartIter = zeros(par.inv.niter,1); % Vector of starting times of each iteration. 
timeVeryStart = now(); % Starting time of inversion
tFact = (24 * 60 * 60); % Factor to convert from serial time to seconds. 
%%% End initializing several varibales. 


%% Prep posterior structure
[ misfits,allmodels,savedat ] = b0_RESULTS_SETUP(par);

%% initiate model
[ifpass, numFails, model0, model,...
    Pm_prior, par, Kbase] = initiate_model(...
        par, trudata, chainstr, fail_chain, iii); 
model0_perchain{iii} = model0;


%% ========================================================================
%% ------------------------- Start iterations -----------------------------
%% ========================================================================
ptb = cell({});
nchain = 0;
fail_chain = 0; fail_reset = 0;
ifaccept=true;
if isfield(trudata,'SW_Ray')
    preSW = zeros(length(trudata.SW_Ray.periods),ceil(par.inv.niter./par.inv.saveperN));
end
% reset_likelihood;
log_likelihood = -Inf;
predata=[]; predat_save = []; misfit = [];
% not parfor

fprintf('\n =========== STARTING ITERATIONS %s ===========\n',chainstr)
ii = 0;
time0 = now;

par.hkResetInfo = struct('timesReset',0, 'timesSWKernelsReset',0); % Structure to keep track of when to reset HK stacks. 
while ii < par.inv.niter
ii = ii+1; par.ii = ii; % bb2022.10.18 Keep track of ii in par for easily plotting how inversion is changing with ii. 

% Maintenance in case folder was abused during last iterations. 
[~,~]=clean_ram(ii,par,'clearRamInterval',200) % Prevent large amount of files from building up in ram, which can dramatically slow down matlab. If the ram is getting very full, then there is probably a problem!
if ~strcmp(pwd,par.res.chainExecFold); 
    warning('Somehow you left the ram folder! This can tremendously slow down code execute. Changing back (note that cd is slow)... Fix it. brb2022.04.02'); 
    cd(par.res.chainExecFold); 
end


% Quickly check to make sure inversion isn't going way to slow. 
% absTimeIter(ii) = now() * tFact; 

timeStartIter(ii) = (now() - timeVeryStart) * tFact;  % Seconds since the inversion very first started. 
littleIter=max(1,ii-100); smallTime = timeStartIter(littleIter); bigTime = timeStartIter(ii); 
avgTime = (bigTime-smallTime) / (ii-littleIter); % Average computational time over past min([1,100]) iterations. 
if ii ~= 1; absTimeIter(ii) = datetime(); end; % Dont replace first time. Want that in case fail_chain is a PITA
% Finish getting computational time. 


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
    
%     %%% Kill chain if too slow. 
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
%     %%% End kill chain if too slow. 
   
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
    
    delay_reject_bool = false; if (mod(ii, 100)==0) && ~delay_reject_bool; warning('brb2022.04.12 non_acceptk=0. NO DELAYED REJECTION.'); end; % brb2022.04.12. 
    if ~ delay_reject_bool; 
        non_acceptk = 0; 
    end
    
    [model1,ptbnorm,ifpass,p_bd,Pm_prior1,...
    ptb,modptb,nchain,breakTrue,non_acceptk]...
    = delay_reject(model, Pm_prior, ptb, ii, par, temp, Kbase,nchain,...
        model1,ptbnorm,p_bd,Pm_prior1k,non_acceptk); 
%     fprintf('\nptbnorm=%1.3f, non_acceptk=%1.0f\n',ptbnorm, non_acceptk) % Temporary
    
    if breakTrue; break; end;

%% ===========================  FORWARD MODEL  ===========================
	% don't re-calc if the only thing perturbed is the error, or if there
	% is zero probability of acceptance!
    if ~strcmp('sig',ptb{ii}(1:3)) || isempty(predata) % brb2022.03.06 If not recalculating, we get errors later. 
        % make random run ID (to avoid overwrites in parfor)
		ID = [chainstr,num2str(ii,'%05.f'),num2str(randi(99),'_%02.f')];

        try
            [predata,laymodel1] = b3__INIT_PREDATA(model1,par,trudata,0 );
            [predata,par] = b3_FORWARD_MODEL_BW(model1,laymodel1,par,predata,ID,0,predataPrev);
            predata = b3_FORWARD_MODEL_RF_ccp(   model1,laymodel1,par,predata,ID,0 );
            predata = b3_FORWARD_MODEL_SW_kernel(model1,Kbase,par,predata );
            predataPrev = predata; % Keep track of last predata. To keep the previous complete HK stack. 
        catch e
            fail_chain=fail_chain+1;
            fprintf('\nbrb Forward model error, fail_chain %.0f Report below:\n',fail_chain)  
            fprintf('\n\n%s\n\n',getReport(e))
            break;
        end

        % continue if any Sp or PS inhomogeneous or nan or weird output
        if ifforwardfail(predata,par)
            fail_chain=fail_chain+1;
            fprintf('Forward model error, fail_chain %.0f\n',fail_chain)  
            break;
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
        warning('This forward fail forces ifpass = 0. Im not sure, maybe this can make your chain get stuck. brb2022.04.12.')
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
    stuckIfNoChange = 30; % p of not changing 13 times in row is about 0.0001. if there is 50% chance of acceptance each time and nothing has gone wrong.  
    if (ii>stuckIfNoChange+1); 
        % Temporary coding. accept_info.ifaccept doesn't exist until after
        % first iteration. In either case, don't check if stuck until we
        % try at least a few iterations. 
        accptArr = [accept_info.ifaccept]; 
        if all( ~accptArr(ii-stuckIfNoChange:ii-1) ); % last stuckIfNoChange have been rejected. Somethings wrong! 
% % %             ifaccept = true; % Just accept this model. Hopefully will get us away from the model region that breaks the code.  Temporary solution. 
% % %             newK=true; % We might have got stuck because kernels needed to be reset. Let's reset them. 
% % %             fprintf('\nbrb2022.04.06 Stuck! Didnt accept model for %1.0f iterations. Simply accepted a new model to get unstuck. ii=%1.0f\n',stuckIfNoChange,ii)
            fprintf(['\nbrb2022.04.06 Stuck! Didnt accept model for %1.0f iterations.\n  ',...
                     'What errors came before this?. ii=%1.0f\n  ',...
                     'log like prev and current: %f3.2f, -> %f3.2f'],...
                     stuckIfNoChange,ii,log_likelihood,log_likelihood1)
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
    accept_info(ii).misfit = misfit1; 
    accept_info(ii).log_likelihood = log_likelihood1; 
    accept_info(ii).model = model1; 
    accept_info(ii).predat_save = predat_save1; 
    accept_info(ii).iter = ii; 
    accept_info(ii).temp = temp; 
    accept_info(ii).ptbnorm = ptbnorm; 
    accept_info(ii).Pm_prior1 = Pm_prior1; 
    accept_info(ii).p_bd = p_bd; 
    accept_info(ii).non_acceptk = non_acceptk; 
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
            [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,0,SW_precise);
            nchain = 0;
            % need to also reset likelihood and misfit to the new, precise data (likelihood may have been artificially high due to kernel forward  calc. approximation - if so, need to undo this, or chain will get stuck once we reset kernels).
            [log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,model.datahparm);
            Pm_prior = calc_Pm_prior(model,par);
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
    resetK = false;
    mightReset = (newK == false) && (ifaccept == true); 
    if mightReset && ii == par.inv.burnin % reset kernel at end of burn in (and we didn't just reset it tacitly)
            % bb2022.10.12 TODO Need to find way to run kernel on first accepted model after burnin
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
%     fprintf(1,'\n\nLine %1.0f\nMain master par loop error. Error message:\n%s\n\n',...
%         e.stack.line,e.message);
%     if par.inv.verbose, fprintf('  --SOME ERROR--\n'); end
    fprintf('\nbrb Main master par catch: %s\n',getReport(e)); 
    fail_chain = fail_chain+1;
end % on try-catch

if newK||resetK, delete_mineos_files(ID,'R'); end
if newK||resetK, delete_mineos_files(ID,'L'); end

end % on iterations
%% -------------------------- End iteration  ------------------------------
end % on the fail_chain while...
% ----------
fprintf('\n ================= ENDING ITERATIONS %s =================\n',chainstr)
% save([resdir,'/chainout/',chainstr],'model0','misfits','allmodels')

try; 
    plot_h_kappa_progress2(trudata, allmodels, par.res.resdir, iii, accept_info, ...
        par, trudata.HKstack_P.Esum)
catch e 
    fprintf('\nError in plot_h_kappa_progress2. \n  %s  \n',getReport(e)); 
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