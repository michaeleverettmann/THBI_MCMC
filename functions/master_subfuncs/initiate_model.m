function [ifpass, numFails, model0, model,...
    Pm_prior, par, Kbase] = initiate_model(...
        parOrig, tdvalue, chainstr, fail_chain, iii); 
% Make the inversions starting model. 
    
par = parOrig; 

ifpass = 0;
numFails = 0; 
% only let starting model out of the loop if it passes conditions
while ifpass==0
    while ifpass==0 % first make sure the starting model satisfies conditions

        model0 = b1_INITIATE_MODEL_simple(par,[],[],tdvalue); % Use the simple option here, to start with same reasonable model each time. 
        % figure(31), clf, hold on, set(gca,'ydir','reverse')
        % plot(model0.VS,model0.z,'b',model0.VP,model0.z,'r',model0.rho,model0.z,'k')

        model_start = model0; 
        disp('Making init model')
        for ipert = 1:(40*model0.M); % On average, perturb each parameter some number of times. 
            model0 = b2_PERTURB_MODEL(model0, par, 0.125); % Small temperature, so perturbations are small. 
        end 

        model = model0;
        model = expand_dathparms_to_data( model,tdvalue,par );
        ifpass = a1_TEST_CONDITIONS( model, par );
        Pm_prior = calc_Pm_prior(model,par);
        
        % plot(model0.VS,model0.z,'b.',model0.VP,model0.z,'r.',model0.rho,model0.z,'k.')
        % [model.crustmparm.xi; model.mantmparm.xi]; 

        % If we have a bunch of knots, toss it. Too complicated for a starting model. 
        if (model_start.crustmparm.Nsp+1 < model.crustmparm.Nsp) || ...
            (model_start.mantmparm.Nsp+1 < model.mantmparm.Nsp); 
            ifpass = false; 
        end 
    end

    %% starting model kernel
    fprintf('\nCreating starting kernels %s, have tried %1.0f times, fail_chain=%1.0f\n',...
        chainstr, numFails, fail_chain)
    try 
        [Kbase] = make_allkernels(model,[],tdvalue,['start',chainstr],par,...
            'maxrunN',5); % TODO this code caused a loop that would have failed infinitely. SHould not use a while loop like this. bb2021.09.14 
        %bb2021.09.23 Adding argument 5 where if run_mineos fails 5 (or some other small number of)times on the starting model, just toss it. No need to frett with getting a model to work with Mineos if no time is yet invested into inverting that model. 
    catch e 
        fprintf('\n\n%s\n\n',getReport(e))
        ifpass = false; 
        numFails = numFails + 1; 
        continue;
    end

end % now we have a starting model!

% figure(22); clf; hold on
% contourf(trudata.HKstack_P.K,trudata.HKstack_P.H,trudata.HKstack_P.Esum',20,'linestyle','none')
% set(gca,'ydir','reverse')
% figure(23); clf; hold on
% set(gca,'ydir','reverse')

%%% See the model. 
% figure(1); clf; hold on; 
% tiledlayout(1,3); 
% 
% nexttile; hold on; box on; set(gca, 'YDir', 'reverse'); 
% title('Rho'); 
% plot(model.rho, model.z); 
% 
% nexttile; hold on; box on; set(gca, 'YDir', 'reverse'); 
% title('Vs'); 
% plot(model.VS, model.z); 
% 
% nexttile; hold on; box on; set(gca, 'YDir', 'reverse'); 
% title('Vp'); 
% plot(model.VP, model.z); 

end