function [ifpass, numFails, model0, model,...
    Pm_prior, par, Kbase, model0_perchain] = initiate_model(...
        parOrig, TD, chainstr, fail_chain, model0_perchain, iii); 
    
par = parOrig; 


ifpass = 0;
numFails = 0; 
% only let starting model out of the loop if it passes conditions
while ifpass==0
    while ifpass==0 % first make sure the starting model satisfies conditions
        model0 = b1_INITIATE_MODEL(par,[],[],TD.Value);
        model = model0;
        model = expand_dathparms_to_data( model,TD.Value,par );
        ifpass = a1_TEST_CONDITIONS( model, par );
        Pm_prior = calc_Pm_prior(model,par);
    end

    %% starting model kernel
    fprintf('\nCreating starting kernels %s, have tried %1.0f times, fail_chain=%1.0f\n',...
        chainstr, numFails, fail_chain)
    try 
        [Kbase] = make_allkernels(model,[],TD.Value,['start',chainstr],par,...
            'maxrunN',5); % TODO this code caused a loop that would have failed infinitely. SHould not use a while loop like this. bb2021.09.14 
        %bb2021.09.23 Adding argument 5 where if run_mineos fails 5 (or some other small number of)times on the starting model, just toss it. No need to frett with getting a model to work with Mineos if no time is yet invested into inverting that model. 
    catch e 
        disp('Error in while - try,catch loop in MASTER_par.m')
        fprintf(1,'The identifier was:\n%s',e.identifier);
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        ifpass = false; 
        numFails = numFails + 1; 
        continue;
    end

    model0_perchain{iii} = model0;

end % now we have a starting model!
% diary off % Ending diary early for debugging. 
% figure(22); clf; hold on
% contourf(trudata.HKstack_P.K,trudata.HKstack_P.H,trudata.HKstack_P.Esum',20,'linestyle','none')
% set(gca,'ydir','reverse')
% figure(23); clf; hold on
% set(gca,'ydir','reverse')

end