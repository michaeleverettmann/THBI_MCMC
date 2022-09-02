function [model, ii, predata, predataPrev, laymodel, ...
    log_likelihood, misfit, Pm_prior, nchain] = ...
    reverse_inversion_to_last_kbase(...
    Kbase, par, trudata, predata, ID, predataPrev);  
    
model = Kbase.modelk; 
ii = Kbase.itersave; 
[predata,laymodel] = b3__INIT_PREDATA(model,par,trudata,0 );
predata = b3_FORWARD_MODEL_BW( model,laymodel,par,predata,ID,0,predataPrev); % brb2022.04.12 The arguments to forward_model_bw were in the wrong order. Probaly an old version of the code. 
predata = b3_FORWARD_MODEL_RF_ccp( model,laymodel,par,predata,ID,0 );
predata = b3_FORWARD_MODEL_SW_kernel( model,Kbase,par,predata );
predataPrev = predata; % Keep track of last predata. To keep the previous complete HK stack.            

for idt = 1:length(par.inv.datatypes) % In case length of receiver functions changes. brb2022.08.11. 
    predata = predat_process( predata,par.inv.datatypes{idt},par);
end

[log_likelihood,misfit] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,model.datahparm);
Pm_prior = calc_Pm_prior(model,par);
nchain = 0;

end