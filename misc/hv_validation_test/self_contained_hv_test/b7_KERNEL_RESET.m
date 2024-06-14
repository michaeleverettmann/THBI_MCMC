function [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,iter,par,redo_phV,SW_precise)
% [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,iter,par,redo_phV,HVK_new)
%   Function to reset the surface wave kernels using the current model.
%   This function is called when the accepted model has diverged
%   sufficiently from the Kbase model that the kernel-computed
%   approximations to the data may be unreliable, and the full calculation
%   of data needs to be done, as well as updating the Kbase model and
%   computation of new kernels from this base. 
% 
%   In the case that the precise data were already calculated from this
%   model this iteration, the computationally intensive step of calculating
%   new data (and frechet files etc.) using MINEOS can be skipped, so we
%   set redo_phV to FALSE and import the HV kernels from above, where they
%   were calculated. In this case, only the MINEOS kernels need be computed
%   in this function, and the MINEOS files should then be deleted here. 
% 
%   In the case that this is a recomputation of the new data and kernels
%   from scratch, the full forward model has to be done here. In this case
%   redo_phV=TRUE, and the MINEOS (and HV-Tanimoto) calculations are redone
%   via the b3_FORWARD_MODEL_SW_precise function (the predicted data
%   are inserted directly into predata, whence they are procured for
%   insertion into Kbase).

if isempty(Kbase)
    Kbase = initiate_Kbase;
end

ran_ray_vel = false; % Have we ran this data type? Only do calculation for it once. 
ran_lov_vel = false; 
ran_hv = false; 

if nargin < 6 || isempty(redo_phV) || redo_phV==true
    % recalculate the data (phase/group velocities and HVratios)
    % using the precise model, if this has not been done already this
    % iteration
    [ predata,SW_precise ] = b3_FORWARD_MODEL_SW_precise( model,par,predata,ID );
end

Kbase.modelk = model;
Kbase.itersave = iter;
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; 
    pdtyp = parse_dtype(dtype); 

    % check the data type is a surface wave, otherwise skip - not doing BW
    % here
    if ~strcmp(pdtyp{1},'SW'), continue; end
    
    if strcmp(pdtyp{2},'HV')
        if ran_hv; continue; end; 
        HVr_new = predata.SW_HV.HVr;
        Kbase = populate_Kbase( Kbase, dtype, HVr_new, [], {SW_precise.HV.HVK_new},...
            predata.(dtype).for_mod_info.periods_calc ); 
        
        if ~ran_hv; ran_hv = true; end; 
    else
        if ran_ray_vel && strcmp(pdtyp{2},'Ray'); continue; end; 
        if ran_lov_vel && strcmp(pdtyp{2},'Lov'); continue; end; 
        
        MINEOS_file_delete = 1;
        ifplot = 0;
        
        grV = SW_precise.(pdtyp{2}).grV; 
        phV = SW_precise.(pdtyp{2}).phV; 
        
        par_mineos = struct('R_or_L',pdtyp{2},'phV_or_grV',pdtyp{3},'ID',ID);
        
        K = run_kernels(predata.(dtype).for_mod_info.periods_calc, ...
            par_mineos,SW_precise.(pdtyp{2}).eigfiles, MINEOS_file_delete, ifplot, par.inv.verbose);
        Kbase = populate_Kbase( Kbase, dtype, phV, grV, {K}, ...
            predata.(dtype).for_mod_info.periods_calc );
        
        if ~ran_ray_vel && strcmp(pdtyp{2},'Ray'); ran_ray_vel = true; end; 
        if ~ran_lov_vel && strcmp(pdtyp{2},'Lov'); ran_lov_vel = true; end; 
        
    end
    
end



end

