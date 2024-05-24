function [Kbase] = make_allkernels(model,Kbase,data,ID,par,options)
    arguments
        model
        Kbase
        data
        ID
        par
        options.maxrunN
    end
% [Kbase] = make_allkernels(model,Kbase,periods,ID,par)
%   Function to make all surface wave kernels

if isempty(Kbase)
    Kbase = initiate_Kbase;
end

ran_ray_vel = false; % Have we ran this data type? Only do calculation for it once. 
ran_lov_vel = false; 
ran_hv = false; 

Kbase.modelk = model;
Kbase.itersave=0;
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id}; 
    pdtyp = parse_dtype(dtype); 

    if ~strcmp(pdtyp{1},'SW'), continue; end
    
    if strcmp(pdtyp{2},'HV')
        if ran_hv; continue; end; 
        
        for_mod_info = data.(dtype).for_mod_info;
        periods_calc = for_mod_info.periods_calc; % This "for_mod_info" is supposed to be consistent between all versions of Rayleigh wave grv/phv, regardless of author/dataset. This way, we only do kernel forward modelling for one author, but we address all periods at once. I designed it that way in load_data.m. brb2022.06.24
        [HVr_new,HVK_new] = run_HVkernel(model,periods_calc,ID,1,0,par.inv.verbose);
        Kbase = populate_Kbase( Kbase,dtype,HVr_new,[],{HVK_new}, periods_calc );    
        
        if ~ran_hv; ran_hv = true; end; 
        
    else
        if ran_ray_vel && strcmp(pdtyp{2},'Ray'); continue; end; 
        if ran_lov_vel && strcmp(pdtyp{2},'Lov'); continue; end; 
        
        par_mineos = struct('R_or_L',pdtyp{2},'phV_or_grV',pdtyp{3},'ID',ID);
        for_mod_info = data.(dtype).for_mod_info; 
        periods_calc = for_mod_info.periods_calc; % This "for_mod_info" is supposed to be consistent between all versions of Rayleigh wave grv/phv, regardless of author/dataset. This way, we only do kernel forward modelling for one author, but we address all periods at once. I designed it that way in load_data.m. brb2022.06.24
        [phV,grV,eigfiles] = run_mineos(model,periods_calc,par_mineos,0,0,par.inv.verbose,options.maxrunN); %MINEOS_REPLACE
        K = run_kernels(periods_calc,par_mineos,eigfiles,1,0,par.inv.verbose); %MINEOS_REPLACE
        Kbase = populate_Kbase( Kbase,dtype,phV,grV,{K}, periods_calc );
        
        if ~ran_ray_vel && strcmp(pdtyp{2},'Ray'); ran_ray_vel = true; end; 
        if ~ran_lov_vel && strcmp(pdtyp{2},'Lov'); ran_lov_vel = true; end; 

    end
    
end

end

