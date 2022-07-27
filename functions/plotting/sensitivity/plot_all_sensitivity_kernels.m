function [Kbase] = plot_all_sensitivity_kernels(final_model,trudata,par,resdir,options)
    arguments
        final_model
        trudata
        par
        resdir
        options.Kbase = []; 
    end
    
model = final_model; 
model.VS = model.VSav; 
model.VP = model.VPav; 
model.rho = model.rhoav; 
model.Panis = model.Panisav; 
model.Sanis = model.Sanisav; 
model.z = model.Z; 

if isempty(options.Kbase); % Get Kbase if not provided. brb2022.07.15 - I don't think this is made any more in run_all_chains. 
    [Kbase,~] = b7_KERNEL_RESET(model,[],trudata,'junk',0,par,1); 
else
    Kbase = options.Kbase; 
end

if any( contains(string(par.inv.datatypes), 'SW_HV'    ))
    plot_sensitivity_kernel_hv(Kbase.HV.KHV,'dat',Kbase.HV,'model',model,...
        'filename',sprintf('%s/sensitivity_hv.pdf',resdir) ); 
end

if any( contains(string(par.inv.datatypes), 'SW_Ray_phV'))
    plot_sensitivity_kernel_ray(Kbase.Ray.Kph,'dat',Kbase.Ray,'model',model,...
        'filename',sprintf('%s/sensitivity_ray.pdf',resdir) ); 
end

if any( contains(string(par.inv.datatypes), 'SW_Lov_phV'))
    plot_sensitivity_kernel_love(Kbase.Lov.Kph,'dat',Kbase.Lov,'model',model,...
        'filename',sprintf('%s/sensitivity_love.pdf',resdir) ); 
end


end