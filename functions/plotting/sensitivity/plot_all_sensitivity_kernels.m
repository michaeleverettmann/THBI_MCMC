function plot_all_sensitivity_kernels(final_model,trudata,par,Kbase,resdir)
model = final_model; 
model.VS = model.VSav; 
model.VP = model.VPav; 
model.rho = model.rhoav; 
model.Panis = model.Panisav; 
model.Sanis = model.Sanisav; 
model.z = model.Z; 
[Kbase,~] = b7_KERNEL_RESET(model,[],trudata,'junk',0,par,1); 
plot_sensitivity_kernel_hv(Kbase.HV.KHV,'model',model,'predata',final_predata,...
    'filename',sprintf('%s/sensitivity_hv.pdf',resdir) ); 
plot_sensitivity_kernel_ray(Kbase.Ray.Kph,'model',model,'predata',final_predata,...
    'filename',sprintf('%s/sensitivity_ray.pdf',resdir) ); 
plot_sensitivity_kernel_love(Kbase.Lov.Kph,'model',model,'predata',final_predata,...
    'filename',sprintf('%s/sensitivity_love.pdf',resdir) ); 
end