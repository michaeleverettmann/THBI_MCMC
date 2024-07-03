function [mode] = f_mineos_run_phv_grv(parm, swperiods)

%% Wrap Josh's Matlab wrapper for getting dispursion curves from Mineos. 
parameter_FRECHET_save(parm, swperiods);  
parameter_FRECHET
a1_run_mineos_check
a3_pull_dispersion

end