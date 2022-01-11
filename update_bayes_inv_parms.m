function [par, inv] = update_bayes_inv_parms(parOrig, stamp);
% Run default bayes_inv_parms, then update with this function. 
% There should be a unique stamp for each unique test. 
% The stamp should determine what folder the results are saved in. 

par = parOrig; 
str_temp = split(stamp, '/'); 

if     strcmp(stamp, 'ENAM_trial'); 
    disp('Using default parameters') 
    

% % % 'h_kappa_tests/degree_of_freedom/*'
% % elseif strcmp(str_temp{1}, 'h_kappa_tests'               ) && ...
% %        strcmp(str_temp{2}, 'degree_of_freedom'           ); % Set stamp to this to run this test
% %     
% %     disp('Testing degrees of freedom with h-kappa stacks')
% %     things = split(stamp, '/'); 
% %     value = str2num(things{end}); % Degrees of freedom is last part of stamp
% %     
% %     par.mod.data.deg_of_freedom.h_kappa = value; 
% %     
% %     
% % % 'h_kappa_tests/h_kappa_solo/sigma_constant/*'    
% % elseif strcmp(str_temp{1}, 'h_kappa_tests'   ) && ...
% %        strcmp(str_temp{2}, 'h_kappa_solo'    ) && ...
% %        strcmp(str_temp{3}, 'sigma_constant'  );     
% %    
% %     disp('h_kappa test. Solo. Constant sigma. ')
% %     things = split(stamp, '/'); 
% %     value = str2num(things{end}); % Value is last part of stamp
% % 
% %     par.mod.data.logstd_sigma = 0.0; 
% %     par.mod.data.prior_sigma.HKstack = struct('P',value);
% %     par.inv.datatypes = {'HKstack_P'}; 
% % 
% % 
% % 
% % 
% % elseif strcmp(str_temp{1}, 'h_kappa_tests'                         ) && ...
% %        strcmp(str_temp{2}, 'observe_convergence'                   ) && ...
% %        strcmp(str_temp{3}, 'invert_sigma'                          );     
% %    
% % %     disp('h_kappa test. Solo. Constant sigma. ')
% %     things = split(stamp, '/'); 
% %     value = str2num(things{end}); % Value is last part of stamp
% % 
% %     par.mod.data.prior_sigma.HKstack = struct('P',value); 
% %     
% %     par.mod.data.logstd_sigma = 0.05; 
% %     par.inv.datatypes = {'HKstack_P'}; 
% %     par.inv.niter = 12000; 
% %     par.inv.burnin = 4000; 
% %     par.inv.cooloff = 3000; 
% %     par.inv.nchains = 16; 
% % elseif strcmp(str_temp{1}, 'h_kappa_tests'                         ) && ...
% %        strcmp(str_temp{2}, 'observe_convergence'                   ) && ...
% %        strcmp(str_temp{3}, 'invert_sigma_false'                          );     
% %    
% % %     disp('h_kappa test. Solo. Constant sigma. ')
% %     things = split(stamp, '/'); 
% %     value = str2num(things{end}); % Value is last part of stamp
% % 
% %     par.mod.data.prior_sigma.HKstack = struct('P',value);     
% %     
% %     par.mod.data.logstd_sigma = 0.00; 
% %     par.mod.data.prior_sigma.HKstack = struct('P',value); 
% %     par.inv.datatypes = {'HKstack_P'};     
% %     par.inv.niter = 12000; 
% %     par.inv.burnin = 4000; 
% %     par.inv.cooloff = 3000; 
% %     par.inv.nchains = 16; 
elseif strcmp(str_temp{1}, 'h_kappa_tests'                         ) && ...
       strcmp(str_temp{2}, 'observe_convergence'                   ) && ...
       strcmp(str_temp{3}, 'invert_sigma'                          ) && ...
       strcmp(str_temp{4}, 'large_cooldown'                        );     
   
%     disp('h_kappa test. Solo. Constant sigma. ')
    things = split(stamp, '/'); 
    value = str2num(things{end}); % Value is last part of stamp

    par.mod.data.prior_sigma.HKstack = struct('P',value);     
    
    par.mod.data.logstd_sigma = 0.05; 
    par.mod.data.prior_sigma.HKstack = struct('P',value); 
    par.inv.datatypes = {'HKstack_P'};     
    par.inv.niter = 16000; 
    par.inv.burnin = 8000; 
    par.inv.cooloff = 6000;
    par.inv.tempmax = 9; 
    par.inv.nchains = 16; 
    
    
    
    
    
    
    
    
    
    
end



inv = par.inv; 

end


