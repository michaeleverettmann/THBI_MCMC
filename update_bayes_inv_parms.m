function [par, inv] = update_bayes_inv_parms(parOrig, stamp);
% Run default bayes_inv_parms, then update with this function. 
% There should be a unique stamp for each unique test. 
% The stamp should determine what folder the results are saved in. 

par = parOrig; 
str_temp = split(stamp, '/'); 


% { Different data types alone
% niter_only   = 12000; 
% burnin_only  = 4000; 
% cooloff_only = 3000; 
% nchains_only = 16; 
niter_only   = 300; 
burnin_only  = 100; 
cooloff_only = 80 ; 
nchains_only = 1  ; 
% niter_only   = 600; 
% burnin_only  = 100; 
% cooloff_only = 80 ; 
% nchains_only = 2  ; 

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
% % elseif strcmp(str_temp{1}, 'h_kappa_tests'                         ) && ...
% %        strcmp(str_temp{2}, 'observe_convergence'                   ) && ...
% %        strcmp(str_temp{3}, 'invert_sigma'                          ) && ...
% %        strcmp(str_temp{4}, 'large_cooldown'                        );     
% %    
% % %     disp('h_kappa test. Solo. Constant sigma. ')
% %     things = split(stamp, '/'); 
% %     value = str2num(things{end}); % Value is last part of stamp
% % 
% %     par.mod.data.prior_sigma.HKstack = struct('P',value);     
% %     
% %     par.mod.data.logstd_sigma = 0.05; 
% %     par.mod.data.prior_sigma.HKstack = struct('P',value); 
% %     par.inv.datatypes = {'HKstack_P'};     
% %     par.inv.niter = 16000; 
% %     par.inv.burnin = 8000; 
% %     par.inv.cooloff = 6000;
% %     par.inv.tempmax = 9; 
% %     par.inv.nchains = 16; 
% %     
    
elseif strcmp(stamp, 'lower_prior_sigma')
%     disp('Do stuff')
    par.mod.data.prior_sigma.HKstack.P = 0.1; 
        
    
elseif strcmp(stamp, 'prior_sigma_0.2_min_sigma_0.1')
%     disp('Do stuff')
    par.mod.data.prior_sigma.HKstack.P = 0.2;   
    par.mod.data.min_sigma.HKstack.P   = 0.1;     

    
elseif strcmp(stamp, 'prior_sigma_1_min_sigma_0.5')
%     disp('Do stuff')
    par.mod.data.prior_sigma.HKstack.P = 1;   
    par.mod.data.min_sigma.HKstack.P   = 0.5;    
    
    
elseif strcmp(stamp, 'prior_sigma_1_min_sigma_0.1_err')
%     disp('Do stuff')
    par.mod.data.prior_sigma.HKstack.P = 1;   
    par.mod.data.min_sigma.HKstack.P   = 0.1;    
    
    
elseif strcmp(stamp, 'prior_sigma_10_min_sigma_3')
%     disp('Do stuff')
    par.mod.data.prior_sigma.HKstack.P = 10;   
    par.mod.data.min_sigma.HKstack.P   = 3;    
    
    
elseif strcmp(stamp, 'HKappa_001'); 
    par.mod.data.prior_sigma.HKstack.P = .6; 
    par.mod.data.min_sigma.HKstack.P   = .1; 
    par.datprocess.HKappa.min_error    = 0.1; 
    par.datprocess.HKappa.scale_error  = 1; 
elseif strcmp(stamp, 'HKappa_002'); 
    par.mod.data.prior_sigma.HKstack.P = 1; 
    par.mod.data.min_sigma.HKstack.P   = .2; 
    par.datprocess.HKappa.min_error    = 0.1; 
    par.datprocess.HKappa.scale_error  = 1; 
elseif strcmp(stamp, 'HKappa_003'); 
    par.mod.data.prior_sigma.HKstack.P = 1; 
    par.mod.data.min_sigma.HKstack.P   = .01; 
    par.datprocess.HKappa.min_error    = 0.1; 
    par.datprocess.HKappa.scale_error  = 1; 
elseif strcmp(stamp, 'HKappa_004'); 
    par.mod.data.prior_sigma.HKstack.P = 1; 
    par.mod.data.min_sigma.HKstack.P   = .5; 
    par.datprocess.HKappa.min_error    = 0.01; 
    par.datprocess.HKappa.scale_error  = 100; 
    
elseif strcmp(stamp, 'HKappa_005'); 
    par.datprocess.HKappa.weightDistanceMax = 1; 
    par.datprocess.HKappa.weightDistanceMin = 1; 
elseif strcmp(stamp, 'HKappa_006'); 
    par.datprocess.HKappa.weightDistanceMax = 1; 
    par.datprocess.HKappa.weightDistanceMin = 0.5;
    
elseif strcmp(stamp, 'HKappa_007'); 
    par.mod.starting.HKappa.startAtHK = true; 
    
    
    
elseif strcmp(stamp, 'all_001'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 16000; 
    par.inv.burnin                    = 6000; 
    par.inv.cooloff                       = 4500; 
    par.inv.nchains = 16; 
    
    
% { --- Different data types alone
elseif strcmp(stamp, 'SW_Ray_phV_only'); 
    par.inv.datatypes                 = {'SW_Ray_phV'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only; 
elseif strcmp(stamp, 'SW_Lov_phV_only'); 
    par.inv.datatypes                 = {'SW_Lov_phV'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only; 
elseif strcmp(stamp, 'RF_Sp_ccp_only'); 
    par.inv.datatypes                 = {'RF_Sp_ccp'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only; 
elseif strcmp(stamp, 'HKstack_P_only'); 
    par.inv.datatypes                 = {'HKstack_P'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only; 
elseif strcmp(stamp, 'SW_HV_only'); 
    par.inv.datatypes                 = {'SW_HV'}   ;
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only;     
% } --- End test different data types alone   


elseif strcmp(stamp, 'all_demo'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 2000; 
    par.inv.burnin                    = 700 ; 
    par.inv.cooloff                   = 500 ; 
    par.inv.nchains                   = 16  ; 
elseif strcmp(stamp, 'all_002'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 16000; 
    par.inv.burnin                    = 5000 ; 
    par.inv.cooloff                   = 4000 ; 
    par.inv.nchains                   = 16   ; 
elseif strcmp(stamp, 'all_no_hv'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 16000; 
    par.inv.burnin                    = 5000 ; 
    par.inv.cooloff                   = 4000 ; 
    par.inv.nchains                   = 16   ; 
    
    
end



  


inv = par.inv; 

end


