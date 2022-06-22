function [par, inv] = update_bayes_inv_parms(parOrig, stamp);
% Run default bayes_inv_parms, then update with this function. 
% There should be a unique stamp for each unique test. 
% The stamp should determine what folder the results are saved in. 

par = parOrig; 
str_temp = split(stamp, '/'); 


% { Different data types alone
% niter_only   = 16000; 
% burnin_only  = 4000; 
% cooloff_only = 3000; 
% nchains_only = 12; 
% 
niter_only   = 300; 
burnin_only  = 100; 
cooloff_only =  70; 
nchains_only = 2; warning('Short "only" stamp');  

% niter_only   = 4000; 
% burnin_only  = 1000; 
% cooloff_only = 800 ; 
% nchains_only = 16  ; 
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
% { --- ONE CHAIN Different data types alone
elseif strcmp(stamp, 'SW_Ray_phV_only_one_chain'); 
    par.inv.datatypes                 = {'SW_Ray_phV'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = 1; 
elseif strcmp(stamp, 'SW_Lov_phV_only_one_chain'); 
    par.inv.datatypes                 = {'SW_Lov_phV'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = 1; 
elseif strcmp(stamp, 'RF_Sp_ccp_only_one_chain'); 
    par.inv.datatypes                 = {'RF_Sp_ccp'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = 1; 
elseif strcmp(stamp, 'HKstack_P_only_one_chain'); 
    par.inv.datatypes                 = {'HKstack_P'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = 1; 
elseif strcmp(stamp, 'SW_HV_only_one_chain'); 
    par.inv.datatypes                 = {'SW_HV'}   ;
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = 1;     
% } --- End test different data types alone   

% { --- A couple data types. 
elseif strcmp(stamp, 'HKstack_P___SW_HV_only__only'); 
    par.inv.datatypes                 = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only;  
% } ---


elseif strcmp(stamp, 'all_demo'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 500; 
    par.inv.burnin                    = 100 ; 
    par.inv.cooloff                   = 80 ; 
    par.inv.nchains                   = 1  ; 
elseif strcmp(stamp, 'all_002'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = niter_only; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only ; 
    par.inv.nchains                   = nchains_only   ; 
elseif strcmp(stamp, 'all_no_hv'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P'}; 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 16000; 
    par.inv.burnin                    = 5000 ; 
    par.inv.cooloff                   = 4000 ; 
    par.inv.nchains                   = 16   ; 
    
    
elseif strcmp(stamp, 'HK_fast'); 
    par.inv.datatypes                 = {'HKstack_P'};
    par.inv.niter                     = 500  ; 
    par.inv.burnin                    = 100  ; 
    par.inv.cooloff                   = 80   ; 
    par.inv.nchains                   = 16    ; 
elseif strcmp(stamp, 'SW_HV_fast'); 
    par.inv.datatypes                 = {'SW_HV'}   ;
    par.inv.niter                     = 500  ; 
    par.inv.burnin                    = 100  ; 
    par.inv.cooloff                   = 80   ; 
    par.inv.nchains                   = 16    ;  
elseif strcmp(stamp, 'HK_faster'); 
    par.inv.datatypes                 = {'HKstack_P'};
    par.inv.niter                     = 500  ; 
    par.inv.burnin                    = 100  ; 
    par.inv.cooloff                   = 80   ; 
    par.inv.nchains                   = 2    ;  
elseif strcmp(stamp, 'all_fast'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.inv.niter                     = 300  ; 
    par.inv.burnin                    = 100  ; 
    par.inv.cooloff                   = 70  ; 
    par.inv.nchains                   = 1    ; 
elseif strcmp(stamp, 'all_faster'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P', 'SW_HV'}; 
    par.inv.niter                     = 500  ; 
    par.inv.burnin                    = 100  ; 
    par.inv.cooloff                   = 80   ; 
    par.inv.nchains                   = 4    ; 

    
elseif strcmp(stamp, 'sage_gage'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P'}; % No HV 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 16000; 
    par.inv.burnin                    = 5000 ; 
    par.inv.cooloff                   = 4000 ; 
    par.inv.nchains                   = 12   ; % Less chains. So we can run more stations. 
elseif strcmp(stamp, 'sage_gage_faster'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P'}; % No HV 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 700; 
    par.inv.burnin                    = 200 ; 
    par.inv.cooloff                   = 100 ; 
    par.inv.nchains                   = 12   ; % Less chains. So we can run more stations. 
    
    
    
elseif strcmp(stamp, 'most_fast'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P'}; % No HV 
    par.mod.starting.HKappa.startAtHK = true;   
    par.inv.niter                     = 500 ; 
    par.inv.burnin                    = 100 ; 
    par.inv.cooloff                   = 80  ; 
    par.inv.nchains                   = 12  ; % Less chains. So we can run more stations. 
elseif strcmp(stamp, 'most_003'); 
    par.inv.datatypes = {'SW_Ray_phV', 'SW_Lov_phV', 'RF_Sp_ccp', 'HKstack_P'}; % No HV 
    par.mod.starting.HKappa.startAtHK = true ;   
    par.inv.niter                     = 16000; 
    par.inv.burnin                    = 5000 ; 
    par.inv.cooloff                   = 4000 ; 
    par.inv.nchains                   = 12   ; % Less chains. So we can run more stations. 


% Sp receiver function, upweight at depth. 
elseif strcmp(stamp, 'ccp_weights_1'); 
    par.inv.datatypes                 = {'RF_Sp_ccp'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only; 
    par.datprocess.CCP.weight_depth_val = [-10,1 ; 6371,1]; % Constant weight. 1. 
elseif strcmp(stamp, 'ccp_weights_2'); 
    par.inv.datatypes                 = {'RF_Sp_ccp'};
    par.inv.niter                     = niter_only  ; 
    par.inv.burnin                    = burnin_only ; 
    par.inv.cooloff                   = cooloff_only; 
    par.inv.nchains                   = nchains_only; 
    par.datprocess.CCP.weight_depth_val = [-10,0.3 ; 50,0.3 ; 100,1 ; 6371,1]; % Upweight at depth. This gets normalized.  
end

inv = par.inv; 

end


