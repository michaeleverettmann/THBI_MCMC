
%% Description
%  This script sets all parameters that define the model and inversion
%  for the true inversion with the real data
% 
% % Inversion parms
% inv = struct(    'verbose',false                 ,... % option to spit out more information+plots
%                  'niter',1000                   ,... % Number of iterations
%                  'burnin',7000                   ,... % don't record results before burnin iterations
%                  'cooloff',4000                  ,... % # of iterations over which temperature declines as erf
%                  'tempmax',5                     ,... % maximum multiple of all standard deviations
%                  'saveperN',30                   ,... % save only every saveperN iterations       
%                  'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
%                  'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
%                  'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'nchains',10                    ,... % number of chains to start in parallel
%                  'Nsavestate',2500               ,... % Niter per which the state of the parallel inversion is saved in .mat file
%                  'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
%                  'BWclust',1                     ,... % option to use only one c x             
%                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV'}});  
%                                 % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
%                                 %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
%                                 %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
%                                 %          'HKstack_x' with x='P'

% disp('NOT REAL SYNTHETIC\nUsing fast, debuging options (few iterations). See: bays_inv_parms.m')                                
% inv = struct(    'verbose',true                 ,... % option to spit out more information+plots
%                  'niter',200                    ,... % Number of iterations
%                  'burnin',20                    ,... % don't record results before burnin iterations
%                  'cooloff',100                    ,... % # of iterations over which temperature declines as erf
%                  'tempmax',5                     ,... % maximum multiple of all standard deviations
%                  'saveperN',1                   ,... % save only every saveperN iterations    % bb2021.09.14 savig each one, since I have 100 iterations, this way we can still do probability math (taking the 5 most poorly performing models... otherwise, we get code errors later on).    
%                  'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
%                  'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
%                  'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'nchains',1                    ,... % number of chains to start in parallel
%                  'Nsavestate',25               ,... % Niter per which the state of the parallel inversion is saved in .mat file
%                  'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
%                  'BWclust',1                     ,... % option to use only one c x             
%                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV'}})  
%                                 % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
%                                 %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
%                                 %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
%                                 %          'HKstack_x' with x='P'
                                
% % % disp('100 iterations, 4 chains. See: bays_inv_parms.m')                                
% inv = struct(    'verbose',false                 ,... % option to spit out more information+plots
%                  'niter',8000                    ,... % Number of iterations
%                  'burnin',2000                    ,... % don't record results before burnin iterations
%                  'cooloff',1500                    ,... % # of iterations over which temperature declines as erf
%                  'tempmax',5                     ,... % maximum multiple of all standard deviations
%                  'saveperN',10                   ,... % save only every saveperN iterations    % bb2021.09.14 savig each one, since I have 100 iterations, this way we can still do probability math (taking the 5 most poorly performing models... otherwise, we get code errors later on).    
%                  'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
%                  'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
%                  'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'nchains',15                    ,... % number of chains to start in parallel
%                  'Nsavestate',100               ,... % Niter per which the state of the parallel inversion is saved in .mat file
%                  'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
%                  'BWclust',1                     ,... % option to use only one c x             
%                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV', 'SW_HV'}})  
%                                 % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
%                                 %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
%                                 %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
%                                 %          'HKstack_x' with x='P'

% disp('NOT REAL SYNTHETIC\nUsing fast, debuging options (few iterations). See: bays_inv_parms.m')                                
% inv = struct(    'verbose',false                 ,... % option to spit out more information+plots
%                  'niter',300                    ,... % Number of iterations
%                  'burnin',50                    ,... % don't record results before burnin iterations
%                  'cooloff',40                    ,... % # of iterations over which temperature declines as erf
%                  'tempmax',5                     ,... % maximum multiple of all standard deviations
%                  'saveperN',10                   ,... % save only every saveperN iterations    % bb2021.09.14 savig each one, since I have 100 iterations, this way we can still do probability math (taking the 5 most poorly performing models... otherwise, we get code errors later on).    
%                  'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
%                  'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
%                  'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'nchains',1                    ,... % number of chains to start in parallel
%                  'Nsavestate',10               ,... % Niter per which the state of the parallel inversion is saved in .mat file
%                  'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
%                  'BWclust',1                     ,... % option to use only one c x  
%                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV', 'SW_HV'}})  
% %                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV'}})  
%                                 % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
%                                 %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
%                                 %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
%                                 %          'HKstack_x' with x='P'


% disp('NOT REAL SYNTHETIC\nUsing fast, debuging options (few iterations). See: bays_inv_parms.m')                                
% inv = struct(    'verbose',false                 ,... % option to spit out more information+plots
%                  'niter',300                    ,... % Number of iterations
%                  'burnin',50                    ,... % don't record results before burnin iterations
%                  'cooloff',40                    ,... % # of iterations over which temperature declines as erf
%                  'tempmax',5                     ,... % maximum multiple of all standard deviations
%                  'saveperN',10                   ,... % save only every saveperN iterations    % bb2021.09.14 savig each one, since I have 100 iterations, this way we can still do probability math (taking the 5 most poorly performing models... otherwise, we get code errors later on).    
%                  'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
%                  'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
%                  'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'nchains',3                    ,... % number of chains to start in parallel
%                  'Nsavestate',10               ,... % Niter per which the state of the parallel inversion is saved in .mat file
%                  'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
%                  'BWclust',1                     ,... % option to use only one c x  
%                  'datatypes',{{'SW_HV'}})  
% %                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV'}})  
%                                 % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
%                                 %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
%                                 %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
%                                 %          'HKstack_x' with x='P'
%                                 

disp('NOT REAL SYNTHETIC\nUsing fast, debuging options (few iterations). See: bays_inv_parms.m')                                
inv = struct(    'synthTest',true                ,...
                 'verbose',false                 ,... % option to spit out more information+plots
                 'niter',500                    ,... % Number of iterations
                 'burnin',100                    ,... % don't record results before burnin iterations
                 'cooloff',100                    ,... % # of iterations over which temperature declines as erf
                 'tempmax',5                     ,... % maximum multiple of all standard deviations
                 'saveperN',10                   ,... % save only every saveperN iterations    % bb2021.09.14 savig each one, since I have 100 iterations, this way we can still do probability math (taking the 5 most poorly performing models... otherwise, we get code errors later on).    
                 'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
                 'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
                 'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'nchains',16                   ,... % number of chains to start in parallel
                 'Nsavestate',10               ,... % Niter per which the state of the parallel inversion is saved in .mat file
                 'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
                 'BWclust',1                     ,... % option to use only one c x  
                  'datatypes',{{'HKstack_P','RF_Sp','SW_Ray_phV','SW_Lov_phV', 'SW_HV'}})  
%                  'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV'}})  
                                % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
                                %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
                                %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
                                %          'HKstack_x' with x='P'
                                
% inv = struct(    'synthTest',true                ,...
%                  'verbose',false                 ,... % option to spit out more information+plots
%                  'niter',16000                     ,... % Number of iterations
%                  'burnin',5000                    ,... % don't record results before burnin iterations
%                  'cooloff',4000                    ,... % # of iterations over which temperature declines as erf
%                  'tempmax',5                     ,... % maximum multiple of all standard deviations
%                  'saveperN',30                   ,... % save only every saveperN iterations    % bb2021.09.14 savig each one, since I have 100 iterations, this way we can still do probability math (taking the 5 most poorly performing models... otherwise, we get code errors later on).    
%                  'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
%                  'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
%                  'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
%                  'nchains',16                    ,... % number of chains to start in parallel
%                  'Nsavestate',100                ,... % Niter per which the state of the parallel inversion is saved in .mat file
%                  'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
%                  'BWclust',1                     ,... % option to use only one c x             
%                   'datatypes',{{'HKstack_P','RF_Sp','SW_Ray_phV','SW_Lov_phV', 'SW_HV'}})  
%                                 % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
%                                 %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
%                                 %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
%                                 %          'HKstack_x' with x='P'


profileRun = true; if profileRun; fprintf('\n\nDoing an mpi profile run.\n\n'), end
                                
%% Model parms
modl = struct([]);

modl(1).nstas = 1;
modl.maxz = 300;                                      % maximum depth in model from ref ellipsoid, km
modl.maxkz = 250;                                     % maximum depth of deepest non-basal knot, km
modl.dz = 2;                                        % depth spacing of model, km

modl.starting.HKappa.startAtHK = true; % Starting model at maximum of hkappa stack

modl.sed = struct(...
    ... thickness of the sediments
                     'hmax',0.0                  ,... %5 max sed layer thickness, km
                     'hmin',0.0                  ,... %0 min sed layer thickness, km
                     'hstd',0.0                  ,... % std of sed layer thickness for perturbation, km
    ... shear velocity of the sediments
                     'vsmax',3.3                 ,... % max sed velocity, km/s
                     'vsmin',2.8                 ,... % min sed velocity, km/s
                     'vsstd',0.00                 );  % std of sed velocity for perturbation, km/s

modl.crust = struct(...
    ... thickness of the crust
                     'hmax',70                   ,... %60 max xtal crust thickness, km % Shen and Ritzwoller 2016 Don't show anything above ~55 km in eastern US. bb2021.10.26
                     'hmin',25                   ,... %10 min xtal crust thickness, km % bb2021.10.26 Going VERY shallow because we will go offshore, expecting down to about 7 km thickness (Shuck et al., 2019). We don't want priors to thicken the offshore crust. 
                     'hstd',2.5                    ,... % std of xtal crust thickness, for perturbation, km
                ... gaussian prior probability for crust thickness - mean=30, std=10
                     'h_pprior',@(h) 1,...exp(-(h-30).^2/4.^2),... % prior probability 
    ... shear velocity in the crust
                     'vsmax',4.4                 ,...4.5 % max crust spline velocity, km/s bb2021.10.26 changes from 4.3 to 4.4... Shen and Ritzwoller Fig 12 shows lower crustal velocity going above 4.2 
                     'vsmin',2.5                 ,...3.3 % min crust spline velocity, km/s bb2021.10.26 ?? Not sure about this one. There must be a theoretical limit. Shen2016 shows SUPER low velocities in Gulf of Mexico, which is clearly a consequence of sediment. North/East of that, values are higher than 2.8. 
                     'vsstd',0.08                ,... % std of crust spline velocity for perturbation, km/s
    ... Vp/Vs in the crust
                     'vpvsmax',2.1               ,...1.9 % max crust vpvs ratio bb2021.10.26 Changed from 1.9 because station LSCT seemed to want a very high Vp/Vs! Looks like Jon also went with 2.1. 
                     'vpvsmin',1.6               ,...1.65 % min crust vpvs ratio
                     'vpvsstd',0.01              ,... % std of crust vpvs ratio for perturbation, km/s
                ... gaussian prior probability for VpVs - mean=1.8, std=0.05
                     'vpvs_pprior',@(vpvs) 1,...exp(-(vpvs-1.7).^2/0.03.^2),... % prior probability 
    ... Xi in the crust
                     'ximax',1.2                 ,...1.05 % min crust Vs radial anis value % bb2021.10.26 Changed from 1.1 because Station LSCT suggests we want wider psi bounds. 
                     'ximin',0.75                 ,...1.00 % min crust Vs radial anis value % bb2021.10.26 Changed from 0.9 to 0.75 because station LSCT wanted SUPER low psi. THIS MIGHT BE A CONSEQUENCE OF HAVING MANTLE PSI = 1. .9 to 1.1 was suggested by GLOBAL and thus lower resolution compilation of Porrit et al., 2021. 
                     'xistd',0.01                ,... % std of crust Vs radial anis value
    ... knots in the crust
                     'kdstd',2                   ,... % std of knot movement, for perturbation, km
                     'kmax',7                    ,... % max number of spline knots in crust (inc ends)
                     'kmin',2                    );  % min number of spline knots in crust (inc ends)

modl.mantle = struct(...
    ... shear velocity in the mantle
                     'vsmax',5.1                 ,...4.9 % max mantle spline velocity, km/s bb2021.10.26 I haven't seen more than about 4.9 in Shen and Ritzwoller or even full waveform models. BUT the pdfs in Shen and Ritzwoller did clip at 4. 9 (Figure 8). So 5.1 seems good. 
                     'vsmin',3.7                 ,...3.7 % min mantle spline velocity, km/s bb2021.10.26 Long et al., 2021 harrisonburg anomaly compilation shows velocities all above about 4.2 km/s... so 3.7 can definatley handle mantle anomalies like this. 
                     'vsstd',0.08                ,... % std of mantle spline velocity for perturbation, km/s
    ... Xi in the mantle
                     'ximax',1.0                 ,...1.05 % min mantle Vs radial anis value. bb2021.10.26 I don't understand why we set this min and max to zero. I guess we need Long period Love and Rayleigh wave data to constrain it. 
                     'ximin',1.0                 ,...1.00 % min mantle Vs radial anis value
                     'xistd',0                   ,... % std of mantle Vs radial anis value
    ... knots in the mantle
                     'kdstd',4                   ,... % std of knot movement, for perturbation, km
                     'kmax',15                   ,... % max number of spline knots in mantle (inc ends)
                     'kmin',5                    );  % max number of spline knots in mantle (inc ends)

modl.data = struct('prior_sigma',struct(                 ... % PRIOR
                  	 'BW',struct(                 ... %  Body waves
                    	'Ps',struct(              ... %   P-s data
                           'def',0.3             ,... %    default
                           'lo',0.2              ,... %    low-f
                           'cms',0.3)            ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',0.2             ,... %    default
                           'lo',0.1))            ,... %    low-f
                  	 'RF',struct(                 ... %  Receiver functions
                    	'Ps',struct(              ... %   P-s data
                           'def',0.2             ,... %    default
                           'lo',0.1              ,... %    low-f
                           'cms',0.3)            ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',0.2             ,... %    default
                           'ccp',0.1             ,... %    ccp stack
                           'lo',0.1))            ,... %    low-f
                  	 'HKstack',struct(            ... %  H-K stack
                    	   'P',.3)              ,... %    P combination
                  	 'SW',struct(                 ... %  Surface waves
                    	'Ray',struct(             ... %   Rayleigh waves
                           'phV',0.05            ,... %    phase velocities
                           'grV',0.06)           ,... %    group velocities
                    	'HV',struct(             ... %   Rayleigh wave ellipticity
                           'HVr',0.06)           ,... %    HV ratio
                    	'Lov',struct(             ... %   Love waves
                           'phV',0.05            ,... %    phase velocities
                           'grV',0.06)))         ,... %    group velocities
                                                  ...  
                  'min_sigma',struct(             ... % PRIOR
                  	 'BW',struct(                 ... %  Body waves
                    	'Ps',struct(              ... %   P-s data
                           'def',1e-2            ,... %    default
                           'lo',5e-3             ,... %    low-f
                           'cms',1e-3)           ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',1e-3            ,... %    default
                           'lo',1e-3))           ,... %    low-f
                  	 'RF',struct(                 ... %  Receiver functions
                    	'Ps',struct(              ... %   P-s data
                           'def',1e-2            ,... %    default
                           'lo',1e-2             ,... %    low-f
                           'cms',1e-2)           ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',1e-2            ,... %    default
                           'ccp',1e-2             ,... %    ccp stack
                           'lo',1e-2))           ,... %    low-f
                  	 'HKstack',struct(            ... %  H-K stack
                    	   'P',0.06)             ,... %    P combination
                  	 'SW',struct(                 ... %  Surface waves
                    	'Ray',struct(             ... %   Rayleigh waves
                           'phV',1e-4            ,... %    phase velocities
                           'grV',1e-4)           ,... %    group velocities
                    	'HV',struct(             ... %   Rayleigh wave ellipticity
                           'HVr',1e-4)           ,... %    HV ratio
                    	'Lov',struct(             ... %   Love waves
                           'phV',1e-4            ,... %    phase velocities
                           'grV',1e-4)))         ,... %    group velocities
                                                  ...  
                  'logstd_sigma',0.05,            ...
                  'deg_of_freedom',struct(       ...
                      'h_kappa', 15)) ;   % LOGSTD

                 
%% Forward calc. parms
forc = struct(      'mindV',0.05                 ,... % min delta Vs for layerising
                    'nsamps',2^11                ,... % number of samples (more means longer run time)
                    'PSVorZR','PSV'             ,... % whether to rotate data into PSV or keep in ZR
                    'synthperiod',2              );  % period for propmat response
                
%% Data processing parms
datprocess=struct( 'normdata',true               ,... % normalise data in processing
                   'decdata',false               ,... % decimate data in processing
                   'clipmain',false              ,... % whether to clip the main phase with a taper
                   'clipdaughter',true           ,... % whether to clip the daughter component at the timing of the main phase (to account for improper rotation) == CHEAT
                   'PSV',true                    ,... % PSV (if tru) or ZR
                   'Ps',struct(                   ... % P-s data
                      'Twin'                     ,... %   time window    
                      struct('def',[-2 10],       ... %     default	 
                             'cms',[-2 25])      ,... %     crustal multiples 	  
                      'filtf'                    ,... %   filter frequencies    
                      struct('def',[1e3 1e-3]    ,... %     default	[fhi flo] 
                             'lo',[1/3  1e-3]))  ,... %     low-f [fhi flo]
                   'Sp',struct(                   ... % P-s data
                      'Twin'                     ,... %   time window    
                      struct('def',[-35 4])      ,... %     default 	  
                      'filtf'                    ,... %   filter frequencies    
                      struct('def',[1e3 1e-3]    ,... %     default	[fhi flo] 
                             'lo',[1/5  1e-3]))  ,...
                   'CCP',struct(                  ... % CCP data from Emily
                      'parent_zw',30             ,... %   depth-width of "parent" pulse for CCP    
                      'rayp_S',12.303            ,... %   S wave ray parameter (s/deg) for gcarc=0, edep=0 in iasp91  
                      'surv_Vp_vs',[6.1 3.55]    ,... %   [VP, VS] surface velocity values for P-SV rotation
                      'taperz',10               ,... %   taper width at the edges of the Zwin
                      'Zwin'                     ,... %   depth window    
                      struct('def',[20 250]))    ,...
                   'HKappa',struct(              ...
                       'min_error', 0.002,           ... % Add this much "error" to h-kappa stacks (error of 0 can result in sigma inverting improperly)
                       'scale_error', 1,           ... % Multiply h-kappa error by this constant. Sigma needs to be scaled accordingly. If using 100, we can think of it like percent. 
                       'weightDistanceMax', 0,   ... % At start of burnin, gives 0 to 1 weight toward the (scaled) Euclidian distance from HKappa energy maximum. In otherwords, this tends toward disregarding the actual energy value, and pays attention to its position. 
                       'weightDistanceMin', 0));     % At end of burnin, give this much weight 0 to 1 to distance from h,k where E is max.     	 
%                     
         
%% Model Conditions
cond = struct(  'pos_moho',         true         ,... % No negative moho jumps
                'pos_sed2basement', true         ,... % No negative sed bottom jumps
                'pos_seddV',        true         ,... % Monotonic increase of V(p/s) in sediments
                'pos_crustdV',      false         ,... % Monotonic increase of V(p/s) in crust
                'nobigdVmoh',       true         ,... % No Vs moho jumps exceeding 30%
                'no_moho_sharpgrad',true         ,... % No strong or negative gradients on either side of the Moho
                'noVSgt49',         true         );  % No VS exceeding 4.9 km/s
            
%% Synthetic data parms 
synth = struct( 'gcarcs',[70]                 ,... % average gcarc
                'samprate',10                    ,... % sample rate.
                'noisetype','gaussian'           ,... % noisetype - 'gaussian' or 'real' (in which case drawn from station speficied
                'noise_sigma_SW_Ray',0.015        ,... %0.03 std for random added noise for SWs
                'noise_sigma_SW_Lov',0.015        ,... %0.03 std for random added noise for SWs
                'noise_sigma_SW_HV',0.005        ,... %0.03 std for random added noise for SWs
                'noise_sigma_BW_Sp',0.009        ,... %0.02 std for random added noise for SpRFs
                'noise_sigma_BW_Ps',0.012        ,... %0.02 std for random added noise for PsRFs 0.012
                'noise_sigma_RF_Sp',0.009        ,... %0.02 std for random added noise for SpRFs
                'noise_sigma_RF_Ps',0.012        ,... %0.02 std for random added noise for PsRFs 0.012
                'surf_Vp_Vs',[6.1 3.55]          ,... % [VP, VS] surface velocity values - if empty, uses True vals % bb2022.02.08 Not sure what these are. Different sets of values are used in z0_SYNTH_MODEL...
                'SW_Ray_phV_periods',logspace(log10(6),log10(167),22)',...  % Rayleigh wave phV periods
                'SW_Ray_grV_periods',logspace(log10(6),log10(40),10)',...  % Rayleigh wave phV periods
                'SW_Lov_phV_periods',logspace(log10(6),log10(40),10)',...  % Love wave phV periods
                'SW_HV_periods',logspace(log10(8),log10(90),11)',...  % Rayleigh wave HV periods
                'synthperiod',2                  ,...  % period for propmat synth
                'nsamps',[]                      );  % number of samples. 
warning('No ps noise')
%---RF parameters---%
RFparms = struct([]);
% 'IDRF' for iterative time domain method, or 'ETMTM' for Extended time multitaper method
RFparms(1).method = 'IDRF';
%input: time domain deconvolution
RFparms.gauss_t = 0.5;
RFparms.accept_mis=1e-10; %accepted misfit
RFparms.itmax=50; %number of iterations
%input: freq domain deconvolution
RFparms.TB = 1.5;%1.5; %period*bandwidth
RFparms.NT = 2;%2; %number of tapers
RFparms.tag = 'synth'; %data or synth (synthetic) - Zach uses 'synth' & this gives better results
RFparms.Poverlap = 0.95; %fractional window overlap
%=================================%
synth.RFparms = RFparms;


%% Results parms
res = struct('zatdep',[5:5:modl.maxz]');
            
%% Bundle together
par = struct('inv',inv,'mod',modl,'conditions',cond,'forc',forc,'synth',synth,'datprocess',datprocess,'res',res);
