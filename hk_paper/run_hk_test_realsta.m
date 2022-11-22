% close all 
% clear all % Cannot clear all, because then we allways reset network_manual and station_manual below. 
%% Setup
run('../a0_STARTUP_BAYES.m')

addpath('/Users/brennanbrunsvik/Documents/repositories/hk_anis'); %brb TODO need to move this. 

fig_path = sprintf('%s/../figures/hk_paper', pwd()) ; 
fhand_figname = @(zmoh, k, thisfig, frmt)sprintf(...
    '%s/z%3.0f_k%3.0f_%s.%s',fig_path, zmoh*10, k*100, thisfig, frmt); % Convenient function to make figure names. Get rid of decimals. 

%%
stafold = sprintf('%s/ta_kmsc/',pwd()); 

model = load(sprintf('%sfinal_model.mat',stafold)).final_model; 
par   = load(sprintf('%spar.mat'  ,stafold)).par; 
predata = load(sprintf('%sfinal_predata.mat'  ,stafold)).final_predata; 
waves = predata.HKstack_P.waves; 

%% Extract relevant variables. 
tt = waves.tt; 
rf = waves.rf; 

z     = model.Z; 
vs    = model.VSav; 
vp    = model.VPav; 
xi    = model.Sanisav/100 + 1; 
phi   = model.Panisav/100 + 1; 
eta   = ones(size(xi)); 
rho   = model.rhoav; 
rayp  = predata.HKstack_P.waves.rayParmSecDeg; 

zmoh  = model.zmohav; 
vpvs  = model.vpvsav; 

save(sprintf('%smodel_hktest.mat',stafold), 'z', 'vs', 'vp', 'xi', 'phi', 'eta', ...
    'rho', 'zmoh', 'vpvs', 'tt', 'rf', 'rayp'); 