 %%
clc; clear;

%%% Misc startup stuff shouldn't have to change. 
addpath('matlab_to_hv_kernel') ; 
addpath('functions') ;  
getPaths=@()struct('HV_ellipticity', [pwd '/HV_ellipticity'],...  
                   'timeout',         '/usr/local/bin/timeout',... % TOSHIRO: you may need to change your path to the shell "timeout" function.
                   'THBIpath',       [pwd],...
                   'execPath',       [pwd],... 
                   'ramDrive',        [pwd]); 
paths = getPaths(); 
matlab.io.saveVariablesToScript('pathsAutoGen.m', 'paths');
global prem_anisotropic prem_isotropic
prem_anisotropic = prem_perfect('SPVW',0.25);
prem_isotropic = prem;
resdir_fig = pwd; % Any junk folder. 
%%% End Misc startup stuff shouldn't have to change. 


%% Load a model and corresponding variables that illustrate kernel problems. I can send more if you want. 
load('TA.R54A_C_ii_1346_likelihood_drop_hv_workspace.mat'); % TOSHIRO Can change this to see other examples of models/kernels. 
par.inv.datatypes = {'SW_HV'}; % Not sending code for other datatypes. 
par.res.resdir = resdir_fig; 

model0 = Kbase.modelk; % The old model, associated with the old kernels i.e. Kbase. 
[Kbase,predata] = b7_KERNEL_RESET(model0,Kbase,predata,ID,ii,par,1); % Make kernels for old model. 
predata = b3_FORWARD_MODEL_SW_kernel(model,Kbase,par,predata ); % This is where we use the new model "model" and the kernels from old model "Kbase" to estimate HV ratio

[KbaseNew,predataNew] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,1); % Making the kernel for HV here. "model" is the second, or new, model.
if ifforwardfail(predataNew,par) 
    error('Failed making new kernels! Fail_chain %.0f, ii=%1.0f\n',fail_chain,ii)
end
KbasePrev = Kbase; % Old kernels. 




%%%% Make plots. 
fns = {'SW_HV'}; % fns = fieldnames(misfitNew.E2); 
ifn = 1; 

LW = 1.5; 
figure(3); clf; hold on; set(gcf, 'pos', [-1215 280 894 638]); 
tiledlayout(1,3,'TileSpacing','compact'); 

%%% Plot HV ratio
nexttile(1, [1,1]); cla; hold on; box on; 
set(gca, 'ydir', 'reverse', 'LineWidth', 1.5);  
xlim([0.6, 1.2]); 
ylim([min(trudata.SW_HV.periods-1.5), max(trudata.SW_HV.periods+5)]); 
set(gca, 'yscale', 'log'); 
grid off; grid on;  % Log scale requires turning grid off before on

plot(trudata.SW_HV.HVr, trudata.SW_HV.periods, ...
    'DisplayName', 'Goal', 'LineWidth', LW*3/2, ...
    'Color', 'Green'); 

sw_hv = Kbase.HV; 
plot(sw_hv.HVr, sw_hv.periods, ...
    'DisplayName', 'Old HV (full calc)', 'LineWidth', LW, ...
    'Color', 'k');

sw_hv = KbaseNew.HV; 
plot(sw_hv.HVr, sw_hv.periods, ...
    'DisplayName', 'New HV (full calc)', 'LineWidth', LW, ...
    'Color', 'blue');

sw_hv = predata.(fn); 
plot(sw_hv.HVr, sw_hv.periods, ...
    'DisplayName', 'New HV (from kernels)', 'LineWidth', LW, ...
    'Color', 'r');

legend('Location', 'Best'); 
xlabel('H/V ratio'); ylabel('Period');

mismatch = rms(predata.(fn).HVr - KbaseNew.HV.HVr); % Mistmach between kernel and new full forward model
dC_est   = rms(predata.(fn).HVr - Kbase   .HV.HVr); % dC predicted from kernel
dC_full  = rms(Kbase   .HV.HVr  - KbaseNew.HV.HVr); % dC based on full calculations 
%%% End HV ratio

%%% Plot models. 
nexttile(2, [1,2]); hold on; box on; 
set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
plot([Kbase.modelk.VS, Kbase.modelk.VP], Kbase.modelk.z, ...
    'DisplayName', 'Previous kernel model', 'LineWidth', LW, ...
    'Color', 'k');
plot([KbaseNew.modelk.VS,KbaseNew.modelk.VP], KbaseNew.modelk.z, ...
    'DisplayName', 'New kernel model', 'LineWidth', LW, ...
    'Color', 'blue');
leg=legend('Location', 'Best'); 
xlabel('V'); ylabel('Depth (km)'); 
ylim([-5, max(Kbase.modelk.z)]); 
%%% End models

%%% Title showing some info about how much the kernel and new full forward modelled HV ratios agree. 
sgtitle(sprintf(...
    ['HV kernel versus full calculation.\nIteration=%6.0f. ',...
    'Model pertubation norm = %3.2f. ',...
    'Change in log-likelihood = %3.5f\n',...
    'RMS of: Kernel error: %1.4f -- dC kernel: %1.4f -- dC full model: %1.4f.'],... 
    par.ii, ptbnorm, logL_diff,...
    mismatch, dC_est, dC_full)); 

exportgraphics(gcf, sprintf('%s/%s_ii_%1.0f_likelihood_drop_hv.pdf',...
    par.res.resdir, par.res.chainstr, ii ) ); % %!%! Temporary 

%%% Plot sensitivity kernels from old and new Kbase. 
plot_sensitivity_kernel_hv(Kbase.HV.KHV   , 'dat', Kbase.HV   ,... 
    'model',Kbase.modelk   ,...
    'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_old.pdf', ...
    par.res.resdir, par.res.chainstr, par.ii) ); % Old kernels

plot_sensitivity_kernel_hv(KbaseNew.HV.KHV, 'dat', KbaseNew.HV, ...
    'model',KbaseNew.modelk,...
    'filename',sprintf('%s/%s_ii_%1.0f_sensitivity_hv_new.pdf',...
    par.res.resdir, par.res.chainstr, par.ii) ); % New kernels
%%% End sensitivity kernels.

save( sprintf('%s/%s_ii_%1.0f_likelihood_drop_hv_workspace.mat',...
    par.res.resdir, par.res.chainstr, ii ) ); 
