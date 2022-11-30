% close all 
% clear all % Cannot clear all, because then we allways reset network_manual and station_manual below. 
%% Setup
run('../a0_STARTUP_BAYES.m')

addpath('/Users/brennanbrunsvik/Documents/repositories/hk_anis'); %brb TODO need to move this. 

fig_path = sprintf('%s/../figures/hk_paper', pwd()) ; 
fhand_figname = @(zmoh, k, thisfig, frmt)sprintf(...
    '%s/z%3.0f_k%3.0f_%s.%s',fig_path, zmoh*10, k*100, thisfig, frmt); % Convenient function to make figure names. Get rid of decimals. 

proj = struct('name', 'SYNTHETICS'); % bb2021.08.04 changed from EXAMPLE because I don't have the example data files. %,'EXAMPLE');
paths = getPaths(); 
proj.STAinversions = paths.STAinversions; 
proj.dir = [fileparts(mfilename('fullpath')),'/'];
proj.STAinversions = paths.STAinversions; ; % [proj.dir,'inversion_results/'];
save([proj.dir,'project_details.mat'],'proj')

wd = pwd; addpath(wd);
cd(proj.dir);

%% specify details of this run
generation = 0; % generation of solution and data processing
gc = '';
BWclust = '';
STAMP = 'hk_paper_1';
onesta = '';

%% put parameters in place 
global run_params

% you can obviously adapt these (likely in some loop over stations etc.) to be appropriate for your dataset
run_params.projname = proj.name; % from above
run_params.gc = gc; % great circle distance of body wave data, if relevant
run_params.BWclust = BWclust; % cluster of BW data, if relevant
run_params.datN = generation; % processing iteration, if relevant
run_params.STAMP = STAMP; % NEED - some identifier for this inversion run
run_params.overwrite = 1; % do you want to overwrite previous results?
% % % if ~ (exist('network_manual', 'var') && exist('station_manual', 'var')) ; 
network_manual = 'testnwk'; 
station_manual = 'simple_layers_1'; %
fprintf('\nReseting to %s.%s\n',network_manual,station_manual)
% % % end
run_params.sta = station_manual; % name of station
run_params.nwk = network_manual; % name of network


global run_params
paths = getPaths(); 

projname = run_params.projname;
sta = run_params.sta;
nwk = run_params.nwk;
gc = run_params.gc;
BWclust = run_params.BWclust;
datN = run_params.datN;
STAMP = run_params.STAMP;
overwrite = run_params.overwrite;
global projdir TRUEmodel
projdir = [paths.THBIpath,'/',projname,'/'];
cd(projdir);
run([paths.THBIpath,'/a0_STARTUP_BAYES']);
load('project_details'); %TODO_STATION_NETWORK bb2021.11.12
addpath([proj.dir,'matguts/']);

%% PARMS
run parms/bayes_inv_parms
[par, inv] = update_bayes_inv_parms(par, STAMP); % Modify this function to make different tests. 


if strcmp(projname,'SYNTHETICS')
    par.stadeets = struct('sta',sta','nwk',nwk'); 
	noisesta = 'RSSD';
	noisenwk = 'IU';
	noisegcarcs = [73,38];
	noiseshape = 'real'; % 'white' or 'real'
	noiseup = 0.5; % factor to increase real noise
end

if strcmp(projname,'SYNTHETICS') || strcmp(projname,'LAB_tests')
    par.synth.noise_sta_deets = struct('datadir',['/Volumes/data/THBI/US/STAsinv/',noisesta,'_dat20/'],...
                         'sta',noisesta,'nwk',noisenwk,'gc',noisegcarcs,'noiseup',noiseup,'noiseshape',noiseshape);
end

par.inv.BWclust = BWclust;
ifsavedat = false;

%% get saving things ready
par.proj = proj;
avardir = sprintf('%s%s_%s_dat%.0f/',par.proj.STAinversions,sta,nwk,datN);
resdir = [avardir,STAMP];
if ~exist(resdir,'dir'), try mkdir(resdir); catch, error('Looks like no path to output directory - is it mounted?'); end, end

par.data = struct('stadeets',struct('sta',sta,'nwk',nwk,'Latitude',[],'Longitude',[]),...
                  'gc',gc,'datN',datN,'avardir',avardir);

par.res.STAMP = STAMP;
par.res.resdir= resdir;
par.res = orderfields(par.res,{'STAMP','resdir','zatdep'});

%% Get some directories ready. 
% Switch to execution folder, to make synthetic data. 
prev_dir = pwd(); 
cd(paths.ramDrive); % Execute everything from a folder in ram for major speedup. 
mkdir([nwk '_' sta]); cd([nwk '_' sta]); % Go to station specific folder to keep things clean . TODO just to cd once. 

%% HK tests, analysis, starts here. 
xi_a = [0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15]'; 
xi_true = 0.85; 
xi_a = sort(unique([xi_a; xi_true])); 
nxi = length(xi_a); 
i_xi_true = find(xi_a == xi_true); 

% % % par.datprocess.kNum = 1001; 
% % % par.datprocess.hNum = 1000; 
par.datprocess.kNum = 101; 
par.datprocess.hNum = 100; 
warning('Less hk resolution right now')

Exi_all = cell([length(xi_a), 1]); 
E00_all = cell([length(xi_a), 1]); 
% t_predxi_all = cell([length(xi_a), 1]); 
% t_pred00_all = cell([length(xi_a), 1]); 
t_pred_xi_best_all = zeros(length(xi_a), 3); 
t_pred_xi_noan_all = zeros(length(xi_a), 3); 

hmax_all_noan = zeros(nxi,1); 
kmax_all_noan = zeros(nxi,1); 

rf_all = cell(length(xi_a),1); 

herr = zeros(nxi, 1); 
kerr = zeros(nxi, 1); 

global TRUEmodel %Unfortunately this was already used as a global model
each_model = {}; 

for ixi = 1:length(xi_a);  
%     for ixi = i_xi_true;  

    fprintf('Do something about xi_true.\n')

%     ztrue = 45; 
%     ktrue = 1.75; 
    xitruei = xi_a(ixi); 

    [trudata,par] = a2_LOAD_DATA_hk_test(par, 'nwk', nwk, 'sta', sta, ...
        'xi_crust', xitruei );

    each_model{ixi, 1} = TRUEmodel; 
    ztrue = TRUEmodel.zmoh; 
    ktrue = TRUEmodel.vpvs; 
    % plot_TRU_WAVEFORMS(trudata);

    
    
    H = trudata.HKstack_P.H;
    K = trudata.HKstack_P.K; 
    Exi = trudata.HKstack_P.Esum; 
    E00 = trudata.HKstack_P_noan.Esum; 
    waves = trudata.HKstack_P.waves; 
    t_predxi = trudata.HKstack_P.t_pred; 
    t_pred00 = trudata.HKstack_P_noan.t_pred; 

%     [Exi_max, iemax] = max(Exi, [], 'all'); %linear index. 
    [ikmax, ihmax] = find(E00 == max(E00 ,[], 'all')); 
    warning('I think I had k and h indicies backwards in the above line')
    kmax_noan = K(ikmax); 
    hmax_noan = H(ihmax); 
    herr(ixi) = hmax_noan - ztrue; 
    kerr(ixi) = kmax_noan - ktrue; 

    [ikmax, ihmax] = find(Exi == max(Exi ,[], 'all')); 
    kmax_best = K(ikmax); 
    hmax_best = H(ihmax); 
%     herr(ixi) = hmax_noan - ztrue; 
%     kerr(ixi) = kmax_noan - ktrue; 
    
    t_pred_xi_best  = zeros(1, 3); % Get the times from anisotropic stack at true parameters
    t_pred_xi_noan  = zeros(1, 3); % Get the times from ISOTROPIC stack at true parameters. 
%     t_pred_xi_noan = zeros(1, 3);
    for it = 1:length(t_pred_xi_best); 
        t_pred_xi_best(1, it) = interpn(H, K, ...
            reshape(t_predxi(it,:,:), size(t_predxi,2), size(t_predxi,3)) , ...
            ztrue, ktrue, 'cubic'); 
        t_pred_xi_noan(1, it) = interpn(H, K, ...
            reshape(t_pred00(it,:,:), size(t_predxi,2), size(t_predxi,3)) , ...
            ztrue, ktrue, 'cubic'); 
    end
    t_pred_xi_best_all(ixi, :) = t_pred_xi_best; 
    t_pred_xi_noan_all(ixi, :) = t_pred_xi_noan; 
    
    rf_all{ixi} = waves.rf; %  + ixi; 
    Exi_all{ixi} = Exi; 
    E00_all{ixi} = E00; 
    hmax_all_noan(ixi) = hmax_noan; 
    kmax_all_noan(ixi) = kmax_noan; 

    hmax_all_best(ixi) = hmax_best; 
    kmax_all_best(ixi) = kmax_best; 

    %%% How much error in the non-anisotropic HK stack? 
%     interp2(H, K, E00, hmax, kmax)
    %%%
end
 
%%
figure(201); clf; hold on; 
subplot(1,2,1); hold on; 
set(gca,'ydir', 'reverse'); 
contourf(K, H, Exi_all{i_xi_true}', 30, 'EdgeAlpha', 0.1); 
ylim([25, 55]); 
xlim([1.6, 1.9]); 

%%% Percentage contours
max_hk = max(Exi_all{i_xi_true},[], 'all'); 
contour(K, H, Exi_all{i_xi_true}', [0.68, 0.95].* max_hk, 'k', 'LineWidth',1); 
%%%

size_scat = 40; 

scatter(kmax_all_noan, hmax_all_noan, size_scat, 'red', 'filled'); 
scatter(kmax_all_noan(i_xi_true), hmax_all_noan(i_xi_true), size_scat*2, 'red', 'diamond', 'filled'); 
plot(kmax_all_noan, hmax_all_noan, 'red', 'LineWidth', 1.);

scatter(kmax_all_best, hmax_all_best, size_scat, 'blue', 'filled'); 
scatter(kmax_all_best(i_xi_true), hmax_all_best(i_xi_true), size_scat*2, 'blue', 'diamond', 'filled'); 
plot(kmax_all_best, hmax_all_best, 'blue', 'LineWidth', 1.);

% text(kmax_all_noan([1,nxi])+0.005, hmax_all_noan([1,nxi]) - 0.75, string(xi_a([1,nxi])) )
text(kmax_all_noan+0.005, hmax_all_noan - 0.75, string(xi_a) )

%% 
fhand_norm = @(inval)inval ./ max(max(inval)); % Return normalized inval 

figure(301); clf; hold on; set(gcf, 'pos', [-1089 329 364 218]); 
subplot(1,1,1); hold on; 
set(gca,'ydir', 'reverse', 'LineWidth', 1.5);
grid on; 
box on; 
xlabel('\kappa'); 
ylabel('H (km)'); 
title('H-\kappa stack ignoring \xi', 'fontweight', 'normal'); 
% contourf(K, H, Exi_all{i_xi_true}', 30, 'EdgeAlpha', 0.1); 

plt_ylim = [40, 50]; 
plt_xlim = [1.6, 1.9]; 
ylim(plt_ylim); 
xlim(plt_xlim); 

% plot([ktrue, ktrue], plt_ylim + [-1, 1])
hnd_true = scatter(ktrue, ztrue, 200, [252, 3, 252]./256, 'pentagram', 'filled', ...
    'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
    'DisplayName', 'True'); 

lvl_cnt = [0.7, 0.95]; 
LW = 1; 
[~,hnd_xistart] = contour(K, H, fhand_norm(E00_all{1      }'),...
    lvl_cnt, 'r', 'LineWidth', LW, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(1      ) ),...
    'LineStyle','--'); 
[~,hnd_xi1    ] = contour(K, H, fhand_norm(E00_all{end    }'),...
    lvl_cnt, 'b', 'LineWidth', LW, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(end    ) ),...
    'LineStyle','--'); 
[~,hnd_xiend  ] = contour(K, H, fhand_norm(E00_all{xi_a==1}'),...
    lvl_cnt, 'k', 'LineWidth', LW*1.5, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(xi_a==1) ) ); 


% scatter(kmax_all_noan, hmax_all_noan, size_scat, 'red', 'filled'); 
% scatter(kmax_all_noan(i_xi_true), hmax_all_noan(i_xi_true), size_scat*2, 'red', 'diamond', 'filled'); 
% plot(kmax_all_noan, hmax_all_noan, 'red', 'LineWidth', 1.);
% text(kmax_all_noan+0.005, hmax_all_noan - 0.75, string(xi_a) )


legend([hnd_xistart, hnd_xiend, hnd_xi1, hnd_true], 'Location', 'best'); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'multiple_contour', 'pdf'), 'ContentType', 'vector'); 


%%
figure(302); clf; hold on; set(gcf, 'pos', [-1089 329 364 218]); 
subplot(1,1,1); hold on; 
set(gca,'ydir', 'reverse', 'LineWidth', 1.5);
grid on; 
box on; 
xlabel('\kappa'); 
ylabel('H (km)'); 
title('H-\kappa stack \xi corrected', 'fontweight', 'normal'); 
% contourf(K, H, Exi_all{i_xi_true}', 30, 'EdgeAlpha', 0.1); 

plt_ylim = [40, 50]; 
plt_xlim = [1.6, 1.9]; 
ylim(plt_ylim); 
xlim(plt_xlim); 

% plot([ktrue, ktrue], plt_ylim + [-1, 1])
hnd_true = scatter(ktrue, ztrue, 200, [252, 3, 252]./256, 'pentagram', 'filled', ...
    'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
    'DisplayName', 'True'); 

lvl_cnt = [0.7, 0.95]; 
LW = 1; 
[~,hnd_xistart] = contour(K, H, fhand_norm(Exi_all{1      }'),...
    lvl_cnt, 'r', 'LineWidth', LW, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(1      ) ),...
    'LineStyle','--'); 
[~,hnd_xi1    ] = contour(K, H, fhand_norm(Exi_all{end    }'),...
    lvl_cnt, 'b', 'LineWidth', LW, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(end    ) ),...
    'LineStyle','--'); 
[~,hnd_xiend  ] = contour(K, H, fhand_norm(Exi_all{xi_a==1}'),...
    lvl_cnt, 'k', 'LineWidth', LW*1.5, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(xi_a==1) ) ); 


% scatter(kmax_all_noan, hmax_all_noan, size_scat, 'red', 'filled'); 
% scatter(kmax_all_noan(i_xi_true), hmax_all_noan(i_xi_true), size_scat*2, 'red', 'diamond', 'filled'); 
% plot(kmax_all_noan, hmax_all_noan, 'red', 'LineWidth', 1.);
% text(kmax_all_noan+0.005, hmax_all_noan - 0.75, string(xi_a) )


% legend([hnd_xistart, hnd_xiend, hnd_xi1, hnd_true], 'Location', 'best'); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'multiple_contour_hkcor', 'pdf'), 'ContentType', 'vector'); 



%%
figure(202); clf; hold on; 
set(gcf, 'pos', [1060 564 445 235]); 
set(gca, 'LineWidth', 1.5, 'XGrid', 'on', 'XMinorTick', 'on'); box on; %grid on; 
xlabel('Time (s)'); 
title('Phase timing', 'FontWeight','normal'); 
set(gca, 'YTick', []); 
xlim([-3, 30])
yshift_const = 0.075; 

for ixi = 1:nxi
    yshift = ixi * yshift_const; 
    rf = rf_all{ixi}; 

    t_pred_xi_best = t_pred_xi_best_all(ixi,:)'; 
    t_pred_xi_noan = t_pred_xi_noan_all(ixi,:)'; 

    hnd_t_xi = scatter(...
        t_pred_xi_best', yshift + interp1(waves.tt, rf, t_pred_xi_best, 'cubic'),...
        40, 'blue', 'filled') % If using true parameters and anisotropic stack
    hnd_t_00 = scatter(...
        t_pred_xi_noan', yshift + interp1(waves.tt, rf, t_pred_xi_noan, 'cubic'),...
        40, 'red', 'filled') % If using true parameters and isotropic stack
%     hnd_t_xi = scatter(...
%         t_pred_xi_best', yshift + interp1(waves.tt, rf, t_pred_xi_best, 'cubic'),...
%         100, '+blue') % If using true parameters and anisotropic stack
%     hnd_t_00 = scatter(...
%         t_pred_xi_noan', yshift + interp1(waves.tt, rf, t_pred_xi_noan, 'cubic'),...
%         100, '+red') % If using true parameters and isotropic stack
    hnd_rf = plot(waves.tt, yshift+rf, 'k', 'linewidth', 1.5);

    if ixi == nxi; 
        xilabel = '\xi = '; 
    else; 
        xilabel = "      "; 
    end

    text(1, yshift + yshift_const * .5, sprintf('%s%1.2f', xilabel, xi_a(ixi) ) )

end

% scatter(0, -2*yshift_const, 0.00001); 
ylim([-2*yshift_const, yshift_const * (nxi+2.5)])
lgd = legend([hnd_rf, hnd_t_00, hnd_t_xi], ...
    'Receiver function', 't ignore \xi', 't with \xi'); 
set(lgd, 'Orientation', 'horizontal', 'Location', 'south'); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'rftiming', 'pdf'), 'ContentType', 'vector'); 

%%
figure(203); clf; hold on; 
set(gcf, 'pos', [-826 509 291 168]); 
box on; 
grid on; 
set(gca,'LineWidth', 1.5); 
xlabel('\xi true')

yyaxis left; 
ylabel('H error (km)'); 
plot(xi_a, herr, 'o')
plot(xi_a, herr, '-')

yyaxis right; 
ylabel('\kappa error')
plot(xi_a, kerr, 'o'); 
plot(xi_a, kerr, '-'); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'rf_error', 'pdf'), 'ContentType', 'vector'); 

%% Try to plot objective function, or something like that. 

figure(401); clf; hold on; box on; grid on; set(gca,'LineWidth', 1.5); 



%%
zylim = [0, 60]; 
model_cor = each_model{i_xi_true}; %model we did the correction to 
LW = 1.5; 

figure(204); clf; hold on; 
set(gcf, 'pos', [-800 377 317 240]); 
tiledlayout(1,2,'TileSpacing','compact'); 
% sgtitle('Model'); 

nexttile(); hold on; set(gca, 'LineWidth', LW, 'YDir', 'reverse'); box on; ylim(zylim); 
xlabel('Velocity (km/s)'); 
grid on; 
plot(model_cor.VS, model_cor.z, 'DisplayName', 'VS', 'LineWidth', LW); 
plot(model_cor.VP, model_cor.z, 'DisplayName', 'Vp', 'LineWidth', LW); 
legend(); 
xlim([min(model_cor.VS-0.75), max(model_cor.VP+0.75)])

ylabel('Depth (km)'); 

nexttile(); hold on; set(gca, 'LineWidth', LW, 'YDir', 'reverse'); box on; ylim(zylim); 
xlabel('% Anisotropy'); 
grid on; 
set(gca, 'yticklabel', []); 
plot(   model_cor.Sanis, model_cor.z, 'DisplayName', '+ \xi', ...
    'LineWidth', LW, 'LineStyle','-'); 
plot( - model_cor.Panis, model_cor.z, 'DisplayName', '- \phi', ...
    'LineWidth', LW*1.5, 'LineStyle','--'); 
% xticks = string(get(gca, 'XTicklabel')) ; 
% xticks(xticks~="0") = ""; 
% set(gca, 'xticklabels', xticks); 
xlim([-18, 0 + 3]); 
legend(); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'synth_model', 'pdf'), 'ContentType', 'vector'); 
