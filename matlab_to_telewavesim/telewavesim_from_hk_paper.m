restoredefaultpath(); 
run_hk_test_setup_bs; % This does some obnoxious setup. % Run it then comment it out if you want to save some time. 
addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/hk_paper'); % For a2_LOAD_DATA_hk_test
% addpath('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_telewavesim'); 
pyenv('Version', '~/opt/anaconda3/envs/tws/bin/python', ... % Use anaconda environment where telewavesim is installed. This affects the entire Matlab session. TODO define this path in somewhere more obvious. 
    'ExecutionMode','OutOfProcess'); % ERROR ALERT Could not import numpy if using an anconda environment. Matlab would simply crash. However, setting executionMode=OutOfProcess fixed that for me. https://www.mathworks.com/matlabcentral/answers/502458-why-can-py-numpy-array-not-be-resolved
if count(py.sys.path, '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_telewavesim') == 0
    insert(py.sys.path, int32(0), '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_telewavesim');
end


%% HK tests, analysis, starts here. 
xi_a = [0.85, 1]'; 
xi_true = 0.85; 
if ~ any(xi_a == xi_true); error('Pick xi true that is in xi_a'); end
i_xi_true = find(xi_a == xi_true); 
nxi = length(xi_a); 

ztrue_a = [45  ]; 
ktrue_a = [1.75];

nzi = length(ztrue_a); 
nki = length(ktrue_a); 


par.datprocess.kNum = 201; 
par.datprocess.hNum = 200; 

Exi_all = cell([length(xi_a), 1]); 
E00_all = cell([length(xi_a), 1]); 

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
    xitruei = xi_a(ixi); 

    [trudata,par] = a2_LOAD_DATA_hk_test(par, 'nwk', nwk, 'sta', sta, ...
        'xi_crust', xitruei );

    each_model{ixi, 1} = TRUEmodel; 
    ztrue = TRUEmodel.zmoh; 
    ktrue = TRUEmodel.vpvs; 
    
    H = trudata.HKstack_P.H;
    K = trudata.HKstack_P.K; 
    Exi = trudata.HKstack_P.Esum; 
    E00 = trudata.HKstack_P.HKstack_P_noan.Esum; 
    waves = trudata.HKstack_P.waves; 
    t_predxi = trudata.HKstack_P.t_pred; 
    t_pred00 = trudata.HKstack_P.HKstack_P_noan.t_pred; 

    [ikmax, ihmax] = find(E00 == max(E00 ,[], 'all')); 
    warning('I think I had k and h indicies backwards in the above line')
    kmax_noan = K(ikmax); 
    hmax_noan = H(ihmax); 
    herr(ixi) = hmax_noan - ztrue; 
    kerr(ixi) = kmax_noan - ktrue; 

    [ikmax, ihmax] = find(Exi == max(Exi ,[], 'all')); 
    kmax_best = K(ikmax); 
    hmax_best = H(ihmax); 
    
    t_pred_xi_best  = zeros(1, 3); % Get the times from anisotropic stack at true parameters
    t_pred_xi_noan  = zeros(1, 3); % Get the times from ISOTROPIC stack at true parameters. 
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
    hnd_rf = plot(waves.tt, yshift+rf, 'k', 'linewidth', 1.5);

    if ixi == nxi; 
        xilabel = '\xi = '; 
    else; 
        xilabel = "      "; 
    end

    text(1, yshift + yshift_const * .5, sprintf('%s%1.2f', xilabel, xi_a(ixi) ) )

end

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
% For this, need to make E that is nxi x nh x nk, (or something like that).
Eob = zeros(nxi, par.datprocess.kNum, par.datprocess.hNum); 
for ixi = 1:nxi; 
    Eob(ixi,:,:) = Exi_all{ixi}; 
end

ind_closest = @(H, ztrue)find(abs(H - ztrue) == min(abs(H - ztrue))); % Index closest to H in ztrue. 

ztruei = ind_closest(H, ztrue); 
ktruei = ind_closest(K, ktrue); 
xitruei= ind_closest(xi_a, xi_true); 

Ebest = Eob(xitruei, ktruei, ztruei); 

Etmp = Eob(:, ktruei, :); 

flatar = @(M)reshape(M,[],1); % flat array

ylim_man = [-0.01, 0.11]; 
figure(401); clf; hold on; box on; grid on; 
set(gcf, 'pos', [2235 809 407 359]); 
tiledlayout(3,1,'TileSpacing','compact'); 
LW = 1.25; 

nexttile(), hold on, ylim(ylim_man), box on, grid on; set(gca,'LineWidth', LW); 
xlabel('H (km)'); xlim([25, 60]); 
plot( H, flatar(Eob(xitruei, ktruei, :)) ,...
    'k', 'LineWidth', LW); 
sct = scatter(ztrue, Ebest, 150, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k'); 
text(0.03, 0.1, sprintf('$\\kappa=%1.2f$ km, $\\xi = %1.2f$', ktrue, xi_true), ...
    'Units', 'normalized', 'VerticalAlignment', 'bottom',...
    'Interpreter','latex'); 

nexttile(), hold on, ylim(ylim_man), box on, grid on; set(gca,'LineWidth', LW); 
xlabel('\kappa'); % xlim([25, 60]); 
ylabel('E')
plot( K, Eob(xitruei, :, ztruei) ,...
    'k', 'LineWidth', LW); 
sct = scatter(ktrue, Ebest, 150, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k'); 
text(0.03, 0.1, sprintf('$H=%1.1f$ km, $\\xi = %1.2f$', ztrue, xi_true), ...
    'Units', 'normalized', 'VerticalAlignment', 'bottom',...
    'Interpreter','latex'); 

nexttile(), hold on, ylim(ylim_man), box on, grid on; set(gca,'LineWidth', LW); 
xlabel('\xi'); 
plot( xi_a, Eob(:,ktruei, ztruei) ,...
    'k', 'LineWidth', LW); 
sct = scatter(xi_true, Ebest, 150, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k'); 
% xlim([0.8, 1.2]); 
text(0.03, 0.1, sprintf('$H=%1.1f$ km, $\\kappa = %1.2f$', ztrue, ktrue), ...
    'Units', 'normalized', 'VerticalAlignment', 'bottom',...
    'Interpreter','latex'); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'HK_objective_function', 'pdf'), 'ContentType', 'vector'); 


%%
zylim = [0, 60]; 
model_cor = each_model{i_xi_true}; %model we did the correction to 
LW = 1.5; 

figure(204); clf; hold on; 
set(gcf, 'pos', [-800 377 317 240]); 
tiledlayout(1,2,'TileSpacing','compact'); 

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
xlim([-18, 0 + 3]); 
legend(); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'synth_model', 'pdf'), 'ContentType', 'vector'); 
