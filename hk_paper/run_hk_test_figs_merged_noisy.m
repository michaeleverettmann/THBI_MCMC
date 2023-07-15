
run_hk_test_setup_bs; % This does some obnoxious setup. % Run it then comment it out if you want to save some time. 
par.inv.datatypes = {'HKstack_P'}

%% HK tests, analysis, starts here. 
xi_a = [0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15]'; 
xi_true = 0.85; 
if ~ any(xi_a == xi_true); error('Pick xi true that is in xi_a'); end
i_xi_true = find(xi_a == xi_true); 
nxi = length(xi_a); 

ztrue_a = [45  ]; 
ktrue_a = [1.75];

noise_amt = 0; % 4*0.012; 
par.synth.noise_sigma_BW_Ps = noise_amt; 
par.synth.noise_sigma_RF_Ps = noise_amt; 
rng(1); % Set random number seed in case we want to add noise. 

name_modifier = ""; 
if noise_amt > 0; 
    name_modifier = name_modifier + "_noisy"; 
end

nzi = length(ztrue_a); 
nki = length(ktrue_a); 


% par.datprocess.kNum = 101; 
% par.datprocess.hNum = 100; 
par.datprocess.HKappa.kNum = 501; 
par.datprocess.HKappa.hNum = 500; 

% Exi_all = zeros(nzi, nki, nxi, ) % Nope, not worth it

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
    E00 = trudata.HKstack_P.HKstack_P_noan.Esum; 
    waves = trudata.HKstack_P.waves; 
    t_predxi = trudata.HKstack_P.t_pred; 
    t_pred00 = trudata.HKstack_P.HKstack_P_noan.t_pred; 

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

    fprintf('\nXi = %1.3f. \n    H = %1.2f, K = %1.3f.\n    dH = %1.3f. percent dH = %1.3f.\n    dK = %1.4f. Percent dK = %1.4f\n',...
        xitruei, hmax_noan, kmax_noan, hmax_noan-ztrue, (hmax_noan-ztrue)/ztrue*100, ...
        kmax_noan-ktrue, (kmax_noan-ktrue)/ktrue*100) % brb2023.03.31 Working on this. 

    %%% How much error in the non-anisotropic HK stack? 
%     interp2(H, K, E00, hmax, kmax)
    %%%
end


% % % %%
% % % % Pleft = 0.05; 
% % % % Pbot = 0.05; 
% % % % P
% % % % Ptop = 1-Pbot;
% % % % Pright = 1-Pleft; 
% % % 
% % % Sbuf = 0.08
% % % 
% % % MWIDTH = 0.12; 
% % % 
% % % figmain = figure(500); clf; %hold on; 
% % % ax1 = axes('Position', [Sbuf               , 0.5+Sbuf, MWIDTH           , 0.5-Sbuf*2]); 
% % % ax2 = axes('Position', [Sbuf*2+MWIDTH      , 0.5+Sbuf, MWIDTH           , 0.5-Sbuf*2]); 
% % % ax3 = axes('Position', [MWIDTH * 2 + Sbuf*3, 0.5+Sbuf, 1-2*MWIDTH-4*Sbuf, 0.5-Sbuf*2]); 
% % % 
% % % ax4 = axes('Position', [Sbuf               , Sbuf    , 0.5-Sbuf*3/2       , 0.5-Sbuf*2]);
% % % ax5 = axes('Position', [0.5+Sbuf*.5        , Sbuf    , 0.5-Sbuf*3/2       , 0.5-Sbuf*2]);
% % % 
% % % each_ax = [ax1, ax2, ax3, ax4, ax5]; 
% % % 
% % % for ieach_ax = 1:length(each_ax); 
% % %     axes(each_ax(ieach_ax)); 
% % %     set(gca, 'LineWidth', 1.5); 
% % %     box on; 
% % %     grid on; 
% % %     xlabel('Test'); 
% % %     ylabel('Test'); 
% % %     title('Test', 'Fontweight', 'norma'); 
% % % end

%% Set figure layout
c_xilow = [179, 54, 214]./255; % Low values of xi
c_xihi = [201, 135, 2]./255; 
c_mod_true = 'y'; % Yellow, star? 
c_t_with_xi = [3, 15, 252]./255; ;% Anisotropic 
c_t_no_xi = [230, 2, 14]./255; 
c_t_iso = [10, 247, 30]./255; 
cnt_mod = 1/1.5; % Multiply colors by this when doing contours. Make them a bit darker. 
xlbl = 0; 
ylbl = 1.05; 


figmain = figure(501); clf; 
set(figmain, 'pos', [2284 706 735 561]); 
nxt = 10; 
nyt = 10;
xy_to_t = @(x,y)(y-1)*nyt+x; 
tiledlayout(nxt, nyt, 'TileSpacing','compact');
ax1 = nexttile(xy_to_t(1,1), [5,2]); hold on; 
ax2 = nexttile(xy_to_t(3,1), [5,2]); hold on; 
ax3 = nexttile(xy_to_t(5,1), [5,6]); hold on; 
ax4 = nexttile(xy_to_t(1,6), [5,5]); hold on; 
ax5 = nexttile(xy_to_t(6,6), [5,5]); hold on; 

each_ax = [ax1, ax2, ax3, ax4, ax5]; 
for ieach_ax = 1:length(each_ax); 
    axes(each_ax(ieach_ax)); 
    set(gca, 'LineWidth', 1.5); box on; grid on; 
end

axes(ax1); 
xlabel('Velocity (km/s)'); 
ylabel('Depth (km)'); 
% ttl12 = title('Synthetic Model Example', 'fontweight', 'normal'); 
% set(ttl12, 'Position', [1.30 1.0190 0.5000], 'Units', 'normalized') % Manually adjust title position to cover ax1 and ax2
title('Velocity', 'fontweight', 'normal'); 
text(xlbl, ylbl, '(a)', 'Units','normalized'); 

axes(ax2); 
xlabel('Anisotropy'); 
ylabel('Depth (km)'); % yticklabels([]); 
title('Anisotropy', 'fontweight', 'normal');  
text(xlbl, ylbl, '(b)', 'Units','normalized'); 


axes(ax3); 
xlabel('Time (s)'); 
yticklabels([]); 
ttl3 = title('Synthetic Receiver Functions', 'FontWeight','normal'); 
text(xlbl, ylbl, '(c)', 'Units','normalized'); 

axes(ax4); 
xlabel('\kappa'); 
ylabel('H (km)'); 
title('H\kappa stack: ignoring \xi', 'fontweight', 'normal'); 
text(xlbl, ylbl, '(d)', 'Units','normalized'); 

axes(ax5); 
xlabel('\kappa'); 
ylabel('H (km)'); 
% yticklabels([]); 
title('H\kappa stack: \xi corrected', 'fontweight', 'normal'); 
text(xlbl, ylbl, '(e)', 'Units','normalized'); 


% %%%%%%%% Synth model plots
zylim = [0, 60]; 
model_cor = each_model{i_xi_true}; %model we did the correction to 
LW = 1.5; 

axes(ax1); hold on; set(gca, 'YDir', 'reverse'); 
ylim(zylim); 
plot(model_cor.VS, model_cor.z, 'k', 'DisplayName', 'Vs', 'LineWidth', LW); 
plot(model_cor.VP, model_cor.z, 'k', 'DisplayName', 'Vp', 'LineWidth', LW,...
    'LineStyle','--'); 
legend('Location','northeast'); 
xlim([min(model_cor.VS-0.75), max(model_cor.VP+0.75)])

axes(ax2); 
set(gca, 'YDir', 'reverse'); 
ylim(zylim); 
axes(ax2); ylim(zylim); 

plot(   model_cor.Sanis/100+1, model_cor.z, 'k', 'DisplayName', '\xi', ...
    'LineWidth', LW, 'LineStyle','-'); 
plot( - model_cor.Panis/100+1, model_cor.z, 'DisplayName', '1/\phi', ...
    'LineWidth', LW*1.5, 'LineStyle','--', 'Color', [4, 184, 139]./255); 
xlim([0.8, 1.05]); 
legend('Location','northeast'); 

%%%%%%%% Receiver function waveform plots
axes(ax3); 
xlim([-3, 30])
yshift_const = 0.075; 

for ixi = 1:nxi
    yshift = ixi * yshift_const; 
    rf = rf_all{ixi}; 

    t_pred_xi_best = t_pred_xi_best_all(ixi,:)'; 
    t_pred_xi_noan = t_pred_xi_noan_all(ixi,:)'; 

    hnd_t_xi = scatter(...
        t_pred_xi_best', yshift + interp1(waves.tt, rf, t_pred_xi_best, 'cubic'),...
        40, 'filled', 'MarkerFaceColor', c_t_with_xi); % If using true parameters and anisotropic stack
    hnd_t_00 = scatter(...
        t_pred_xi_noan', yshift + interp1(waves.tt, rf, t_pred_xi_noan, 'cubic'),...
        40, 'filled', 'MarkerFaceColor', c_t_no_xi); % If using true parameters and isotropic stack
    hnd_rf = plot(waves.tt, yshift+rf, 'k', 'linewidth', 1.5);

    if ixi == nxi; 
        xilabel = '\xi = '; 
    else; 
        xilabel = "      "; 
    end

    text(0.75, yshift + yshift_const * .5, sprintf('%s%1.2f', xilabel, xi_a(ixi) ) );

end

ylim([-2*yshift_const, yshift_const * (nxi+2.5)])
lgd = legend([hnd_rf, hnd_t_00, hnd_t_xi], ...
    'Receiver function', 't ignore \xi', 't with \xi'); 
set(lgd, 'Orientation', 'horizontal', 'Location', 'south'); 

%%%%%%%% HK stack, no anis

% First, a couple things that will apply to corrected, and uncorrected plots. 
fhand_norm = @(inval)inval ./ max(max(inval)); % Return normalized inval 
lvl_cnt = [0.7, 0.95]; 
LW = 1; 

axes(ax4); 
set(gca,'ydir', 'reverse');

plt_ylim = [40, 50]; 
plt_xlim = [1.6, 1.9]; 
ylim(plt_ylim); 
xlim(plt_xlim); 

LW_scale = 1.75; 
[~,hnd_xiend  ] = contour(K, H, fhand_norm(E00_all{xi_a==1}'),...
    lvl_cnt, 'k', 'LineWidth', LW*LW_scale, ... 'color', c_t_with_xi.*cnt_mod, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(xi_a==1) ) ); 
[~,hnd_xistart] = contour(K, H, fhand_norm(E00_all{1      }'),...
    lvl_cnt, 'LineWidth', LW*LW_scale, 'color', c_xilow,...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(1      ) ),...
    'LineStyle','-'); 
[~,hnd_xi1    ] = contour(K, H, fhand_norm(E00_all{end    }'),...
    lvl_cnt, 'LineWidth', LW*LW_scale, 'color', c_xihi,...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(end    ) ),...
    'LineStyle','-'); 

hnd_true = scatter(ktrue, ztrue, 200, c_mod_true, 'pentagram', 'filled', ...
    'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
    'DisplayName', 'True'); 

%%%%%%%% HK stack with anis
axes(ax5); 
set(gca,'ydir', 'reverse');%, 'LineWidth', 1.5);

ylim(plt_ylim); 
xlim(plt_xlim); 

[~,hnd_xiend  ] = contour(K, H, fhand_norm(Exi_all{xi_a==1}'),...
    lvl_cnt, 'k', 'LineWidth', LW*LW_scale, ...'color', c_t_with_xi.*cnt_mod, ...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(xi_a==1) ) ); 
[~,hnd_xistart] = contour(K, H, fhand_norm(Exi_all{1      }'),...
    lvl_cnt, 'LineWidth', LW*LW_scale, 'color', c_xilow,...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(1      ) ),...
    'LineStyle','-'); 
[~,hnd_xi1    ] = contour(K, H, fhand_norm(Exi_all{end    }'),...
    lvl_cnt, 'LineWidth', LW*LW_scale, 'color', c_xihi,...
    'DisplayName', sprintf('\\xi = %1.2f', xi_a(end    ) ),...
    'LineStyle','-'); 


hnd_true = scatter(ktrue, ztrue, 200, c_mod_true, 'pentagram', 'filled', ...
    'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
    'DisplayName', 'True'); 

legend([hnd_xistart, hnd_xiend, hnd_xi1, hnd_true], 'Location', 'best'); 


exportgraphics(gcf, fhand_figname(ztrue, ktrue, "rf_synth_paper"+name_modifier, 'pdf'), 'ContentType', 'vector');  
%%% End main figure! 

%% Some calculations on change in HK stacks with anisotropy
% Relies on external function, get_hk_best.m 
% [~, ~, H_bst_low, K_bst_low] = get_hk_best(E00_all, 1, H, K); 
% [~, ~, H_bst_hi , K_bst_hi ] = get_hk_best(E00_all, size(E00_all, 1), H, K); 
% [~, ~, H_bst_00 , K_bst_00 ] = get_hk_best(E00_all, xi_a == 1, H, K); 
% 
% (H_bst_hi - H_bst_00)
% [~, ~, H_bst_00 , K_bst_00 ] = get_hk_best(E00_all, xi_a == 1, H, K); 
% fprintf('\n\n')
% for ixi = 1:length(xi_a); 
%     [~, ~, H_bst_xi , K_bst_xi ] = get_hk_best(E00_all,ixi, H, K); 
%     dH = -(H_bst_xi - H_bst_00); 
%     dH_prct = dH/H_bst_00*100; 
%     dK = -(K_bst_xi - K_bst_00); 
%     dK_prct = dK/K_bst_00*100; 
%     fprintf('\nXi = %1.3f. \n    H = %1.2f, K = %1.3f.\n    dH = %1.3f. percent dH = %1.3f.\n    dK = %1.4f. Percent dK = %1.4f\n',...
%         xi_a(ixi), H_bst_xi, K_bst_xi, dH, dH_prct, dK, dK_prct)
% end


%% HK stack error figure
LW = 1.5; 
figure(203); clf; hold on; 
set(gcf, 'pos', [-826 509 291 168]); 
box on; 
grid on; 
set(gca,'LineWidth', 1.5); 
xlabel('\xi true')

yyaxis left; 
ylabel('H offset (km)'); 
scatter(xi_a, herr, 50, 'filled'); 
plot(xi_a, herr, '-', 'LineWidth', LW); 

yyaxis right; 
ylabel('\kappa offset')
scatter(xi_a, kerr, 50, 'filled'); 
plot(xi_a, kerr, '-', 'LineWidth', LW); 

exportgraphics(gcf, fhand_figname(ztrue, ktrue, "rf_error"+name_modifier, 'pdf'), 'ContentType', 'vector'); 

%% Try to plot objective function, or something like that. 

%ztrue, ktrue, xi_true
E_intersect = 0.9;  % Want to know which values intersect this high of E

% For this, need to make E that is nxi x nh x nk, (or something like that).
Eob = zeros(nxi, par.datprocess.HKappa.kNum, par.datprocess.HKappa.hNum); 
for ixi = 1:nxi; 
    Eob(ixi,:,:) = Exi_all{ixi}; 
end

Eob = Eob ./ max(Eob, [], 'all'); % Normalize it. 

ind_closest = @(H, ztrue)find(abs(H - ztrue) == min(abs(H - ztrue))); % Index closest to H in ztrue. 

ztruei = ind_closest(H, ztrue); 
ktruei = ind_closest(K, ktrue); 
xitruei= ind_closest(xi_a, xi_true); 

Ebest = Eob(xitruei, ktruei, ztruei); 

Etmp = Eob(:, ktruei, :); 

flatar = @(M)reshape(M,[],1); % flat array

% ylim_man = [-0.01, 0.11]; 
ylim_man = [-.1, 1.1]; 
% ylim_man = [0.75, 0.76]; 

figure(401); clf; hold on; box on; grid on; 
set(gcf, 'pos', [2235 809 251 366]); 
tiledlayout(3,1,'TileSpacing','compact'); 
LW = 1.25; 

nexttile(), hold on, ylim(ylim_man), box on, grid on; set(gca,'LineWidth', LW); 
xlabel('H (km)'); xlim([25, 60]); 
toplt = flatar(Eob(xitruei, ktruei, :)); 
plot( H, toplt ,...
    'k', 'LineWidth', LW); 
sct = scatter(ztrue, Ebest, 150, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k'); 
text(0.03, 0.1, sprintf('$\\kappa=%1.2f$ km, $\\xi = %1.2f$', ktrue, xi_true), ...
    'Units', 'normalized', 'VerticalAlignment', 'bottom',...
    'Interpreter','latex'); 
intersecting_x = funct_find_intersection(H, toplt, E_intersect) ; 
fprintf('H crosses %1.2f at [%s]\n\n', E_intersect, sprintf('%1.3f, ', intersecting_x))


nexttile(), hold on, ylim(ylim_man), box on, grid on; set(gca,'LineWidth', LW); 
xlabel('\kappa'); % xlim([25, 60]); 
ylabel('E'); 
toplt = Eob(xitruei, :, ztruei); 
plot( K, toplt ,...
    'k', 'LineWidth', LW); 
sct = scatter(ktrue, Ebest, 150, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k'); 
text(0.03, 0.1, sprintf('$H=%1.1f$ km, $\\xi = %1.2f$', ztrue, xi_true), ...
    'Units', 'normalized', 'VerticalAlignment', 'bottom',...
    'Interpreter','latex'); 
intersecting_x = funct_find_intersection(K, toplt, E_intersect) ; 
fprintf('K crosses %1.2f at [%s]\n\n', E_intersect, sprintf('%1.3f, ', intersecting_x))

nexttile(), hold on, ylim(ylim_man), box on, grid on; set(gca,'LineWidth', LW); 
xlabel('\xi'); 
toplt = Eob(:,ktruei, ztruei); 
plot( xi_a, toplt ,...
    'k', 'LineWidth', LW); 
sct = scatter(xi_true, Ebest, 150, 'yellow', 'filled', 'pentagram', 'MarkerEdgeColor', 'k'); 
% xlim([0.8, 1.2]); 
text(0.03, 0.1, sprintf('$H=%1.1f$ km, $\\kappa = %1.2f$', ztrue, ktrue), ...
    'Units', 'normalized', 'VerticalAlignment', 'bottom',...
    'Interpreter','latex'); 
intersecting_x = funct_find_intersection(xi_a, toplt, E_intersect) ; 
fprintf('xi crosses %1.2f at [%s]\n\n', E_intersect, sprintf('%1.3f, ', intersecting_x))

exportgraphics(gcf, fhand_figname(ztrue, ktrue, "HK_objective_function"+name_modifier, 'pdf'), 'ContentType', 'vector'); 



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







