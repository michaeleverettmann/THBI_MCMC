run_hk_test_setup_bs; % This does some obnoxious setup. % Run it then comment it out if you want to save some time. 

%% HK tests, analysis, starts here. 
xi_a = [0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15]'; 
xi_true = 0.85; 
if ~ any(xi_a == xi_true); error('Pick xi true that is in xi_a'); end
i_xi_true = find(xi_a == xi_true); 
nxi = length(xi_a); 

ztrue_a = [45  ]; 
ktrue_a = [1.75];

nzi = length(ztrue_a); 
nki = length(ktrue_a); 


% par.datprocess.kNum = 101; 
% par.datprocess.hNum = 100; 
par.datprocess.kNum = 201; 
par.datprocess.hNum = 200; 
warning('Less hk resolution right now')

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
ttl12 = title('Synthetic Model Example', 'fontweight', 'normal'); 
set(ttl12, 'Position', [1.30 1.0190 0.5000], 'Units', 'normalized') % Manually adjust title position to cover ax1 and ax2

axes(ax2); 
xlabel('\xi'); 
ylabel('Depth (km)'); % yticklabels([]); 

axes(ax3); 
xlabel('Time (s)'); 
yticklabels([]); 
ttl3 = title('Phase Timing', 'FontWeight','normal'); 

axes(ax4); 
xlabel('\kappa'); 
ylabel('H (km)'); 
title('HK stack: ignoring \xi', 'fontweight', 'normal'); 

axes(ax5); 
xlabel('\kappa'); 
ylabel('H (km)'); 
% yticklabels([]); 
title('HK stack: \xi corrected', 'fontweight', 'normal'); 


%%%%%%%% Synth model plots
zylim = [0, 60]; 
model_cor = each_model{i_xi_true}; %model we did the correction to 
LW = 1.5; 

axes(ax1); hold on; set(gca, 'YDir', 'reverse'); 
ylim(zylim); 
plot(model_cor.VS, model_cor.z, 'DisplayName', 'VS', 'LineWidth', LW); 
plot(model_cor.VP, model_cor.z, 'DisplayName', 'Vp', 'LineWidth', LW); 
legend('Location','northeast'); 
xlim([min(model_cor.VS-0.75), max(model_cor.VP+0.75)])

axes(ax2); 
set(gca, 'YDir', 'reverse'); 
ylim(zylim); 
axes(ax2); ylim(zylim); 

plot(   model_cor.Sanis/100+1, model_cor.z, 'DisplayName', '+ \xi', ...
    'LineWidth', LW, 'LineStyle','-'); 
plot( - model_cor.Panis/100+1, model_cor.z, 'DisplayName', '- \phi', ...
    'LineWidth', LW*1.5, 'LineStyle','--'); 
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

    text(0.75, yshift + yshift_const * .5, sprintf('%s%1.2f', xilabel, xi_a(ixi) ) )

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

hnd_true = scatter(ktrue, ztrue, 200, [252, 3, 252]./256, 'pentagram', 'filled', ...
    'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
    'DisplayName', 'True'); 


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

%%%%%%%% HK stack with anis
axes(ax5); 
set(gca,'ydir', 'reverse');%, 'LineWidth', 1.5);

ylim(plt_ylim); 
xlim(plt_xlim); 

hnd_true = scatter(ktrue, ztrue, 200, [252, 3, 252]./256, 'pentagram', 'filled', ...
    'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
    'DisplayName', 'True'); 

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

legend([hnd_xistart, hnd_xiend, hnd_xi1, hnd_true], 'Location', 'best'); 


exportgraphics(gcf, fhand_figname(ztrue, ktrue, 'rf_synth_paper', 'pdf'), 'ContentType', 'vector');  
%%% End main figure! 


%% HK stack error figure
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

%ztrue, ktrue, xi_true

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
set(gcf, 'pos', [2235 809 251 366]); 
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







